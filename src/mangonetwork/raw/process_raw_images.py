#!/usr/bin/env python3
"""MANGO Level 1 Image Processing"""

##########################################################################
#
#   MANGO Raw Image Post Processor
#
#   Create level-1 data files
#
#   2022-xx-xx  Leslie Lamarche and Asti Bhatt
#               Initial implementation
#
#   2023-03-04  Todd Valentic
#               Refactor code to support parallel processing
#               Confirm to PEP8 formatting
#
##########################################################################

import argparse
import configparser
import logging
import multiprocessing
import pathlib
import os
import sys

import h5py
import numpy as np
from scipy.interpolate import griddata

from . import imageops

if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources

RE = 6371.0  # Earth radius (m)

# -------------------------------------------------------------------------
# Image Processor
# -------------------------------------------------------------------------


class ImageProcessor:
    """Process a single image"""

    def __init__(self, config):
        self.config = config

        ha = self.config.getfloat("PROCESSING", "ALTITUDE")
        self.REha = RE + ha

        # Results

        self.image = None
        self.metadata = None
        self.azimuth = None
        self.elevation = None
        self.latitude = None
        self.longitude = None
        self.image_mask = None
        self.new_x_grid = None
        self.new_y_grid = None
        self.trans_x_grid = None
        self.trans_y_grid = None

    def run(self, filename):
        """Run processing algorithm"""

        logging.debug(filename)

        raw_image = h5py.File(filename)["image"]

        self.metadata = self.get_metadata(raw_image)
        self.metadata["filename"] = filename

        self.create_transform_grids(raw_image)
        self.create_position_arrays(raw_image)

        self.image = self.process(raw_image)

        logging.debug("Processing finished")

    def get_metadata(self, image):
        """Extract metadata"""

        metadata = {}

        start_time = image.attrs["start_time"]
        exposure_time = float(image.attrs["exposure_time"])

        metadata["start_time"] = start_time
        metadata["end_time"] = start_time + exposure_time
        metadata["ccd_temp"] = image.attrs["ccd_temp"]
        metadata["code"] = image.attrs["station"]
        metadata["label"] = image.attrs["label"]
        metadata["site_lat"] = image.attrs["latitude"]
        metadata["site_lon"] = image.attrs["longitude"]

        return metadata

    def create_transform_grids(self, image):
        """Create transformation grids"""

        # pylint: disable=too-many-locals

        i_max = image.attrs["width"]
        j_max = image.attrs["height"]

        new_i_max = self.config.getint("PROCESSING", "NEWIMAX")
        new_j_max = self.config.getint("PROCESSING", "NEWJMAX")

        x0 = self.config.getfloat("CALIBRATION_PARAMS", "X0")
        y0 = self.config.getfloat("CALIBRATION_PARAMS", "Y0")
        rl = self.config.getfloat("CALIBRATION_PARAMS", "RL")
        theta = self.config.getfloat("CALIBRATION_PARAMS", "THETA")

        A = self.config.getfloat("CALIBRATION_PARAMS", "A")
        B = self.config.getfloat("CALIBRATION_PARAMS", "B")
        C = self.config.getfloat("CALIBRATION_PARAMS", "C")
        D = self.config.getfloat("CALIBRATION_PARAMS", "D")

        elev_cutoff = self.config.getfloat("PROCESSING", "ELEVCUTOFF")

        # Create a grid and find the transformed coordinates of each grid point

        x_grid, y_grid = np.meshgrid(np.arange(i_max), np.arange(j_max))

        self.trans_x_grid, self.trans_y_grid = self.transform(
            x_grid, y_grid, x0, y0, rl, theta, A, B, C, D
        )

        # Find unwarped distance to elevation cutoff and create new grid

        d = self.unwarp(elev_cutoff * np.pi / 180.0)

        self.new_x_grid, self.new_y_grid = np.meshgrid(
            np.linspace(-d, d, new_i_max), np.linspace(-d, d, new_j_max)
        )

    # pylint: disable=too-many-arguments, too-many-locals

    def transform(self, x, y, x0, y0, rl, theta, A, B, C, D):
        """Coordinate transform"""

        x1 = (x - x0) / rl
        y1 = (y - y0) / rl

        t = theta * np.pi / 180.0
        x2 = np.cos(t) * x1 - np.sin(t) * y1
        y2 = np.sin(t) * x1 + np.cos(t) * y1

        r = np.sqrt(x2**2 + y2**2)
        lam = A + B * r + C * r**2 + D * r**3
        d = self.unwarp(lam)

        x3 = d * x2 / r
        y3 = d * y2 / r

        return x3, y3

    def unwarp(self, lam):
        """Convert elevation angle to unwarped distance"""

        b = np.arcsin(RE / self.REha * np.sin(lam + np.pi / 2))  # Law of Sines
        psi = np.pi / 2 - lam - b
        d = self.REha * psi

        return d

    def create_position_arrays(self, image):
        """Create position arrays"""

        site_lat = image.attrs["latitude"] * np.pi / 180
        site_lon = image.attrs["longitude"] * np.pi / 180

        elev_cutoff = self.config.getfloat("PROCESSING", "ELEVCUTOFF")

        psi = np.sqrt(self.new_x_grid**2 + self.new_y_grid**2) / self.REha

        # Law of Cosines

        c = np.sqrt(RE**2 + self.REha**2 - 2 * RE * self.REha * np.cos(psi))
        b = self.REha - RE * np.cos(psi)

        elv = np.pi / 2.0 - np.arccos(b / c) - psi
        azm = np.arctan2(self.new_x_grid, self.new_y_grid)
        self.image_mask = elv < elev_cutoff * np.pi / 180.0

        # Use Haversine equations to find lat/lon of each point

        lat = np.arcsin(
            np.sin(site_lat) * np.cos(psi)
            + np.cos(site_lat) * np.sin(psi) * np.cos(azm)
        )

        lon = site_lon + np.arctan2(
            np.sin(azm) * np.sin(psi) * np.cos(site_lat),
            np.cos(psi) - np.sin(site_lat) * np.sin(lat),
        )

        self.azimuth = azm * 180.0 / np.pi
        self.elevation = elv * 180.0 / np.pi
        self.latitude = lat * 180.0 / np.pi
        self.longitude = lon * 180.0 / np.pi


    def atmospheric_correction(self, image, vanrhijn=True, extinction=True):
        """Apply atmospheric correction"""

        # Atmospheric corrections are taken from Kubota et al., 2001
        # Kubota, M., Fukunishi, H. & Okano, S. Characteristics of medium- and
        #   large-scale TIDs over Japan derived from OI 630-nm nightglow observation.
        #   Earth Planet Sp 53, 741–751 (2001). https://doi.org/10.1186/BF03352402

        # calculate zenith angle

        za = np.pi / 2 - self.elevation * np.pi / 180.0

        if vanrhijn:

            # Kubota et al., 2001; eqn. 6

            vanrhijn_factor = np.sqrt(1.0 - np.sin(za) ** 2 * (RE / self.REha) ** 2)

        else:
            
            vanrhijn_factor = np.ones(za.shape)


        if extinction:

            # Kubota et al., 2001; eqn. 7,8

            a = 0.2
            F = 1.0 / (np.cos(za) + 0.15 * (93.885 - za * 180.0 / np.pi) ** (-1.253))
            extinction_factor = 10.0 ** (0.4 * a * F)

        else:

            extinction_factor = np.ones(za.shape)

        correction = vanrhijn_factor * extinction_factor

        return image * correction

    def process(self, raw_image):
        """Processing algorithm"""

        elev_cutoff = self.config.getfloat("PROCESSING", "ELEVCUTOFF")
        remove_background = self.config.getboolean("PROCESSING", "REMOVE_BACKGROUND")
        contrast = self.config.getfloat("PROCESSING", "CONTRAST", fallback=100)
        histequal = self.config.getboolean("PROCESSING", "EQUALIZATION")
        vanrhijn = self.config.getboolean("PROCESSING", "VANRHIJN")
        extinction = self.config.getboolean("PROCESSING", "EXTINCTION")
        uint8_out = self.config.getboolean("PROCESSING", "UINT8_OUT", fallback=False)

        cooked_image = np.array(raw_image)

        # Does it matter which of these operations is performed first?
        if remove_background:
            cooked_image = imageops.background_removal(cooked_image)

        if histequal:
            cooked_image = imageops.equalize(cooked_image, contrast)

        new_image = griddata(
            (self.trans_x_grid.flatten(), self.trans_y_grid.flatten()),
            cooked_image.flatten(),
            (self.new_x_grid, self.new_y_grid),
            fill_value=0,
        )

        # Apply atmopsheric correction

        new_image = self.atmospheric_correction(
            new_image, vanrhijn=vanrhijn, extinction=extinction
        )

        # Apply mask outside elevation cutoff

        new_image[self.elevation < elev_cutoff] = 0.0

        # Renormalize each image and convert to int

        if uint8_out:
           new_image = (new_image * 255 / np.nanmax(new_image)).astype("uint8")

        return new_image


# -------------------------------------------------------------------------
# Application methods
# -------------------------------------------------------------------------

# pylint: disable=too-many-statements, too-many-locals


def write_to_hdf5(output_file, config, results):
    """Save results to an HDF5 file"""

    site_name = config.get("SITE_INFO", "SITE_NAME")
    ha = config.getfloat("PROCESSING", "ALTITUDE")
    elevcutoff = config.getfloat("PROCESSING", "ELEVCUTOFF")
    backremove = config.getboolean("PROCESSING", "REMOVE_BACKGROUND")
    contrast = config.getfloat("PROCESSING", "CONTRAST", fallback=100)
    histequal = config.getboolean("PROCESSING", "EQUALIZATION")
    vanrhijn = config.getboolean("PROCESSING", "VANRHIJN")
    extinction = config.getboolean("PROCESSING", "EXTINCTION")

    start_time = [rec.metadata["start_time"] for rec in results]
    end_time = [rec.metadata["end_time"] for rec in results]
    ccd_temp = [rec.metadata["ccd_temp"] for rec in results]

    rec = results[0]

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(output_file, "w") as f:
        f.create_group("SiteInfo")
        f.create_group("ProcessingInfo")

        t = f.create_dataset(
            "UnixTime",
            data=np.array([start_time, end_time]),
            compression="gzip",
            compression_opts=1,
        )
        t.attrs["Description"] = "unix time stamp"
        t.attrs["Unit"] = "seconds"

        lat = f.create_dataset(
            "Latitude", data=rec.latitude, compression="gzip", compression_opts=1
        )
        lat.attrs["Description"] = "geodetic latitude of each pixel"
        lat.attrs["Size"] = "Ipixels x Jpixels"
        lat.attrs["Projection Altitude"] = ha
        lat.attrs["Unit"] = "degrees"

        lon = f.create_dataset(
            "Longitude", data=rec.longitude, compression="gzip", compression_opts=1
        )
        lon.attrs["Description"] = "geodetic longitude of each pixel"
        lon.attrs["Size"] = "Ipixels x Jpixels"
        lon.attrs["Projection Altitude"] = ha
        lon.attrs["Unit"] = "degrees"

        az = f.create_dataset(
            "Azimuth", data=rec.azimuth, compression="gzip", compression_opts=1
        )
        az.attrs["Description"] = "azimuth of each pixel"
        az.attrs["Size"] = "Ipixels x Jpixels"
        az.attrs["Unit"] = "degrees"

        el = f.create_dataset(
            "Elevation", data=rec.elevation, compression="gzip", compression_opts=1
        )
        el.attrs["Description"] = "elevation of each pixel"
        el.attrs["Size"] = "Ipixels x Jpixels"
        el.attrs["Unit"] = "degrees"

        pc = f.create_dataset(
            "PixelCoordinates",
            data=np.array([rec.new_x_grid, rec.new_y_grid]),
            compression="gzip",
            compression_opts=1,
        )
        pc.attrs[
            "Description"
        ] = "coordinates of each pixel in grid at the airglow altitude"
        pc.attrs["Size"] = "2 (X,Y) x Ipixels x Jpixels"
        pc.attrs["Projection Altitude"] = ha
        pc.attrs["Unit"] = "km"

        images = f.create_dataset(
            "ImageData",
            data=np.array([rec.image for rec in results]),
            compression="gzip",
            compression_opts=1,
        )
        images.attrs["Description"] = "pixel values for images"
        images.attrs["Site Abbreviation"] = rec.metadata["code"]
        images.attrs["Image label"] = rec.metadata["label"]

        mask = f.create_dataset(
            "Mask", data=rec.image_mask, compression="gzip", compression_opts=1
        )
        mask.attrs["Description"] = "mask of where ImageData array is corners ouside of camera FoV"

        ccd = f.create_dataset("CCDTemperature", data=ccd_temp)
        ccd.attrs["Description"] = "Temperature of CCD"
        ccd.attrs["Unit"] = "degrees C"

        name = f.create_dataset("SiteInfo/Name", data=site_name)
        name.attrs["Description"] = "site name"

        code = f.create_dataset("SiteInfo/Code", data=rec.metadata["code"])
        code.attrs["Description"] = "one letter site abbreviation/code"

        coord = f.create_dataset(
            "SiteInfo/Coordinates",
            data=[rec.metadata["site_lat"], rec.metadata["site_lon"]],
        )
        coord.attrs["Description"] = "geodetic coordinates of site; [lat, lon]"
        coord.attrs["Unit"] = "degrees"

        ec = f.create_dataset("ProcessingInfo/ElevationCutoff", data=elevcutoff)
        ec.attrs["Description"] = "elevation angle cutoff [deg]"

        cont = f.create_dataset("ProcessingInfo/Contrast", data=contrast)
        cont.attrs["Description"] = "contrast value used for histogram equalization"

        he = f.create_dataset("ProcessingInfo/HistogramEqualization", data=histequal)
        he.attrs["Description"] = "0 = no histogram equalization applied to image; 1 = histogram equalization applied to image"

        br = f.create_dataset("ProcessingInfo/BackgroundRemoval", data=backremove)
        br.attrs["Description"] = "0 = background has not been subtracted from image; 1 = background has been subtracted from image"

        vr = f.create_dataset("ProcessingInfo/VanRhijnCorrection", data=vanrhijn)
        vr.attrs["Description"] = "0 = Van Rhijn effect has not be removed from the image; 1 = Van Rhijn effect has been removed from the image"

        ef = f.create_dataset("ProcessingInfo/ExtinctionFactor", data=extinction)
        ef.attrs["Description"] = "0 = extinction factor atmospheric correction has not been applied to the image; 1 = extinction factor atmospheric correction has been applied to the image"


def worker(filename):
    """Parallel processing handleer for a single image"""

    # pylint: disable=global-variable-undefined

    processor = ImageProcessor(main_config)
    processor.run(filename)

    return processor


def worker_init(config):
    """Initialize parallel processing handler"""

    # pylint: disable=global-variable-undefined

    global main_config

    main_config = config


# -------------------------------------------------------------------------
# Main application
# -------------------------------------------------------------------------


def parse_args():
    """Command line options"""

    parser = argparse.ArgumentParser(description="Create a quick look movie")

    parser.add_argument(
        "-c", "--config", metavar="FILE", help="Alternate configuration file"
    )
    parser.add_argument(
        "-f",
        "--filelist",
        metavar="FILE",
        help="A file with a list of .hdf5 file names",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="mango.hdf5",
        help="Output filename (default is mango.hdf5)",
    )
    parser.add_argument(
        "-n", "--numproc", type=int, default=1, help="Number of parallel processes"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbose output (repeat for more detail)",
    )

    parser.add_argument("inputfiles", nargs="*")

    return parser.parse_args()


def find_inputfiles(args):
    """Find input filenames"""

    if args.filelist:
        with open(args.filelist, encoding="utf-8") as f:
            filenames = [line.strip() for line in f]
            # filter blank lines
            filenames = [line for line in filenames if line]
    else:
        filenames = args.inputfiles

    return [filename for filename in filenames if os.path.exists(filename)]


def find_config(filename):
    """Load configuration file from package data"""

    with h5py.File(filename) as h5:
        station = h5["image"].attrs["station"]
        instrument = h5["image"].attrs["instrument"]

    name = f"{station}-{instrument}.ini"

    logging.debug("Using package configuration file: %s", name)

    return resources.files("mangonetwork.raw.data").joinpath(name).read_text()


def main():
    """Main application"""

    args = parse_args()

    fmt = "[%(asctime)s] %(levelname)s %(message)s"

    if args.verbose > 0:
        logging.basicConfig(format=fmt, level=logging.DEBUG)
    else:
        logging.basicConfig(format=fmt, level=logging.INFO)

    inputs = find_inputfiles(args)

    if not inputs:
        logging.error("No input files found")
        sys.exit(1)

    logging.debug("Processing %d files", len(inputs))
    logging.debug("Number of processes: %d", args.numproc)

    if args.config:
        logging.debug("Alternate configuration file: %s", args.config)
        if not os.path.exists(args.config):
            logging.error("Config file not found")
            sys.exit(1)
        with open(args.config, encoding="utf-8") as f:
            contents = f.read()
    else:
        contents = find_config(inputs[0])

    config = configparser.ConfigParser()
    config.read_string(contents)

    with multiprocessing.Pool(
        processes=args.numproc, initializer=worker_init, initargs=(config,)
    ) as pool:
        results = pool.map(worker, inputs, chunksize=1)

    output_file = pathlib.Path(args.output)
    write_to_hdf5(output_file, config, results)

    sys.exit(0)


if __name__ == "__main__":
    main()
