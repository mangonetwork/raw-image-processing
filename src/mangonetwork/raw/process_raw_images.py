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
import platform
import io
import datetime as dt

import h5py
import numpy as np
from scipy.interpolate import griddata


import matplotlib.pyplot as plt
# skyfield stuff needed for moon flags
from skyfield import almanac
from skyfield.api import N, S, E, W, load, wgs84

if sys.version_info < (3, 9):
    import importlib_resources as resources
    import importlib_metadata as metadata
else:
    from importlib import resources
    from importlib import metadata

RE = 6371.0  # Earth radius (m)

# -------------------------------------------------------------------------
# Image Processor
# -------------------------------------------------------------------------


class ImageProcessor:
    """Process a set image to produce a nightly file"""

    def __init__(self, config):
        self.config = config

        ha = self.config.getfloat("PROCESSING", "ALTITUDE")
        self.REha = RE + ha
        self.elev_cutoff = self.config.getfloat("PROCESSING", "ELEVCUTOFF")
        self.remove_background = self.config.getboolean("PROCESSING", "REMOVE_BACKGROUND")

        self.x0 = self.config.getfloat("CALIBRATION_PARAMS", "X0")
        self.y0 = self.config.getfloat("CALIBRATION_PARAMS", "Y0")
        self.rl = self.config.getfloat("CALIBRATION_PARAMS", "RL")
        self.theta = self.config.getfloat("CALIBRATION_PARAMS", "THETA")

        self.A = self.config.getfloat("CALIBRATION_PARAMS", "A")
        self.B = self.config.getfloat("CALIBRATION_PARAMS", "B")
        self.C = self.config.getfloat("CALIBRATION_PARAMS", "C")
        self.D = self.config.getfloat("CALIBRATION_PARAMS", "D")


    def setup(self, filename):
        """Setup the processing algorithm by collecting the time independent metadata and calculating the position arrays that will be constant over the entire night"""

        logging.debug("Setup")

        logging.debug(f"Extracting grid and metadata from: {filename}")

        raw_image = h5py.File(filename)["image"]

        self.metadata = self.get_metadata(raw_image)

        self.create_transform_grids(raw_image)
        self.create_position_arrays()

        self.image_mask = self.create_mask(raw_image)


    def process(self, filename):
        """Process individual images and extract time-dependent data"""

        raw_image = h5py.File(filename)["image"]

        metadata = self.get_time_metadata(raw_image)
        metadata["filename"] = filename

        background = self.estimate_background(raw_image)

        image = self.regrid_image(raw_image)

        if self.remove_background:
            image = image - background

        image[self.elevation<self.elev_cutoff] = 0.

        return image, metadata, background


    def run(self, filelist, numproc=2, seq=False):
        """Run processing on all input files"""

        logging.debug("Processing Images")

        if seq:
            # Process files sequentially in a for loop (primarilly for development)
            results = list()
            for filename in filelist:
                results.append(self.process(filename))
        else:
            # Process files with the multiprocessing module (much faster)
            with multiprocessing.Pool(processes=numproc) as pool:
                results = pool.map(self.process, filelist, chunksize=1)

        self.image_data = np.array([r[0] for r in results])
        self.background = np.array([r[2] for r in results])
        metadata = [r[1] for r in results]
        self.time = np.array([[md["start_time"], md["end_time"]] for md in metadata])
        self.ccd_temp = np.array([md["ccd_temp"] for md in metadata])
        self.metadata["filelist"] = [md["filename"] for md in metadata]

        ut = np.mean(self.time, axis=1)
        self.moon_phase, self.moon_az, self.moon_el = self.moon_position(ut)

        logging.debug("Processing finished")


    def get_metadata(self, image):
        """Extract metadata"""

        metadata = {}
        metadata["station"] = image.attrs["station"]
        metadata["instrument"] = image.attrs["instrument"]
        metadata["label"] = image.attrs["label"]
        metadata["site_lat"] = image.attrs["latitude"]
        metadata["site_lon"] = image.attrs["longitude"]

        return metadata


    def get_time_metadata(self, image):
        """Extract time-dependent metadata"""

        metadata = {}

        start_time = image.attrs["start_time"]
        exposure_time = float(image.attrs["exposure_time"])

        metadata["start_time"] = start_time
        metadata["end_time"] = start_time + exposure_time
        metadata["ccd_temp"] = image.attrs["ccd_temp"]

        return metadata


    def moon_position(self, time):
        """Calculate the Moon phase and position"""

        site_lat = self.metadata["site_lat"]
        site_lon = self.metadata["site_lon"]

        eph = load('de421.bsp')
        ts = load.timescale()
        t = ts.utc(1970, 1, 1, 0, 0, time)

        # Find phase
        phase = almanac.moon_phase(eph, t)

        # Find azimuth and elevation for every time
        earth = eph["earth"]
        site = earth + wgs84.latlon(site_lat, site_lon)
        moon = eph["moon"]
        elev, azmt, dist = site.at(t).observe(moon).apparent().altaz()

        return phase.degrees, azmt.degrees, elev.degrees


    def create_transform_grids(self, image):
        """Create transformation grids"""

        # pylint: disable=too-many-locals

        i_max = image.attrs["width"]
        j_max = image.attrs["height"]

        new_i_max = self.config.getint("PROCESSING", "NEWIMAX")
        new_j_max = self.config.getint("PROCESSING", "NEWJMAX")

        # Create a grid and find the transformed coordinates of each grid point

        x_grid, y_grid = np.meshgrid(np.arange(i_max), np.arange(j_max))

        self.trans_x_grid, self.trans_y_grid = self.transform(
            x_grid, y_grid, self.x0, self.y0, self.rl, self.theta, self.A, self.B, self.C, self.D
        )

        # Find unwarped distance to elevation cutoff and create new grid

        d = self.unwarp(self.elev_cutoff * np.pi / 180.0)

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


    def create_position_arrays(self):
        """Create position arrays"""

        site_lat = self.metadata["site_lat"] * np.pi / 180
        site_lon = self.metadata["site_lon"] * np.pi / 180

        psi = np.sqrt(self.new_x_grid**2 + self.new_y_grid**2) / self.REha

        # Law of Cosines

        c = np.sqrt(RE**2 + self.REha**2 - 2 * RE * self.REha * np.cos(psi))
        b = self.REha - RE * np.cos(psi)

        elv = np.pi / 2.0 - np.arccos(b / c) - psi
        azm = np.arctan2(self.new_x_grid, self.new_y_grid)

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


    def create_mask(self, image):
        """Create mask array"""

        # Mask low elevation angles
        elev_mask = self.elevation < self.elev_cutoff

        # Get azimuth and Elevation of new grids
        elev = np.deg2rad(self.elevation)
        azmt = np.deg2rad(self.azimuth)

        # Calculate r by finding the roots of the cubic lens function
        Delta0 = self.C**2 - 3 * self.D * self.B
        Delta1 = 2 * self.C**3 - 9 * self.D * self.C * self.B + 27 * self.D**2 * (self.A-elev)
        Gamma = ((Delta1 + np.sqrt(Delta1**2 - 4 * Delta0**3)) / 2)**(1./3.)
        r = -(self.C + Gamma + Delta0/Gamma)/(3 * self.D)

        # Invert position and rotatation to original ccd coordinates
        x2 = r*np.sin(azmt)
        y2 = r*np.cos(azmt)

        t = -np.deg2rad(self.theta)
        x1 = np.cos(t) * x2 - np.sin(t) * y2
        y1 = np.sin(t) * x2 + np.cos(t) * y2

        x = x1 * self.rl + self.x0
        y = y1 * self.rl + self.y0

        # Find location of pixels in new image that are outside the original CCD
        i_max = image.attrs["width"]
        j_max = image.attrs["height"]

        x_mask = np.logical_or( x < 0, x > i_max )
        y_mask = np.logical_or( y < 0, y > j_max )
        edge_mask = np.logical_or( x_mask, y_mask )

        # Mask where elevation angle is low or pixel outside original CCD
        mask = np.logical_or(elev_mask, edge_mask)

        return mask



    def estimate_background(self, image):
        """Estimate background from dark corners of CCD"""

        i_max = image.attrs["width"]
        j_max = image.attrs["height"]

        x_grid, y_grid = np.meshgrid(np.arange(i_max), np.arange(j_max))

        mask = np.sqrt((x_grid - self.x0)**2 + (y_grid - self.y0)**2) > self.rl

        masked_image = image[mask]

        m = np.nanmean(masked_image)

        return m



    def regrid_image(self, raw_image):
        """Processing algorithm"""

        cooked_image = np.array(raw_image)

        new_image = griddata(
            (self.trans_x_grid.flatten(), self.trans_y_grid.flatten()),
            cooked_image.flatten(),
            (self.new_x_grid, self.new_y_grid),
            fill_value=0,
            method='linear'
        )

#        fig = plt.figure()
#        ax = fig.add_subplot(121)
#        ax.scatter(self.trans_x_grid, self.trans_y_grid, c=cooked_image, s=0.5)
#        d = self.unwarp(self.elev_cutoff * np.pi / 180.0)
#        ax.axvline(-d, color='magenta')
#        ax.axvline(d, color='magenta')
#        ax.axhline(-d, color='magenta')
#        ax.axhline(d, color='magenta')
#        ax = fig.add_subplot(122)
#        ax.scatter(self.new_x_grid, self.new_y_grid, c=new_image, s=0.5)
#        plt.show()


        return new_image


    def write_to_hdf5(self, output_file):
        """Save results to an HDF5 file"""
    
        site_name = self.config.get("SITE_INFO", "SITE_NAME")
        bright_alt = self.config.getfloat("PROCESSING", "ALTITUDE")
    
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with h5py.File(output_file, "w") as f:

            # Image Data
            images = f.create_dataset(
                "ImageData",
                data=self.image_data.astype('uint16'),
                compression="gzip",
                compression_opts=1,
            )
            images.attrs["Description"] = "pixel values for images"
            images.attrs["Size"] = "Nrecords x Ipixels x Jpixels"
            images.attrs["station"] = self.metadata["station"]
            images.attrs["instrument"] = self.metadata["instrument"]
            images.attrs["remove_background"] = self.remove_background
    
            mask = f.create_dataset(
                "Mask", 
                data=self.image_mask, 
                compression="gzip", 
                compression_opts=1,
            )
            mask.attrs["Description"] = "mask of where ImageData array is corners ouside of camera FoV"
            mask.attrs["Size"] = "Ipixels x Jpixels"
    
            back = f.create_dataset("Background", data=self.background)
            back.attrs["Description"] = "Background brightness estimation from corners"
            back.attrs["Size"] = "Nrecords"
 
            # Time
            t = f.create_dataset(
                "UnixTime",
                data=self.time,
                compression="gzip",
                compression_opts=1,
            )
            t.attrs["Description"] = "unix time stamp"
            t.attrs["Size"] = "Nrecords"
            t.attrs["Unit"] = "seconds"
            
            # Coordinates
            f.create_group("Coordinates")

            lat = f.create_dataset(
                "Coordinates/Latitude", 
                data=self.latitude.astype('float32'), 
                compression="gzip", 
                compression_opts=1,
            )
            lat.attrs["Description"] = "geodetic latitude of each pixel"
            lat.attrs["Size"] = "Ipixels x Jpixels"
            lat.attrs["Projection Altitude"] = bright_alt
            lat.attrs["Unit"] = "degrees"
    
            lon = f.create_dataset(
                "Coordinates/Longitude", 
                data=self.longitude.astype('float32'), 
                compression="gzip", 
                compression_opts=1,
            )
            lon.attrs["Description"] = "geodetic longitude of each pixel"
            lon.attrs["Size"] = "Ipixels x Jpixels"
            lon.attrs["Projection Altitude"] = bright_alt
            lon.attrs["Unit"] = "degrees"
    
            az = f.create_dataset(
                "Coordinates/Azimuth", 
                data=self.azimuth.astype('float32'), 
                compression="gzip", 
                compression_opts=1,
            )
            az.attrs["Description"] = "azimuth of each pixel"
            az.attrs["Size"] = "Ipixels x Jpixels"
            az.attrs["Unit"] = "degrees"
    
            el = f.create_dataset(
                "Coordinates/Elevation", 
                data=self.elevation.astype('float32'), 
                compression="gzip", 
                compression_opts=1,
            )
            el.attrs["Description"] = "elevation of each pixel"
            el.attrs["Size"] = "Ipixels x Jpixels"
            el.attrs["Unit"] = "degrees"
    
            pc = f.create_dataset(
                "Coordinates/PixelCoordinates",
                data=np.array([self.new_x_grid, self.new_y_grid]).astype('float32'),
                compression="gzip",
                compression_opts=1,
            )
            pc.attrs["Description"] = "coordinates of each pixel in grid at the airglow altitude"
            pc.attrs["Size"] = "2 (X,Y) x Ipixels x Jpixels"
            pc.attrs["Projection Altitude"] = bright_alt
            pc.attrs["Unit"] = "km"
    
            # Site Info
            f.create_group("SiteInfo")

            name = f.create_dataset("SiteInfo/Name", data=site_name)
            name.attrs["Description"] = "site name"
    
            code = f.create_dataset("SiteInfo/Code", data=self.metadata["station"])
            code.attrs["Description"] = "three letter site abbreviation/code"
    
            code = f.create_dataset("SiteInfo/Instrument", data=self.metadata["label"])
            code.attrs["Description"] = "type of instrument"
    
            coord = f.create_dataset(
                "SiteInfo/Coordinates",
                data=[self.metadata["site_lat"], self.metadata["site_lon"]],
            )
            coord.attrs["Description"] = "geodetic coordinates of site; [lat, lon]"
            coord.attrs["Unit"] = "degrees"
    
            # Data Quality
            f.create_group("DataQuality")

            ccd = f.create_dataset("DataQuality/CCDTemperature", data=self.ccd_temp)
            ccd.attrs["Description"] = "Temperature of CCD"
            ccd.attrs["Size"] = "Nrecords"
            ccd.attrs["Unit"] = "degrees C"

            phase = f.create_dataset("DataQuality/MoonPhase", data=self.moon_phase)
            phase.attrs["Description"] = "Moon Phase in degrees; New = 0, First Quarter = 90, Full = 180, Last Quarter = 270"
            phase.attrs["Size"] = "Nrecords"
            phase.attrs["Unit"] = "degrees"

            az = f.create_dataset("DataQuality/MoonAzimuth", data=self.moon_az)
            az.attrs["Description"] = "Azimuth of Moon position"
            az.attrs["Size"] = "Nrecords"
            az.attrs["Unit"] = "degrees"

            el = f.create_dataset("DataQuality/MoonElevation", data=self.moon_el)
            el.attrs["Description"] = "Elevation of Moon position"
            el.attrs["Size"] = "Nrecords"
            el.attrs["Unit"] = "degrees"

            # Processing Info
            f.create_group("ProcessingInfo")

            ec = f.create_dataset("ProcessingInfo/ElevationCutoff", data=self.elev_cutoff)
            ec.attrs["Description"] = "elevation angle cutoff [deg]"
    
            ha = f.create_dataset("ProcessingInfo/Altitude", data=bright_alt)
            ha.attrs["Description"] = "assumed altitude of brightness layer [km]"

            timestamp = dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            ts = f.create_dataset("ProcessingInfo/TimeStamp", data=timestamp)
            ts.attrs["Description"] = "timestamp when this file was processed"

            ver = f.create_dataset("ProcessingInfo/Version", data=metadata.version("mangonetwork-raw"))
            ver.attrs["Description"] = "version of processing code"

            pv = f.create_dataset("ProcessingInfo/PythonVersion", data=platform.python_version())
            pv.attrs["Description"] = "python version"

            rf = f.create_dataset("ProcessingInfo/RawFileList", data=self.metadata["filelist"])
            rf.attrs["Description"] = "list of raw files used to generate this dataset"

            with io.StringIO() as buffer:
                self.config.write(buffer)
                config_str = buffer.getvalue()
            cfg = f.create_dataset("ProcessingInfo/ConfigFile", data=config_str)
            cfg.attrs["Description"] = "full text of config file"

            # Platform Info
            f.create_group("ProcessingInfo/PlatformInfo")
            f.create_dataset("ProcessingInfo/PlatformInfo/MachineType", data=platform.machine())
            f.create_dataset("ProcessingInfo/PlatformInfo/System", data=platform.system())
            f.create_dataset("ProcessingInfo/PlatformInfo/Release", data=platform.release())
            f.create_dataset("ProcessingInfo/PlatformInfo/Version", data=platform.version())
            f.create_dataset("ProcessingInfo/PlatformInfo/HostName", data=platform.node())


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
        default="mango-l1.hdf5",
        help="Output filename (default is mango-l1.hdf5)",
    )
    parser.add_argument(
        "-n", 
        "--numproc", 
        type=int, 
        default=1, 
        help="Number of parallel processes"
    )
    parser.add_argument(
        "-s", 
        "--sequential", 
        action="store_true", 
        help="Process images sequentially without multiprocessing"
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

    # Placeholder for default config file location
    #   This function can be rewritten later
    config_dir = os.environ['MANGONETWORK_CONFIGS']

    config_file = os.path.join(config_dir, f"{station}-{instrument}.ini")

    logging.debug("Using package configuration file: %s", config_file)

    #name = f"{station}-{instrument}.ini"

    #logging.debug("Using package configuration file: %s", name)

    #return resources.files("mangonetwork.raw.data").joinpath(name).read_text()
    return config_file



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
        if os.path.exists(args.config):
            config_file = args.config
        else:
            logging.error("Config file not found")
            sys.exit(1)
    else:
        config_file = find_config(inputs[0])

    config = configparser.ConfigParser()
    config.read(config_file)

    # Process images
    processor = ImageProcessor(config)
    processor.setup(inputs[0])
    processor.run(inputs, numproc=args.numproc, seq=args.sequential)
    output_file = pathlib.Path(args.output)
    processor.write_to_hdf5(output_file)

    sys.exit(0)


if __name__ == "__main__":
    main()
