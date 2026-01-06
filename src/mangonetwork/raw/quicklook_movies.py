#!/usr/bin/env python3
"""Quicklook movie generator"""

##########################################################################
#
#   Generate quicklook movies from raw images.
#
#   2022-xx-xx  Todd Valentic
#               Refactored from original code.
#
#   2023-03-07  Todd Valentic
#               Refactor for PEP8 compliance
#
##########################################################################

import argparse
import configparser
import datetime
import logging
import os
import pathlib
import sys

import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import skimage.transform



if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources


class QuickLook:
    """Quicklook image processor"""

    def __init__(self, config):
        self.load_config(config)

    def load_config(self, config):
        """Get configuration values"""

        self.site_name = config.get("SITE_INFO", "SITE_NAME")
        self.site_state = config.get("SITE_INFO", "SITE_STATE")
        self.remove_background = config.getboolean("QUICKLOOK", "REMOVE_BACKGROUND", fallback=True)
        self.contrast = config.getfloat("QUICKLOOK", "CONTRAST", fallback=False)

        try:
            self.rotation_angle = config.getfloat("CALIBRATION_PARAMS", "THETA")
        except configparser.NoOptionError:
            self.rotation_angle = config.getfloat("QUICKLOOK", "MANUAL_THETA")
            logging.debug("Using Manual Rotation Angle: %s", self.rotation_angle)

    def process_images(self, input_list, output_file):
        """Process input images"""

        output_file.parent.mkdir(parents=True, exist_ok=True)

        fig, self.ax = plt.subplots(facecolor="black")
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

        frames = []
        for filename in input_list:
            logging.debug(filename)
            f = self.plot(filename)
            frames.append(f)

        ## Process files with the multiprocessing module (much faster)
        ## Doesn't work
        #with multiprocessing.Pool(processes=4) as pool:
        #    frames = pool.map(self.plot, input_list, chunksize=1)

        ani = animation.ArtistAnimation(fig=fig, artists=frames, interval=100)
        ani.save(output_file)


    # pylint: disable=too-many-locals

    def plot(self, filename):
        """Process and annotate image"""

        wavelength = {"": "Unknown", "Green Line": "557.7 nm", "Red Line": "630.0 nm"}

        image = h5py.File(filename, "r")["image"]

        cooked_image = np.array(image)
        if self.remove_background:
            cooked_image = self.background_removal(cooked_image)
        if self.contrast:
            cooked_image = self.equalize(cooked_image, self.contrast)
        cooked_image = self.rotate(cooked_image, self.rotation_angle)

        start_time = datetime.datetime.utcfromtimestamp(image.attrs["start_time"])
        exposure_time = image.attrs["exposure_time"]
        ccd_temp = image.attrs["ccd_temp"]
        code = image.attrs["station"]
        label = image.attrs["label"]
        latlon = f"{image.attrs['latitude']:4.1f} N, {image.attrs['longitude']:5.1f} W"

        frm = self.ax.imshow(cooked_image, cmap="gray")

        # Collect artists to generate annimation
        artists = [frm]

        ny, nx = cooked_image.shape
        dy = 20
        lx = 10
        rx = nx - lx
        by = 470

        t1 = self.ax.annotate("N", xy=(nx / 2, dy), ha="center", color="white")
        t2 = self.ax.annotate("E", xy=(nx - 60, ny / 2), ha="right", color="white")
        t3 = self.ax.annotate(f"{start_time.date()}", xy=(lx, dy), color="white")
        t4 = self.ax.annotate(f"{start_time.time()} UTC", xy=(lx, 2 * dy), color="white")
        artists.extend([t1, t2, t3, t4])

        t1 = self.ax.annotate(label, xy=(rx, dy), color="white", ha="right")
        t2 = self.ax.annotate(wavelength[label], xy=(rx, 2 * dy), color="white", ha="right")
        t3 = self.ax.annotate(f"{ccd_temp:+5.1f} C", xy=(rx, 3 * dy), color="white", ha="right")
        t4 = self.ax.annotate(f"{exposure_time} s", xy=(rx, 4 * dy), color="white", ha="right")
        artists.extend([t1, t2, t3, t4])

        t1 = self.ax.annotate(f"{code.upper()} - {self.site_state}", xy=(lx, by), color="white")
        t2 = self.ax.annotate(latlon, xy=(lx, by + dy), color="white")
        t3 = self.ax.annotate(self.site_name, xy=(lx, by + 2 * dy), color="white")

        t4 = self.ax.annotate(
            "NSF/SRI MANGO DASI", xy=(rx, by + 2 * dy), color="white", ha="right"
        )
        artists.extend([t1, t2, t3, t4])

        return artists

       

    def equalize(self, image, contrast, num_bins=10000):
        """Histogram Equalization to adjust contrast [1%-99%]"""
    
        image_array_1d = image.flatten()
    
        image_histogram, bins = np.histogram(image_array_1d, num_bins)
        image_histogram = image_histogram[1:]
        bins = bins[1:]
        cdf = np.cumsum(image_histogram)
    
        # spliced to cut off non-image area
        # any way to determine this dynamically?  How periminant is it?
        cdf = cdf[:9996]
    
        max_cdf = max(cdf)
        max_index = np.argmin(abs(cdf - contrast / 100 * max_cdf))
        min_index = np.argmin(abs(cdf - (100 - contrast) / 100 * max_cdf))
        vmax = float(bins[max_index])
        vmin = float(bins[min_index])
        low_value_indices = image_array_1d < vmin
        image_array_1d[low_value_indices] = vmin
        high_value_indices = image_array_1d > vmax
        image_array_1d[high_value_indices] = vmax
    
        return image_array_1d.reshape(image.shape)
    
    def background_removal(self, image):
        """Subtract dark current background from image based on corner brightness"""
    
        # Offset from edge of image and size of region for determining the background in each corner of image
        offx = int(0.015*image.shape[1])
        sizex = int(0.05*image.shape[1])
        offy = int(0.015*image.shape[0])
        sizey = int(0.05*image.shape[0])
    
        # Calculate means in the four corners
        m1 = image[offy:offy+sizey, offx:offx+sizex]
        m2 = image[-(offy+sizey):-offy, -(offx+sizex):-offx]
        m3 = image[offy:offy+sizey, -(offx+sizex):-offx]
        m4 = image[-(offy+sizey):-offy, offx:offx+sizex]
        m = np.nanmedian([m1,m2,m3,m4])
        # Subtract mean and set a 0 floor
        image = image - m
        image[image < 0] = 0
    
        return image
    
    
   
    def rotate(self, image, angle):
        """Rotate image"""
    
        # MUST flip image before rotating
        # Images are captured from below, but visualized from above in most standard formats
        img_flip =  np.flipud(image)

        # Rotate
        img_rot = skimage.transform.rotate(img_flip, angle, order=3).astype(float)

        return img_rot

    


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
        default="mango.mp4",
        help="Output filename (default is mango.mp4)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    parser.add_argument("inputfiles", nargs="*")

    return parser.parse_args()


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
        logging.error("No input files listed")
        sys.exit(1)

    logging.debug("Processing %d files", len(inputs))

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

    # Generate quicklook movie
    movie = QuickLook(config)
    output_file = pathlib.Path(args.output)
    movie.process_images(inputs, output_file)

    sys.exit(0)


if __name__ == "__main__":
    main()
