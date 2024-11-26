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
import hcipy
import matplotlib.pyplot as plt
import numpy as np

from . import imageops

if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources


class QuickLook:
    """Quicklook image processor"""

    def __init__(self, config, input_list, args):
        self.args = args
        output_file = pathlib.Path(args.output)

        self.load_config(config)
        self.process_images(input_list, output_file)

    def load_config(self, config):
        """Get configuration values"""

        self.site_name = config.get("SITE_INFO", "SITE_NAME")
        self.site_state = config.get("SITE_INFO", "SITE_STATE")
        self.remove_background = config.getboolean("QUICKLOOK", "REMOVE_BACKGROUND")
        self.contrast = config.getfloat("QUICKLOOK", "CONTRAST", fallback=False)

        try:
            self.rotation_angle = config.getfloat("CALIBRATION_PARAMS", "THETA")
        except configparser.NoOptionError:
            self.rotation_angle = config.getfloat("QUICKLOOK", "MANUAL_THETA")
            logging.debug("Using Manual Rotation Angle: %s", self.rotation_angle)

    def process_images(self, input_list, output_file):
        """Process input images"""

        output_file.parent.mkdir(parents=True, exist_ok=True)

        image_writer = hcipy.plotting.FFMpegWriter(output_file, framerate=10)

        for filename in input_list:
            logging.debug(filename)
            self.plot(filename, image_writer)

        image_writer.close()

    # pylint: disable=too-many-locals

    def plot(self, filename, image_writer):
        """Process and annotate image"""

        wavelength = {"": "Unknown", "Green Line": "557.7 nm", "Red Line": "630.0 nm"}

        image = h5py.File(filename, "r")["image"]

        cooked_image = np.array(image)
        if self.remove_background:
            cooked_image = imageops.background_removal(cooked_image)
        cooked_image = imageops.equalize(cooked_image, self.contrast)
        cooked_image = imageops.invert(cooked_image)
        cooked_image = imageops.rotate(cooked_image, self.rotation_angle)

        start_time = datetime.datetime.utcfromtimestamp(image.attrs["start_time"])
        exposure_time = image.attrs["exposure_time"]
        ccd_temp = image.attrs["ccd_temp"]
        code = image.attrs["station"]
        label = image.attrs["label"]
        latlon = f"{image.attrs['latitude']:4.1f} N, {image.attrs['longitude']:5.1f} W"

        fig, ax = plt.subplots(facecolor="black")
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

        ax.imshow(cooked_image, cmap="gray")

        ny, nx = cooked_image.shape
        dy = 20
        lx = 10
        rx = nx - lx
        by = 470

        ax.annotate("N", xy=(nx / 2, dy), ha="center", color="white")
        ax.annotate("E", xy=(nx - 60, ny / 2), ha="right", color="white")
        ax.annotate(f"{start_time.date()}", xy=(lx, dy), color="white")
        ax.annotate(f"{start_time.time()} UTC", xy=(lx, 2 * dy), color="white")

        ax.annotate(label, xy=(rx, dy), color="white", ha="right")
        ax.annotate(wavelength[label], xy=(rx, 2 * dy), color="white", ha="right")
        ax.annotate(f"{ccd_temp:+5.1f} C", xy=(rx, 3 * dy), color="white", ha="right")
        ax.annotate(f"{exposure_time} s", xy=(rx, 4 * dy), color="white", ha="right")

        ax.annotate(f"{code.upper()} - {self.site_state}", xy=(lx, by), color="white")
        ax.annotate(latlon, xy=(lx, by + dy), color="white")
        ax.annotate(self.site_name, xy=(lx, by + 2 * dy), color="white")

        ax.annotate(
            "NSF/SRI MANGO DASI", xy=(rx, by + 2 * dy), color="white", ha="right"
        )

        image_writer.add_frame(fig)
        plt.close()


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

    name = f"{station}-{instrument}.ini"

    logging.debug("Using package configuration file: %s", name)

    return resources.files("mangonetwork.raw.data").joinpath(name).read_text()


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

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    inputs = find_inputfiles(args)

    if not inputs:
        logging.error("No input files listed")
        sys.exit(1)

    logging.debug("Processing %d files", len(inputs))

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

    QuickLook(config, inputs, args)

    sys.exit(0)


if __name__ == "__main__":
    main()
