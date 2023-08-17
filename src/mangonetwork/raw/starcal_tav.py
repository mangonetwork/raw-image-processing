#!/usr/bin/env python
"""Plot Star Calibration files"""

##########################################################################
#
#   Plot starcal file
#
#   2022-xx-xx  Leslie Lamarche and Asti Bhatt
#               Initial implementation
#
#   2023-03-08  Todd Valentic
#               PEP8 compliance
#
##########################################################################

import argparse
import io
import logging
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt

from . import imageops

if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources


class StarCal:
    """Star calibration"""

    def __init__(self, starcal_file):
        self.plot_stars(starcal_file)

    def plot_stars(self, starcal_file):
        """Plot stars"""

        az, el, i, j, raw_file = self.parse_file(starcal_file)
        el = el * np.pi / 180.0
        az = az * np.pi / 180.0

        contrast = 99.95
        rotation_angle = 0.0

        image = h5py.File(raw_file, "r")["image"]

        cooked_image = np.array(image)
        cooked_image = imageops.equalize(cooked_image, contrast)
        cooked_image = imageops.rotate(cooked_image, rotation_angle)

        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(111)
        ax.imshow(cooked_image.imageData, cmap="gray")
        ax.scatter(i, j, facecolors="none", edgecolors="r")

        plt.show()

    def parse_file(self, starcal_file):
        """Read starcal file"""

        raw_filename = starcal_file.split("\n")[0].split()[-1]

        az, el, i, j = np.loadtxt(
            io.StringIO(starcal_file), usecols=(1, 2, 3, 4), unpack=True
        )

        return az, el, i, j, raw_filename


# ------------------------------------------------------------------------
# Main application
# ------------------------------------------------------------------------


def parse_args():
    """Command line parsing"""

    parser = argparse.ArgumentParser(
        description="Check the star identification for calibration"
    )

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")

    parser.add_argument(
        "-s", "--starcal", metavar="FILE", help="Alternate starcal file"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    return parser.parse_args()


def find_starcal(station, instrument):
    """Find starcal file in package data"""

    starcal_file = f"starcal-{station}-{instrument}.txt"

    logging.debug("Using package starcal file: %s", starcal_file)

    return resources.files("mangonetwork.raw.data").joinpath(starcal_file).read_text()


def main():
    """Main application"""

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.starcal:
        logging.debug("Alternate starcal file: %s", args.starcal)
        if not os.path.exists(args.starcal):
            logging.error("Config file not found")
            sys.exit(1)
        with open(args.starcal, encoding="utf-8") as f:
            contents = f.read()
    else:
        contents = find_starcal(args.station, args.instrument)

    StarCal(contents)

    sys.exit(0)


if __name__ == "__main__":
    main()
