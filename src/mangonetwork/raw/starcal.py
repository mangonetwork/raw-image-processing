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

import datetime as dt
import argparse
import io
import logging
import os
import sys

import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_point_clicker import clicker

from skyfield.api import Star, load, wgs84
from skyfield.data import hipparcos

from . import imageops

if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources


class StarCal:
    """Star calibration"""

    def __init__(self, star_cal_file, output):

        self.display_image(star_cal_file)
        self.lookup_stars()
        self.save_starcal_file(output, star_cal_file)

    def add_star(self, position, klass):
        x, y = position
        print(f"Star at {x=:02f}, {y=:02f}")
        hip = input('HIP #: ')
        self.star_hip.append(hip)

    def prep_image(self, image, contrast=99.95, rotation_angle=0.):
        """Prepare image to display"""

        cooked_image = np.array(image)
        cooked_image = imageops.equalize(cooked_image, contrast)
        cooked_image = imageops.rotate(cooked_image, rotation_angle)

        return cooked_image


    def display_image(self, star_cal_file):
        """Display image and set up klicker"""

        az, el, i, j, raw_file = self.parse_file(star_cal_file)

        image = h5py.File(raw_file, 'r')['image']
        cooked_image = self.prep_image(image)

        self.time = dt.datetime.utcfromtimestamp(image.attrs['start_time'])
        self.site_lat = image.attrs['latitude']
        self.site_lon = image.attrs['longitude']

        # Display image with stars
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)
        ax.imshow(cooked_image, cmap='gray')
        ax.scatter(i, j, facecolors='none', edgecolors='r')

        # Setup clicker object to keep track of identified stars
        self.klicker = clicker(ax, ['stars'])
        self.star_hip = list()
        self.klicker.on_point_added(self.add_star)

        plt.show()

        self.star_pos = self.klicker.get_positions()['stars']
       

    def lookup_stars(self):
        """Look up star position"""
        ts = load.timescale()
        t = ts.utc(self.time.year,self.time.month,self.time.day,self.time.hour,self.time.minute,self.time.second)
        planets = load('de421.bsp')
        earth = planets['earth']
        site = earth + wgs84.latlon(self.site_lat, self.site_lon, elevation_m=0)
        site_ref = site.at(t)

        with load.open(hipparcos.URL) as f:
            df = hipparcos.load_dataframe(f)

        self.star_azel = list()
        for hip in self.star_hip:
            s = Star.from_dataframe(df.loc[float(hip)])

            elev, azmt, _ = site_ref.observe(s).apparent().altaz()
            self.star_azel.append([azmt.degrees, elev.degrees])


    def save_starcal_file(self, output, star_cal_file):
        """ Save output starcal file"""

        with open(output, 'w') as f:
            # copy existing file
            f.write(star_cal_file)

            # add new stars
            for hip, azel, pos in zip(self.star_hip, self.star_azel, self.star_pos):
                f.write(f'{hip}    {azel[0]}    {azel[1]}    {pos[0]}    {pos[1]}\n')

    def parse_file(self, star_cal_file):
        """Read starcal file"""

        raw_filename = star_cal_file.split('\n')[0].split()[-1]

        az, el, i, j = np.loadtxt(io.StringIO(star_cal_file), usecols=(1,2,3,4), unpack=True)

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
        "-sc", "--starcal", metavar="FILE", help="Alternate starcal file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="starcal-mango.txt",
        help="Output starcal filename (default is starcal-mango.txt)",
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

    StarCal(contents, args.output)

    sys.exit(0)


if __name__ == "__main__":
    main()
