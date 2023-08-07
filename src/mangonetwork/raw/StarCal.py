# StarCal.py
# Create/Validate Star Calibration files

import datetime as dt
import h5py
import matplotlib.pyplot as plt
from mpl_point_clicker import clicker
from typing import Tuple
import numpy as np
import argparse
import logging
import io
import os
import sys
if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources

from .MANGOImage import MANGOImage


from skyfield.api import Star, load, wgs84
#from skyfield.api import N,S,E,W, wgs84
from skyfield.data import hipparcos


class StarCal:

    def __init__(self, starCalFile):

        #self.load_image()
        self.display_image(starCalFile)
        self.lookup_stars()
        self.save_starcal_file(starCalFile)
        #self.plot_stars(starCalFile)

    def add_star(self, position: Tuple[float, float], klass: str):
        x, y = position
        print(f"Star at {x=:02f}, {y=:02f}")
        hip = input('HIP #: ')
        self.star_hip.append(hip)

    def prep_image(self, image):

        contrast = 99.95
        rotationAngle = 0.

        cooked_image = MANGOImage(np.array(image))
        cooked_image.equalizeHistogram(contrast)
        cooked_image.rotateImage(rotationAngle)

        return cooked_image


    def display_image(self, starCalFile):

        az, el, i, j, raw_file = self.parse_file(starCalFile)
        #el = el*np.pi/180.
        #az = az*np.pi/180.


        # Function prepare_image() ?
        #contrast = 99.95
        #rotationAngle = 0.

        image = h5py.File(raw_file, 'r')['image']

        #cooked_image = MANGOImage(np.array(image))
        #cooked_image.equalizeHistogram(contrast)
        #cooked_image.rotateImage(rotationAngle)
        cooked_image = self.prep_image(image)

        self.time = dt.datetime.utcfromtimestamp(image.attrs['start_time'])
        self.site_lat = image.attrs['latitude']
        self.site_lon = image.attrs['longitude']

        # Display image with stars
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)
        ax.imshow(cooked_image.imageData, cmap='gray')
        ax.scatter(i, j, facecolors='none', edgecolors='r')

        # Setup clicker object to keep track of identified stars
        self.klicker = clicker(ax, ['stars'])
        self.star_hip = list()

#        def point_added_cb(position: Tuple[float, float], klass: str):
#            x, y = position
#            print(f"Star at {x=:02f}, {y=:02f}")
#            hip = input('HIP #: ')
#            stars.append(hip)

        self.klicker.on_point_added(self.add_star)

        plt.show()

        self.star_pos = self.klicker.get_positions()['stars']
       

        # Function lookup_star() ?
    def lookup_stars(self):
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
            #print(hip, elev.degrees, azmt.degrees)


        # Need to append output to starcalfile
        # may be tricky as this is currently treated as package data

    def save_starcal_file(self, starCalFile):

        with open(starCalFile, 'a') as f:
            for hip, azel, pos in zip(self.star_hip, self.star_azel, self.star_pos):
                f.write(f'{hip}    {azel[0]}    {azel[1]}    {pos[0]}    {pos[1]}\n')

    def parse_file(self, starCalFile):

        with open(starCalFile, 'r') as f:
            line = f.readline()
            raw_filename = line.split()[-1]
        #raw_filename = starCalFile.split('\n')[0].split()[-1]

        az, el, i, j = np.loadtxt(starCalFile, usecols=(1,2,3,4), unpack=True)
        #az, el, i, j = np.loadtxt(io.StringIO(starCalFile), usecols=(1,2,3,4), unpack=True)

        return az, el, i, j, raw_filename


def parse_args():

    parser = argparse.ArgumentParser(description='Check the star identification for calibration')

    parser.add_argument('station', help='Station code')
    parser.add_argument('instrument', help='redline or greenline')

    parser.add_argument('-sc', '--starcal', metavar='FILE',
                        help='Alternate starcal file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')

    return parser.parse_args()

def find_starcal(station, instrument):

    path = os.path.dirname(__file__)

    starcal_file = 'starcal-%s-%s.txt' % (station, instrument)

    logging.debug('Using package starcal file: %s' % starcal_file)

    full_path = os.path.join(path, 'data', starcal_file)

    return full_path
    #return resources.files('mangonetwork.raw.data').joinpath(starcal_file).read_text()

def main():

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    if args.starcal:
        logging.debug('Alternate starcal file: %s' % args.starcal)
        if not os.path.exists(args.starcal):
            logging.error('Config file not found')
            sys.exit(1)
        #contents = open(args.starcal).read()
        starcal_path = args.starcal
    else:
        starcal_path = find_starcal(args.station, args.instrument)

    StarCal(starcal_path)

    sys.exit(0)

if __name__ == '__main__':
    main()
