# StarCal.py
# Create/Validate Star Calibration files

import datetime
import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging
import io
import sys
if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources

from .MANGOImage import MANGOImage


class StarCal:

    def __init__(self, starCalFile):

        self.plot_stars(starCalFile)

    def plot_stars(self, starCalFile):

        az, el, i, j, raw_file = self.parse_file(starCalFile)
        el = el*np.pi/180.
        az = az*np.pi/180.

        contrast = 99.95
        rotationAngle = 0.

        image = h5py.File(raw_file, 'r')['image']

        cooked_image = MANGOImage(np.array(image))
        cooked_image.equalizeHistogram(contrast)
        cooked_image.rotateImage(rotationAngle)


        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)
        ax.imshow(cooked_image.imageData, cmap='gray')
        ax.scatter(i, j, facecolors='none', edgecolors='r')

        plt.show()

    def parse_file(self, starCalFile):

        raw_filename = starCalFile.split('\n')[0].split()[-1]

        az, el, i, j = np.loadtxt(io.StringIO(starCalFile), usecols=(1,2,3,4), unpack=True)

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

    starcal_file = 'starcal-%s-%s.txt' % (station, instrument)

    logging.debug('Using package starcal file: %s' % starcal_file)

    return resources.files('mangonetwork.raw.data').joinpath(starcal_file).read_text()

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
        contents = open(args.starcal).read()
    else:
        contents = find_starcal(args.station, args.instrument)

    StarCal(contents)

    sys.exit(0)

if __name__ == '__main__':
    main()
