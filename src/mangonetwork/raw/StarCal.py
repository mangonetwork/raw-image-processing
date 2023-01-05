# StarCal.py
# Create/Validate Star Calibration files

import datetime
import h5py
import matplotlib.pyplot as plt
import numpy as np
import argparse
import logging
import sys
if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources

from MANGOImage import MANGOImage


class StarCal:

    def __init__(self, starCalFile):

        # raw_file = '/Users/e30737/Desktop/Data/MANGO/blo/greenline/mango-blo-greenline-20211129-061400.hdf5'
        # # starCalFile = '/Users/e30737/Desktop/Projects/MANGO/star_cal_files/BLO_green_20211129_061400.txt'
        # starCalFile = '/Users/e30737/Desktop/Projects/MANGO/raw-image-processing/src/mangonetwork/raw/data/starcal_blo_green.txt'
        # # raw_file = '/Users/e30737/Desktop/Data/MANGO/low/greenline/2021/333/04/mango-low-greenline-20211129-043600.hdf5'
        # # starCalFile = '/Users/e30737/Desktop/Projects/MANGO/star_cal_files/LOW_green_20211129_043600.txt'

        self.plot_stars(starCalFile)

    def plot_stars(self, starCalFile):

        az, el, i, j, raw_file = self.parse_file(starCalFile)
        el = el*np.pi/180.
        az = az*np.pi/180.

        contrast = 99.95
        rotationAngle = 0.

        # raw_file = '/Users/e30737/Desktop/Data/MANGO/blo/greenline/mango-blo-greenline-20211129-061400.hdf5'

        image = h5py.File(raw_file, 'r')['image']
        print(image.shape)

        cooked_image = MANGOImage(np.array(image))
        cooked_image.equalizeHistogram(contrast)
        cooked_image.rotateImage(rotationAngle)

        # start_time = datetime.datetime.utcfromtimestamp(image.attrs['start_time'])
        # exposure_time = image.attrs['exposure_time']
        # ccd_temp = image.attrs['ccd_temp']
        # code = image.attrs['station']
        # label = image.attrs['label']
        # latlon = '%4.1f N, %5.1f W' % (image.attrs['latitude'],image.attrs['longitude'])

        # print(starCalFile)
        # el, az, i, j = np.loadtxt(starCalFile, usecols=(1,2,3,4), unpack=True)
        # print(min(i), max(i), min(j), max(j))


        # fig, ax = plt.subplots()
        fig = plt.figure(figsize=(15,10))
        ax = fig.add_subplot(111)
        # ax.imshow(cooked_image.imageData, cmap='gray', origin='upper')
        ax.imshow(cooked_image.imageData, cmap='gray')
        ax.scatter(695-i, j, facecolors='none', edgecolors='r')

        # ax = fig.add_subplot(122)
        # xp = np.cos(el)*np.sin(az)
        # yp = np.cos(el)*np.cos(az)
        # ax.scatter(xp, yp)
        # ax.set_xlim([-1.0,1.0])
        # ax.set_ylim([-1.0,1.0])
        # ax.set_aspect('equal')

        plt.show()

    def parse_file(self, starCalFile):

        raw_filename = starCalFile.split('\n')[0].split()[-1]

        data = list()
        for line in starCalFile.split('\n'):
            if line.startswith('#') or len(line)==0:
                continue
            data.append([float(x) for x in line.split()[1:5]])

        az, el, i, j = np.array(data).T

        return az, el, i, j, raw_filename


def parse_args():

    parser = argparse.ArgumentParser(description='Check the star identification for calibration')

    parser.add_argument('imager', help='Imager code (ie. CFSred or BLOgreen)')
    parser.add_argument('-c', '--config', metavar='FILE',
                        help='Alternate starcal file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')

    return parser.parse_args()

def find_config(imager):

    # with h5py.File(filename) as h5:
    #     station = h5['image'].attrs['station']
    #     instrument = h5['image'].attrs['instrument']

    name = 'starcal-%s.txt' % (imager)

    logging.debug('Using package starcal file: %s' % name)

    return resources.files('mangonetwork.raw.data').joinpath(name).read_text()

def main():

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    if args.config:
        logging.debug('Alternate starcal file: %s' % args.config)
        if not os.path.exists(args.config):
            logging.error('Config file not found')
            sys.exit(1)
        contents = open(args.config).read()
    else:
        contents = find_config(args.imager)

    logging.debug('Configuration file: %s' % args.config)

    StarCal(contents)

    sys.exit(0)

if __name__ == '__main__':
    main()
