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
        # cooked_image.invertImage()
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
        ax.scatter(i, j, facecolors='none', edgecolors='r')

        xp = 519./2.*np.cos(el)*np.sin(az)+695./2.
        yp = 519./2.*np.cos(el)*np.cos(az)+519./2.
        ax.scatter(xp, yp)

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

        az, el, i, j = np.loadtxt(io.StringIO(starCalFile), usecols=(1,2,3,4), unpack=True)

        # data = list()
        # for line in starCalFile.split('\n'):
        #     if line.startswith('#') or len(line)==0:
        #         continue
        #     data.append([float(x) for x in line.split()[1:5]])
        #
        # az, el, i, j = np.array(data).T

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

    # with h5py.File(filename) as h5:
    #     station = h5['image'].attrs['station']
    #     instrument = h5['image'].attrs['instrument']

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

    # logging.debug('Configuration file: %s' % args.starcal)

    StarCal(contents)

    sys.exit(0)

if __name__ == '__main__':
    main()
