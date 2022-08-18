#!/usr/bin/env python3

import os
import sys
import configparser
import argparse
import logging
import datetime
import pathlib

import h5py
import hcipy
import matplotlib.pyplot as plt
import numpy as np

from .MANGOImage import MANGOImage

if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources

class QuickLook:

    def __init__(self, config, inputList, outputFile):

        self.load_config(config)
        self.process_images(inputList, outputFile)

    def load_config(self, config):

        self.siteName = config.get('PROCESSING','SITE_NAME')
        self.siteState = config.get('PROCESSING','SITE_STATE')
        self.contrast = config.getint('PROCESSING','CONTRAST')
        self.rotationAngle = config.getfloat('QUICKLOOK','ROTATION_ANGLE')

    def process_images(self, inputList, outputFile):

        outputFile = os.path.abspath(outputFile)
        outputPath = os.path.split(outputFile)[0]

        pathlib.Path(outputPath).mkdir(parents=True, exist_ok=True)

        imageWriter = hcipy.plotting.FFMpegWriter(outputFile, framerate = 10)

        for filename in inputList:
            logging.debug(filename)
            self.plot(filename, imageWriter)

        imageWriter.close()

    def plot(self, filename, imageWriter):

        wavelength = {'': 'Unknown', 'Green Line': '557.7 nm', 'Red Line': '630.0 nm'}

        image = h5py.File(filename, 'r')['image']

        cooked_image = MANGOImage(np.array(image))
        cooked_image.equalizeHistogram(self.contrast)
        cooked_image.rotateImage(self.rotationAngle)

        start_time = datetime.datetime.utcfromtimestamp(image.attrs['start_time'])
        exposure_time = image.attrs['exposure_time']
        ccd_temp = image.attrs['ccd_temp']
        code = image.attrs['station']
        label = image.attrs['label']
        latlon = '%4.1f N, %5.1f W' % (image.attrs['latitude'],image.attrs['longitude'])

        fig, ax = plt.subplots(facecolor='black')
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

        ax.imshow(cooked_image.imageData, cmap='gray')

        ny, nx = cooked_image.imageData.shape
        dy = 20
        lx = 10
        rx = nx-lx
        by = 470

        ax.annotate('N', xy=(nx/2, dy), ha='center', color='white')
        ax.annotate('E', xy=(nx-60, ny/2), ha='right', color='white')
        ax.annotate('%s' % start_time.date(), xy=(lx, dy), color='white')
        ax.annotate('%s UTC' % start_time.time(), xy=(lx, 2*dy), color='white')

        ax.annotate(label, xy=(rx,dy), color='white', ha='right')
        ax.annotate(wavelength[label], xy=(rx,2*dy), color='white', ha='right')
        ax.annotate('%+5.1f C' % ccd_temp, xy=(rx, 3*dy), color='white', ha='right')
        ax.annotate('%s s' % exposure_time, xy=(rx,4*dy), color='white', ha='right')

        ax.annotate('%s - %s' % (code.upper(), self.siteState), xy=(lx,by), color='white')
        ax.annotate(latlon, xy=(lx,by+dy), color='white')
        ax.annotate(self.siteName, xy=(lx,by+2*dy), color='white')

        ax.annotate('NSF/SRI MANGO DASI', xy=(rx,by+2*dy), color='white', ha='right')

        imageWriter.add_frame(fig)
        plt.close()

def parse_args():

    parser = argparse.ArgumentParser(description='Create a quick look movie')

    parser.add_argument('-c', '--config', metavar='FILE',
                        help='Alternate configuration file')
    parser.add_argument('-f', '--filelist', metavar='FILE',
                        help='A file with a list of .hdf5 file names')
    parser.add_argument('-o', '--output', default='mango.mp4',
                        help='Output filename (default is mango.mp4)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')

    parser.add_argument('inputfiles', nargs='*')

    return parser.parse_args()

def find_config(filename):

    with h5py.File(filename) as h5:
        station = h5['image'].attrs['station']
        instrument = h5['image'].attrs['instrument']

    name = '%s-%s.ini' % (station, instrument)

    logging.debug('Using package configuration file: %s' % name)

    return resources.files('mangonetwork.raw.data').joinpath(name).read_text()

def main():

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.filelist:
        inputs = [line.strip() for line in open(args.filelist)]
        inputs = [line for line in inputs if line]
    else:
        inputs = args.inputfiles

    if not inputs:
        logging.error('No input files listed')
        sys.exit(1)

    logging.debug('Processing %d files' % len(inputs))

    if args.config:
        logging.debug('Alternate configuration file: %s' % args.config)
        if not os.path.exists(args.config):
            logging.error('Config file not found')
            sys.exit(1)
        contents = open(args.config).read()
    else:
        contents = find_config(inputs[0])

    config = configparser.ConfigParser()
    config.read_string(contents)

    logging.debug('Configuration file: %s' % args.config)

    QuickLook(config, inputs, args.output)

    sys.exit(0)

if __name__ == '__main__':
    main()
