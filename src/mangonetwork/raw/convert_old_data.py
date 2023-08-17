# ConvertOldData.py
# Convert old binary image format to new raw hdf5 image format


import os
import sys
import configparser
import argparse
import logging
import datetime

import h5py
import hcipy
import matplotlib.pyplot as plt
import numpy as np

from .MANGOImage import MANGOImage

if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources


# Do we want to be able to specify output file names, or just use default patterns?
    # Probably both
# Should this file take in a list of input/ouput files, or just do one at a time?

class ConvertOldData:

    def __init__(self, config, inputList, outputFile):

        self.load_config(config)        # Might get rid of this
        self.convert_files(inputList, outputFile)

    def load_config(self, config):
        # From config file, load any "extra" metadata parameters that aren't included in the binary file
        pass
        # self.siteName = config.get('PROCESSING','SITE_NAME')
        # self.siteState = config.get('PROCESSING','SITE_STATE')
        # self.contrast = config.getint('PROCESSING','CONTRAST')
        # self.rotationAngle = config.getfloat('QUICKLOOK','ROTATION_ANGLE')

    def convert_files(self, inputList, outputFile):

        outputFile = os.path.abspath(outputFile)
        outputPath = os.path.split(outputFile)[0]

        if not os.path.exists(outputPath):
            os.makedirs(outputPath)

        for filename in inputList:
            logging.debug(filename)
            self.plot(filename, imageWriter)

        imageWriter.close()

    def convert(filename)


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

    logging.debug('Converting %d files' % len(inputs))

    # CONFIG FILE MAY OR MAY NOT BE NESSISARY
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

    ConvertOldData(config, inputs, args.output)

    sys.exit(0)

if __name__ == '__main__':
    main()
