# QuicklookMovies.py

import os
import warnings
import configparser

import matplotlib.pyplot as plt

from . import MANGOimage
import argparse
import pathlib
import re
import pandas as pd
import datetime as dt
import numpy as np
import h5py
import skimage.transform
import hcipy

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")


class QuickLook:

    def __init__(self, configFile, inputList, outputFile):
        self.configFile = configFile
        self.inputList = inputList
        self.outputFile = outputFile
        self.build_config()
        self.process_images()

    def build_config(self):
        # read in config file
        config = configparser.ConfigParser()
        config.read(self.configFile)

        self.contrast = int(config['Specifications']['contrast'])
        self.calFile = config['Data Locations']['cal_hdf']
        with h5py.File(self.calFile, 'r') as f:
            self.rotationAngle = f['Calibration Angle'][()]


    def process_images(self):
        self.startTime = np.array([])
        self.exposureTime = np.array([])
        self.ccdTemp = np.array([])
        self.rawList = self.inputList
        self.imageArrays = np.array([])
        self.process_general_information(self.rawList[0])
        self.imageWriter = hcipy.plotting.FFMpegWriter(self.outputFile, framerate = 10)
        for file in self.rawList:
            self.process_specific_information(file)
            # self.processSingleImage()
            self.imageArray = MANGOimage.equalizeHistogram(self.imageArray, self.contrast)
            self.imageArray = MANGOimage.rotateImage(self.imageArray, self.rotationAngle)
            self.plot()
        self.imageWriter.close()

    def plot(self):
        wavelength = {'Green Line': '557.7 nm', 'Red Line': '630 nm'}
        fig, ax = plt.subplots()
        # plt.style.use('dark_background')
        ax.imshow(self.imageArray, cmap='gray')
        ax.annotate('N', xy=(350, 20), color='white')
        ax.annotate('E', xy=(620, 520//2), color='white')
        timestring = dt.datetime.utcfromtimestamp(self.currentTime)
        ax.annotate('Date: ' + str(timestring.date()), xy=(20, 20), color='white')
        ax.annotate('Time: ' + str(timestring.time()), xy=(20, 40), color='white')
        ax.annotate('Temp: ' + str(np.round(self.currentTemp,3)) + 'C.', xy=(500, 40), color='white')
        # wavelength = wave[self.label]
        ax.annotate(wavelength[self.label], xy=(500,60), color='white')
        ax.annotate(str(self.currentET)+'s', xy=(500,20), color='white')

        # plt.annotate(filtstring, xy=(800, 690), color='white', xycoords='figure pixels')
        # plt.annotate(expstring, xy=(800, 660), color='white', xycoords='figure pixels')
        ax.annotate(self.code.upper(), xy=(500,470), color='white')
        plt.axis('off')

        self.imageWriter.add_frame(fig)
        plt.close()
        # make figure object and add frame to writer
        # red - 630nm
        # green - 557.7nm
        # print label
        # print exposure time
        # convert unix to utc


# Use shared functions in MANGOimage.py instead
    # def processSingleImage(self):
    #     # Histogram Equalization to adjust contrast [1%-99%]
    #     numberBins = 10000  # A good balance between time and space complexity, and well as precision
    #     imageArray2D_filt = self.imageArray[50:645, 50:469]
    #     imageArray1D = imageArray2D_filt.flatten()
    #     imageHistogram, bins = np.histogram(imageArray1D, numberBins)
    #     imageHistogram = imageHistogram[1:]
    #     bins = bins[1:]
    #     cdf = np.cumsum(imageHistogram)
    #     self.contrast = 99
    #
    #     # spliced to cut off non-image area
    #     cdf = cdf[0:9996]
    #     # don't equalize non-circle area
    #     max_cdf = max(cdf)
    #     contrast = self.contrast
    #     maxIndex = np.argmin(abs(cdf - contrast / 100 * max_cdf))
    #     minIndex = np.argmin(abs(cdf - (100 - contrast) / 100 * max_cdf))
    #     vmax = float(bins[maxIndex])
    #     vmin = float(bins[minIndex])
    #     lowValueIndices = imageArray1D < vmin
    #     imageArray1D[lowValueIndices] = vmin
    #     highValueIndices = imageArray1D > vmax
    #     imageArray1D[highValueIndices] = vmax
    #     self.imageArray1D_to2D = imageArray1D.reshape(469, 419)
    #     # self.imageArray = imageArray1D.reshape(self.imageArray.shape)
    #     self.imageArray[50:645, 50:469] = self.imageArray1D_to2D
    #
    #
    #     # calibration
    #     # self.rotationAngle = 280
    #     self.imageArray = np.fliplr(skimage.transform.rotate(self.imageArray,
    #                                                          self.rotationAngle, order=3)).astype(float)

    def process_general_information(self, file):
        with h5py.File(file, 'r') as hdf5_file:
            data = hdf5_file['image']
            self.code = data.attrs['station']
            self.site_lat = data.attrs['latitude']
            self.site_lon = data.attrs['longitude']
            self.bin_x = data.attrs['bin_x']
            self.bin_y = data.attrs['bin_y']
            self.label = data.attrs['label']
            data_split = re.findall(r'-(\d+)-', file)[0]
            self.siteDate = dt.datetime.strptime(data_split, '%Y%m%d').date()

    def process_specific_information(self, file):
        '''
        Obtains the following attributes:
        1. start_time
        2. exposure_time
        3. ccd_temp
        :param file: input hdf5 file
        :return: None
        '''
        img = h5py.File(file, 'r')['image']
        self.imageArray = np.array(img)
        self.currentTime = img.attrs['start_time']
        self.startTime = np.append(self.startTime, self.currentTime)
        self.currentET = img.attrs['exposure_time']
        self.exposureTime = np.append(self.exposureTime, self.currentET)
        self.currentTemp = img.attrs['ccd_temp']
        self.ccdTemp = np.append(self.ccdTemp, self.currentTemp)


def parse_args():
    """
    Handle the command line arguments.
    Returns:
    Output of argparse.ArgumentParser.parse_args.
    """

    parser = argparse.ArgumentParser(description='Accepting config, input files and output file'
                                                 'to process MANGOImage.')
    parser.add_argument('-c', '--config', dest='config', type=str,
                        help='Config file containing data locations and image specs.')
    parser.add_argument('-i', '--input', dest='inputs', nargs='+', type=str,
                        help='A list of .hdf5 files to process and store in output file.')
    parser.add_argument('-o', '--output', dest='output', type=str,
                        help='Output file to write processed images to.')
    args = parser.parse_args()

    return args



def main():
    command_line_args = parse_args()
    conf = command_line_args.config
    inputs = command_line_args.inputs
    output = command_line_args.output

    QuickLook(conf, inputs, output)

if __name__ == '__main__':
    main()
