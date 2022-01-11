import warnings
import configparser

from .MANGOImage import MANGOImage
import argparse
import re
import pandas as pd
import datetime as dt
import numpy as np
import h5py

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")

# Need docstrings throughout

# PASSING ARGUMENTS vs class attributes???

class ProcessImage:

    def __init__(self, configFile, inputList, outputFile):
        self.inputList = inputList
        self.outputFile = outputFile

        self.calParams = {}

        self.read_config(configFile)
        self.read_calfile()
        self.process_images()
        self.write_to_hdf5()

    def read_config(self, configFile):
        # read in config file
        config = configparser.ConfigParser()
        config.read(configFile)

        self.siteName = config.get('DEFAULT','SITE_NAME')
        self.calFile = config.get('DEFAULT','CALIBRATION_FILE')
        self.calParams['contrast'] = config.getint('DEFAULT','CONTRAST')

    def read_calfile(self):

        with h5py.File(self.calFile, 'r') as f:
            self.calParams['newIMatrix'] = f['New I array'][:]
            self.calParams['newJMatrix'] = f['New J array'][:]
            self.calParams['backgroundCorrection'] = f['Background Correction Array'][:]
            self.calParams['rotationAngle'] = f['Calibration Angle'][()]
            self.latitude = f['Latitude'][:]
            self.longitude = f['Longitude'][:]
        self.longitude[self.longitude<0.] +=360.


    def process_images(self):
        self.startTime = np.array([])
        self.exposureTime = np.array([])
        self.ccdTemp = np.array([])
        self.rawList = self.inputList
        self.imageArrays = []
        self.get_time_independent_information(self.rawList[0])
        for file in self.rawList:
            print(file)
            with h5py.File(file, 'r') as hdf5_file:
                imageData = hdf5_file['image'][:]
                self.get_time_dependent_information(hdf5_file['image'])
                # picture = MANGOimagehdf5.MANGOimage(hdf5_file['image'], self.calParams)

            image = MANGOImage(imageData)
            image.equalizeHistogram(self.calParams['contrast'])
            # self.showImage()
            # self.setLensFunction()
            image.rotateImage(self.calParams['rotationAngle'])
            # self.showImage()
            image.mercatorUnwrap(self.calParams['newIMatrix'], self.calParams['newJMatrix'], self.calParams['backgroundCorrection'])
            # self.showImage()
            # self.writePNG()
            # self.showImage()

            self.imageArrays.append(image.imageData)

        self.imageArrays = np.array(self.imageArrays)
        self.imageMask = image.alphaMask



    def get_time_independent_information(self, file):
        with h5py.File(file, 'r') as hdf5_file:
            data = hdf5_file['image']
            self.code = data.attrs['station']
            self.site_lat = data.attrs['latitude']
            self.site_lon = data.attrs['longitude']
            self.bin_x = data.attrs['bin_x']
            self.bin_y = data.attrs['bin_y']
            self.label = data.attrs['label']
            # get this from starttime rather than file name?
            data_split = re.findall(r'-(\d+)-', file)[0]
            self.siteDate = dt.datetime.strptime(data_split, '%Y%m%d').date()

    # consider appending this in loop?
    # Just a style preference
    def get_time_dependent_information(self, img):
        '''
        Obtains the following attributes:
        1. start_time
        2. exposure_time
        3. ccd_temp
        :param file: input hdf5 file
        :return: None
        '''

        self.startTime = np.append(self.startTime, img.attrs['start_time'])
        self.exposureTime = np.append(self.exposureTime, img.attrs['start_time'])
        self.ccdTemp = np.append(self.ccdTemp, img.attrs['start_time'])



    def write_to_hdf5(self):

        self.endTime = self.startTime + self.exposureTime
        # tstmp_s = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.startTime])
        # tstmp_e = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.endTime])
        tstmp_s = self.startTime
        tstmp_e = self.endTime

        # save hdf5 file
        with h5py.File(self.outputFile, 'w') as f:
            f.create_group('SiteInfo')

            T = f.create_dataset('UnixTime', data=np.array([tstmp_s, tstmp_e]), compression='gzip', compression_opts=1)
            T.attrs['Description'] = 'unix time stamp'
            T.attrs['Unit'] = 'seconds'

            # 250 km should be read in from somewhere - calibration file?
            # this is associated with the locations of the latitude, so that would make sense
            # Consider just copying array with attached attributes from calibration file?
            Lat = f.create_dataset('Latitude', data=self.latitude, compression='gzip', compression_opts=1)
            Lat.attrs['Description'] = 'geodetic latitude of each pixel projected to 250 km'
            Lat.attrs['Size'] = 'Ipixels x Jpixels'
            Lat.attrs['Projection Altitude'] = 250
            Lat.attrs['Unit'] = 'degrees'

            Lon = f.create_dataset('Longitude', data=self.longitude, compression='gzip', compression_opts=1)
            Lon.attrs['Description'] = 'geodetic longitude of each pixel projected to 250 km'
            Lon.attrs['Size'] = 'Ipixels x Jpixels'
            Lon.attrs['Projection Altitude'] = 250
            Lon.attrs['Unit'] = 'degrees'

            I = f.create_dataset('ImageData', data=self.imageArrays, compression='gzip', compression_opts=1)
            I.attrs['Description'] = 'pixel values for images'
            I.attrs['Site Abbreviation'] = self.code
            I.attrs['Image label'] = self.label
            I.attrs['x/y binning'] = np.array([self.bin_x, self.bin_y])

            M = f.create_dataset('Mask', data=self.imageMask, compression='gzip', compression_opts=1)
            M.attrs['Description'] = 'image mask'

            CCD = f.create_dataset('CCDTemperature', data=self.ccdTemp)
            CCD.attrs['Description'] = 'Temperature of CCD'
            CCD.attrs['Unit'] = 'degrees C'

            N = f.create_dataset('SiteInfo/Name', data=self.siteName)
            N.attrs['Description'] = 'site name'
            C = f.create_dataset('SiteInfo/Code', data=self.code)
            C.attrs['Description'] = 'one letter site abbreviation/code'
            L = f.create_dataset('SiteInfo/Coordinates', data=[self.site_lat, self.site_lon])
            L.attrs['Description'] = 'geodetic coordinates of site; [lat, lon]'
            L.attrs['Unit'] = 'degrees'


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
    ProcessImage(conf, inputs, output)

if __name__ == '__main__':
    main()
