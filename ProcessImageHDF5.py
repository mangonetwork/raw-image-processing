import warnings
import configparser

import MANGOimagehdf5
import argparse
import re
import pandas as pd
import datetime as dt
import numpy as np
import h5py

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")

# Need docstrings throughout

class ProcessImage:

    def __init__(self, configFile, inputList, outputFile):
        self.configFile = configFile
        self.inputList = inputList
        self.outputFile = outputFile
        self.build_config()
        self.process_images()
        self.write_to_hdf5()

    def build_config(self):
        # read in config file
        self.config = configparser.ConfigParser()
        self.config.read(self.configFile)
        # self.directories = self.config['Data Locations']
        # is only one value needed from this config file?  Can probably do away with it

    def process_images(self):
        self.startTime = np.array([])
        self.exposureTime = np.array([])
        self.ccdTemp = np.array([])
        self.rawList = self.inputList
        self.imageArrays = np.array([])
        self.process_general_information(self.rawList[0])
        for file in self.rawList:
            self.process_specific_information(file)
            picture = MANGOimagehdf5.MANGOimage(file, self.config, self.imageArrays)
            self.imageArrays = picture.allImagesArray

        #self.imageArrays = np.array(self.imageDict.values(), dtype=np.ndarray)

# potentially combine these two functions?
# depends on overhead of opening single hdf5 file multiple times
    def process_general_information(self, file):
        # with construct here as well
        hdf5_file = h5py.File(file, 'r')
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
        # Use with construct here?  Avoid leaving files opeen
        img = h5py.File(file, 'r')['image']
        self.startTime = np.append(self.startTime, img.attrs['start_time'])
        self.exposureTime = np.append(self.exposureTime, img.attrs['start_time'])
        self.ccdTemp = np.append(self.ccdTemp, img.attrs['start_time'])

    def write_to_hdf5(self):
        sitefile = self.config['Data Locations']['siteInfoFile']    # Not used?
        # create site list from the site file and user input
        self.site_name = 'Capitol Reef Field Station'   # temporary hard-coding - let's make sure to make a note of this

        # Move this to a different function
        # this function should ONLY be writting the output hdf5
        # read lat/lon from where ever Latitude.csv and Longitude.csv are for that site
        latDir = self.config['Data Locations']['latitudeFile']
        lonDir = self.config['Data Locations']['longitudeFile']
        try:
            latitude = np.array(pd.read_csv(latDir, dtype=float, delimiter=','))
            longitude = np.array(pd.read_csv(lonDir, dtype=float, delimiter=','))
            longitude[longitude < 0] += 360.
        except IOError:
            print('Could not process {}!'.format(self.site_name))

        self.endTime = self.startTime + self.exposureTime
        # tstmp_s = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.startTime])
        # tstmp_e = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.endTime])
        tstmp_s = self.startTime
        tstmp_e = self.endTime

        # save hdf5 file
        # with statement here
        f = h5py.File(self.outputFile, 'w')
        f.create_group('SiteInfo')

        T = f.create_dataset('UnixTime', data=np.array([tstmp_s, tstmp_e]), compression='gzip', compression_opts=1)
        T.attrs['Description'] = 'unix time stamp'
        T.attrs['Unit'] = 'seconds'

        # 250 km should be read in from somewhere - calibration file?
        # this is associated with the locations of the latitude, so that would make sense
        # Consider just copying array with attached attributes from calibration file?
        Lat = f.create_dataset('Latitude', data=latitude, compression='gzip', compression_opts=1)
        Lat.attrs['Description'] = 'geodetic latitude of each pixel projected to 250 km'
        Lat.attrs['Size'] = 'Ipixels x Jpixels'
        Lat.attrs['Projection Altitude'] = 250
        Lat.attrs['Unit'] = 'degrees'

        Lon = f.create_dataset('Longitude', data=longitude, compression='gzip', compression_opts=1)
        Lon.attrs['Description'] = 'geodetic longitude of each pixel projected to 250 km'
        Lon.attrs['Size'] = 'Ipixels x Jpixels'
        Lon.attrs['Projection Altitude'] = 250
        Lon.attrs['Unit'] = 'degrees'

        I = f.create_dataset('ImageData', data=self.imageArrays, compression='gzip', compression_opts=1)
        I.attrs['Description'] = 'pixel values for images'
        I.attrs['Site Abbreviation'] = self.code
        I.attrs['Image label'] = self.label
        I.attrs['x/y binning'] = np.array([self.bin_x, self.bin_y])

        CCD = f.create_dataset('CCDTemperature', data=self.ccdTemp)
        CCD.attrs['Description'] = 'Temperature of CCD'
        CCD.attrs['Unit'] = 'degrees C'

        N = f.create_dataset('SiteInfo/Name', data=self.site_name)
        N.attrs['Description'] = 'site name'
        C = f.create_dataset('SiteInfo/Code', data=self.code)
        C.attrs['Description'] = 'one letter site abbreviation/code'
        L = f.create_dataset('SiteInfo/Coordinates', data=[self.site_lat, self.site_lon])
        L.attrs['Description'] = 'geodetic coordinates of site; [lat, lon]'
        L.attrs['Unit'] = 'degrees'

        f.close()


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


if __name__ == '__main__':
    command_line_args = parse_args()
    conf = command_line_args.config
    inputs = command_line_args.inputs
    output = command_line_args.output
    ProcessImage(conf, inputs, output)
