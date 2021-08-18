import sys

import MANGOimage
# import create_hdf5_files
import glob
import os
import warnings
import configparser

import MANGOimagehdf5
import writeHDF5
import sys
import argparse

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")


class ProcessImage:

    def __init__(self, configFile, inputList, outputFile):
        # self.siteName = site
        # self.siteDate = date
        # self.specs = specs
        # self.custom_paths = check_config
        self.configFile = configFile
        self.inputList = inputList
        self.outputFile = outputFile
        self.setup()
        self.read_and_process()

    def setup(self):
        self.check_config()
        self.build_specs()
        self.directories = self.build_paths()

    def read_and_process(self):
        # self.read_in_data()
        self.process_images()
        self.write_to_hdf5()

    def check_config(self):
        # read in config file
        self.config = configparser.ConfigParser()
        self.config.read(self.configFile)
        self.directories = self.config['Data Locations']
        self.specs = self.config['DEFAULT']['any_specification'] == 'Yes'
        self.custom_paths = self.config['DEFAULT']['any_path_customization'] == 'Yes'

    def build_specs(self):
        # if self.specs:
        self.siteName = self.config['Specifications']['siteName']
        self.siteDate = self.config['Specifications']['siteDate']

    def build_paths(self):
        paths = {}
        '''rawFolder = "raw_data"
        siteFiles = "site_files"
        siteName = self.siteName
        siteDate = self.siteDate

        # processedFolder = os.path.join('processed_data', 'processed_images')
        parentFolder = os.path.dirname(os.getcwd())

        paths['parent'] = parentFolder
        paths['rawData'] = os.path.join(parentFolder, rawFolder)
        paths['rawSite'] = os.path.join(parentFolder, rawFolder, siteName)
        paths['rawSiteFiles'] = os.path.join(parentFolder, rawFolder, siteName, siteFiles)

        if self.siteDate is None:
            site_folder = next(os.walk(paths['rawSite']))[1]
            site_folder.remove('site_files')
            for siteDate in site_folder:
                ProcessImage(self.siteName, siteDate)
            sys.exit('All dates for specified site have been processed.')

        paths['rawImages'] = os.path.join(parentFolder, rawFolder, siteName, siteDate)
        paths['processedImages'] = os.path.join(parentFolder, processedFolder, siteName, siteDate)

        # if self.custom_paths:
        
        for i in paths.keys():
            if self.config['Data Locations'][i] != '':
                paths[i] = self.config['Data Locations'][i]
        '''

        return

    def read_in_data(self):
        # raw data --> one site --> images for one day
        # rawImagesPath = self.directories['rawImages']
        # rawSitePath = directories['rawSite']
        # list of raw images

        # processed images location
        # self.processedImagesFolder = self.directories['processedImages']

        # if not os.path.isdir(self.processedImagesFolder):
        # os.makedirs(self.processedImagesFolder)

        # process individual images
        return

    def process_images(self):
        self.imageArrays = []
        self.rawList = self.inputList
        for file in self.rawList:
            # writeAddress = os.path.join(self.processedImagesFolder, rawImage[0:-5] + ".png")
            MANGOimagehdf5.MANGOimage(self.directories, file, self.config, self.imageArrays)

    def write_to_hdf5(self):
        datadict = {}
        siteInfoFile = 'SiteInformation.csv'
        datadict['pathToSiteFile'] = os.path.join(self.directories['rawData'], siteInfoFile)
        datadict['pathToLatLon'] = os.path.join(self.directories['rawSiteFiles'], 'calibration')
        # datadict['pathToProcessedSite'] = os.path.dirname(self.processedImagesFolder)
        datadict['inputList'] = self.inputList
        datadict['outputFile'] = self.outputFile
        datadict['imageArrays'] = self.imageArrays
        writeHDF5.hdf5_file_info(datadict, self.config)


def parse_args():
    """
    Handle the command line arguments.
    Returns:
    Output of argparse.ArgumentParser.parse_args.
    """

    parser = argparse.ArgumentParser(description='Accepting config, input files and output file'
                                                 'to process MANGOImage.')
    parser.add_argument('-c', '--config', dest='config', type=argparse.FileType('r'),
                        help='Config file containing data locations and image specs.')
    parser.add_argument('-i', '--input', dest='inputs', nargs='+', type=argparse.FileType('r'),
                        help='A list of .hdf5 files to process and store in output file.')
    parser.add_argument('-o', '--output', dest='output', type=argparse.FileType('w+'),
                        help='Output file to write processed images to.')
    args = parser.parse_args()

    '''check = lambda filename: filename.lower().endswith('hdf5')
    if not all(check(input_file.name) for input_file in args.inputs):
        sys.stderr.write('All inputs must be .hdf5')
        sys.exit(1)'''

    return args


if __name__ == '__main__':
    command_line_args = parse_args()
    conf = command_line_args.config
    inputs = command_line_args.inputs
    output = command_line_args.output
    print(conf, inputs, output)
    ProcessImage(conf, inputs, output)
