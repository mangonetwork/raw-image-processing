import sys

import MANGOimage
# import create_hdf5_files
import glob
import os
import warnings
import pandas as pd
import configparser
import create_hdf5_files
import sys

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")


class ProcessImage:

    def __init__(self, site, date=None, specs=False, check_config=False):
        self.siteName = site
        self.siteDate = date
        self.specs = specs
        self.custom = check_config

        self.preprocessing()
        self.read_and_process()

    def preprocessing(self):
        self.check_config()
        self.build_specs(self.config)
        self.directories = self.build_paths(self.custom, self.config)

    def read_and_process(self):
        self.read_in_data()
        self.process_image()
        self.write_to_hdf5()

    def check_config(self):
        # read in config file
        self.config = configparser.ConfigParser()
        self.config.read('example_config.ini')
        self.specs = self.config['DEFAULT']['any_customization'] == 'Yes'
        self.custom = self.config['DEFAULT']['any_specification'] == 'Yes'

    def build_specs(self, cf):
        if self.specs:
            self.siteName = cf['Specifications']['siteName']
            self.siteDate = cf['Specifications']['siteDate']

    def build_paths(self, custom, cf):
        paths = {}
        rawFolder = "raw_data"
        siteFiles = "site_files"
        siteName = self.siteName
        siteDate = self.siteDate
        if not custom:
            processedFolder = os.path.join('processed_data', 'processed_images')
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

        for i in paths.keys():
            if cf['Data Locations'][i] != '':
                paths[i] = cf['Data Locations'][i]

        return paths

    def read_in_data(self):
        # raw data --> one site --> images for one day
        rawImagesPath = self.directories['rawImages']
        # rawSitePath = directories['rawSite']
        # list of raw images
        self.rawList = glob.glob(os.path.join(rawImagesPath, '*.[0-9]*'))
        # OR rawList= [file for file in os.listdir(rawPath) if file.endswith('.png')]

        # processed images location
        self.processedImagesFolder = self.directories['processedImages']

        if not os.path.isdir(self.processedImagesFolder):
            os.makedirs(self.processedImagesFolder)

        # process individual images

    def process_image(self):
        for rawImage in self.rawList:
            rawImage = os.path.basename(rawImage)
            writeAddress = os.path.join(self.processedImagesFolder, rawImage[0:8] + ".png")
            if not os.path.exists(writeAddress):
                img = MANGOimage.MANGOimage(self.directories, rawImage)
                img.load_files()
                img.process()

    def write_to_hdf5(self):
        siteInfoFile = "SiteInformation.csv"
        pathToSiteFile = os.path.join(self.directories['rawData'], siteInfoFile)
        pathToLatLon = os.path.join(self.directories['rawSiteFiles'], 'calibration')
        pathToProcessedSite = os.path.dirname(self.processedImagesFolder)
        create_hdf5_files.hdf5_file_info(pathToSiteFile, pathToProcessedSite, pathToLatLon, self.siteName, self.siteDate)


def date_formatter(originalName):
    dates = {
        'Jan': '01',
        'Feb': '02',
        'Mar': '03',
        'Apr': '04',
        'May': '05',
        'Jun': '06',
        'Jul': '07',
        'Aug': '08',
        'Sep': '09',
        'Oct': '10',
        'Nov': '11',
        'Dec': '12'
    }
    month = dates[originalName[:3]]
    day = originalName[3:5]
    year = '20' + originalName[5:]
    converted_str = year + month + day
    return converted_str


if __name__ == '__main__':

    siteName = "EIO"
    date = "20210305"

    ProcessImage(siteName)
