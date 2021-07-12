import MANGOimage
# import create_hdf5_files
import glob
import os
import warnings
import pandas as pd
import configparser

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")


class ProcessImage:

    def __init__(self, site, date, specs=False, check_config=False):
        self.siteName = site
        self.siteDate = date
        self.specs = specs
        self.custom = check_config

        self.preprocessing()
        self.read_and_process()

    def preprocessing(self):
        self.check_config()
        self.build_specs(self.specs, self.config)
        self.directories = self.build_paths(self.custom, self.config)

    def read_and_process(self):
        self.read_in_data(self.directories)
        self.process_image()

    def check_config(self):
        # read in config file
        self.config = configparser.ConfigParser()
        self.config.read('example_config.ini')
        self.specs = self.config['DEFAULT']['any_customization'] == 'Yes'
        self.custom = self.config['DEFAULT']['any_specification'] == 'Yes'

    def build_specs(self, specs, cf):
        if specs:
            self.siteName = cf['Specifications']['siteName']
            self.siteDate = cf['Specifications']['siteDate']

    def build_paths(self, custom, cf):
        list_of_dirs = ['parent', 'rawData', 'rawSite', 'rawSiteFiles', 'rawImages', 'processedImages']
        paths = {}
        rawFolder = "raw_data"
        siteFiles = "site_files"
        siteName = self.siteName
        siteDate = self.siteDate
        if not custom:
            processedFolder = os.path.join('processed_data', 'processed_images')
            parentFolder = os.path.dirname(os.getcwd())
            dirs_to_be_assigned = [parentFolder,
                                   os.path.join(parentFolder, rawFolder),
                                   os.path.join(parentFolder, rawFolder, siteName),
                                   os.path.join(parentFolder, rawFolder, siteName, siteFiles),
                                   os.path.join(parentFolder, rawFolder, siteName, siteDate),
                                   os.path.join(parentFolder, processedFolder, siteName, siteDate)]
            paths = dict(zip(list_of_dirs, dirs_to_be_assigned))
        else:
            for i in list_of_dirs:
                paths[i] = cf['Data Locations'][i]

        return paths

    def read_in_data(self, directories):
        siteInfoFile = "SiteInformation.csv"
        siteInfoFilePath = os.path.join(directories['rawData'], siteInfoFile)
        siteData = pd.read_csv(siteInfoFilePath).to_numpy()

        # raw data --> one site --> images for one day
        rawImagesPath = directories['rawImages']
        # rawSitePath = directories['rawSite']
        # list of raw images
        self.rawList = glob.glob1(rawImagesPath, '*.[0-9]*')
        # OR rawList= [file for file in os.listdir(rawPath) if file.endswith('.png')]

        # processed images location
        self.processedImagesFolder = directories['processedImages']

        if not os.path.isdir(self.processedImagesFolder):
            os.makedirs(self.processedImagesFolder)

        # process individual images

    def process_image(self):
        for rawImage in self.rawList:
            writeAddress = os.path.join(self.processedImagesFolder, rawImage[0:8] + ".png")
            if not os.path.exists(writeAddress):
                MANGOimage.MANGOimage(self.directories, rawImage)


def main():
    siteName = "CFS"
    date = "20210305"
    ProcessImage(siteName, date)

# if __name__ == '__main__':
#    main()
