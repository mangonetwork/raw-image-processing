# process_single_day.py


import MANGOimage
import create_hdf5_files
from pathlib import Path
import glob
import io
import numpy
import os
import warnings
import configparser
warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")
# Read the site information file to get 'PixArray'
siteInfofname = "SiteInformation.csv"
f = open(siteInfofname,'r')
siteText = f.read()
s = io.StringIO(siteText)
siteData = numpy.genfromtxt(s, dtype="|S", delimiter = ',', autostrip=True)

rows = siteData[:,0]

PixArray = []
for i in range(1,len(rows)):
    SiteInfo = siteData[i,:]
    PixArray.append(SiteInfo)
    
siteName = "Capitol Reef Field Station"
siteDir = "May 8 2016 data"
rawFolder = ""
parentFolder = "C:/Users/padma/MANGO SU21/"
#directory where unprocessed data are
rawPath = parentFolder + siteDir + "/" + rawFolder

#list of raw images
rawList = glob.glob1(rawPath, '*.129')
#OR rawList= [file for file in os.listdir(rawPath) if file.endswith('.png')]

# make directory based on folder being processed in the appropriate location
processedFolder =  "Processed data/"
pfoldersdir = rawPath + processedFolder


if not os.path.isdir(pfoldersdir):
    os.makedirs(pfoldersdir)
    
    
# process individual images
for rawImage in rawList:
    writeAddress = pfoldersdir + "/" + rawImage[0:8] + ".png"
    if not os.path.exists(writeAddress):
        MANGOimage.MANGOimage(siteName, siteDir, rawFolder, rawPath, rawImage, PixArray, pfoldersdir)

create_hdf5_files.hdf5_file_info('SiteInformation.csv', pfoldersdir, parentFolder)