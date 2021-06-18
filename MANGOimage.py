#!/usr/bin/env python

#######################################################################
#
#   Image Processing Core Engine for MANGO
#   (Midlatitude All-sky-imager Network for Geophysical Observation)
#   2013-01-23  Rohit Mundra
#               Initial implementation
#
#   2013-02-05  Rohit Mundra
#               Robust Star Removal implemented
#               Histogram Equalization implemented
#
#   2013-02-22  Rohit Mundra
#               Changing unwarping to map to Normal Mercator Projection
#               Added alpha mask and PIL PNG saving
#               
#   2013-02-22  Fabrice Kengne
#               Adding compatibility for FITS read
#
#   2014-10-25  Asti Bhatt
#		Adding compatibility for BU software generated .153 files
#		
#######################################################################


from numpy import genfromtxt, linalg
from math import atan, acos, asin, atan2, sin, cos, tan, radians, isnan

from scipy import misc
from scipy.interpolate import griddata

import numpy as np
import PIL
import matplotlib
matplotlib.use('TkAgg') # to get around tkinter issue
#import pyfits
import copy
import string
import os
import glob
import sys
import logging
import matplotlib.pyplot as plt
import skimage.transform
import pandas as pd

DIRN = "C:/Users/padma/MANGO SU21/"
#DIRN = "~/venv/SERVER/myproject/MANGO/app/static/Sites/"
#DIRN = '/home/abhatt/workspace/WPI-SRI MQP - MANGO SYSTEM1/SERVER/myproject/MANGO/app/static/Sites/'
#DIRN = "C:\\Users\\WPI\\Documents\\WPI-SRI MQP - MANGO SYSTEM1\\SERVER\\myproject\\MANGO\\app\\static\\Sites\\"
#print(DIRN)
class MANGOimage:
    def __init__(self, siteName, siteDir, rawFolder, rawPath, rawImage, PixArray, pfoldersdir):
        self.rawPath = rawPath
        self.siteDir = siteDir
        self.rawFolder = rawFolder
        self.rawImage = rawImage
        self.PixArray = PixArray
        self.rawImageAddress = self.rawPath + self.rawImage 
        self.siteName = siteName
        self.pfoldersdir = pfoldersdir
        try:
            self.loadFITS()
        except ValueError:
            return
        self.process()
	
    def loadFITS(self):
    	# commented version below is for pyfits file.
            #hdulist                 = pyfits.open(self.rawImageAddress)
            #self.imageData          = flipud(array((hdulist[0].data)))
            #self.width              = int(hdulist[0].header['NAXIS1']) #695
            #self.height             = int(hdulist[0].header['NAXIS2']) #519
            #self.siteName           = hdulist[0].header['SITE']
            #self.calibrationFlag    = 0
            #self.writeMode          = 'L'
        # for BU sw generated binary data files. 	
        f = open(self.rawImageAddress, "rb")
        a = np.fromfile(f, dtype='int16')
        #print(a.shape)
        self.width = 695
        self.height = 519
        self.imageData = a[64:360769].reshape([519, 695]).astype('int32')
    	
    #	print self.imageData.shape 
    #	figure1 = self.imageData
    #	plt.imshow(figure1)
    #	plt.title('Raw Image')
    #	plt.gray()
    #	plt.show()
	
    def getSiteName(self):
        for site in self.PixArray:
            if site[1] == self.rawImageAddress[0]:
                self.siteName = site[0]
                

    def loadCalibrationData(self):
        calibrationFilename     = DIRN + self.siteName + "/calibration/" + 'Calibration.csv' # uncomment for windows
#	calibrationFilename     = os.path.dirname(os.getcwd()) + "/Sites/" + self.siteName + "/calibration/" + 'Calibration.csv' # uncomment for linux

        ## This is concatenated into the SiteInfo in getSiteInformation
        #f = open(calibrationFilename,'r')
        #text = f.read()
        #s = StringIO(text)
        #calibrationData         = np.genfromtxt(s, delimiter=',')       
        calibrationData         = genfromtxt(calibrationFilename, delimiter=',')
        self.azimuth            = np.array((calibrationData[1:,1]))
        self.elevation          = np.array((calibrationData[1:,2]))
        self.i                  = np.array((calibrationData[1:,3]))
        self.j                  = np.array((calibrationData[1:,4]))
        self.zenith             = np.array([self.i[0], self.j[0]])
#	print self.zenith
        
    def loadNewIJ(self):
        newIFilename            = DIRN + self.siteName + "/calibration/" + 'newI.csv' # uncomment for windows 
        newJFilename            = DIRN + self.siteName + "/calibration/" + 'newJ.csv' # uncomment for windows
#	newIFilename            = os.path.dirname(os.getcwd()) + "/Sites/" + self.siteName + "/calibration/" + 'newI.csv' # uncomment for linux
#        newJFilename           = os.path.dirname(os.getcwd()) + "/Sites/" + self.siteName + "/calibration/" + 'newJ.csv' # uncomment for linux
        nif_df = pd.read_csv(newIFilename, delimiter = ',', header = None)
        self.newIMatrix = np.array(nif_df)
        
        njf_df = pd.read_csv(newJFilename, delimiter = ',', header = None)
        self.newJMatrix = np.array(njf_df)
        
        #old code:
        #self.newIMatrix = genfromtxt(newIFilename, delimiter = ',')
        #self.newJMatrix = genfromtxt(newJFilename, delimiter = ',')

    def loadBackgroundCorrection(self):
        backgroundCorrectionFilename = DIRN + self.siteName + "/calibration/" + 'backgroundCorrection.csv' # uncomment for windows
        bg_corr_df = pd.read_csv(backgroundCorrectionFilename, delimiter = ',', header = None)
        #backgroundCorrectionFilename    = os.path.dirname(os.getcwd()) + "/Sites/" + self.siteName + "/calibration/" + 'backgroundCorrection.csv' # uncomment for linux
        self.backgroundCorrection = np.array(bg_corr_df)
        #self.backgroundCorrection = genfromtxt(backgroundCorrectionFilename, delimiter = ',')

    def process(self):
        #self.getSiteName()
        self.loadCalibrationData()
        self.loadNewIJ()
        self.loadBackgroundCorrection()
        self.equalizeHistogram()
        #self.removeStars()
        self.setLensFunction()
        self.calibrate()
        self.mercatorUnwarp()
        #self.equalizeHistogram()
        self.writePNG()

    def removeStars(self):
        filteredData = copy.copy(self.imageData).astype(float)
        for i in range(self.width):
            leftPixels = self.imageData[:,max(i-5, 0):max(i-2,1)]
            rightPixels = self.imageData[:,min(i+3, self.width-2):min(i+6, self.width-1)]
            allPixels = np.append(leftPixels, rightPixels, 1)
            pixelData = self.imageData[:, i]
            stdVal = np.std(allPixels, 1)
            meanVal = np.mean(allPixels, 1)
            for j in range(self.height):
                if pixelData[j]>meanVal[j]+3*stdVal[j]:
                    filteredData[j-1:j+2, i-1:i+2] = -1
                (iKnown, jKnown) = np.where(filteredData>=0)
            valuesKnown = np.extract(filteredData>=0,filteredData)
        #	find minimum value from the imageData to find threshold
        m = np.empty(self.width)
        for i in range(self.width):
            m[i] = min(self.imageData[:,i])
        threshold = min(m)
    #   (iAll, jAll) = where(filteredData>=-1) # comment for .153 files 
        (iAll, jAll) = np.where(filteredData>=threshold) # min value from imageData
    #	(iAll, jAll) = where(filteredData>=-32725)
        self.imageData = np.reshape(griddata((iKnown, jKnown), valuesKnown, (iAll, jAll), method='linear', fill_value = 0), filteredData.shape)
#       figure2 = self.imageData
#       plt.imshow(figure2)
#       plt.title('Star Removal Image')
#       plt.gray()
#       plt.show()
    
    def calibrate(self):
        #Spatial Calibration based on star data
        self.zenith = np.array([self.i[0], self.j[0]])
        G_el = 1.0 - [self.getPixelsFromAngle(angle) for angle in self.elevation]/self.fisheyeRadius
        self.f = G_el*np.sin(np.radians(self.azimuth))
        self.g = G_el*np.cos(np.radians(self.azimuth))
        firstColumn = np.array(([1]*len(self.i)))
        oneIJ = np.vstack((firstColumn, self.i, self.j)).transpose()

        intermediate0 = np.dot(oneIJ.transpose(), oneIJ)
        intermediate1 = linalg.pinv(intermediate0)
        intermediate2 = np.dot(intermediate1, oneIJ.transpose())

        aCoefficients = np.dot(intermediate2, self.f)
        bCoefficients = np.dot(intermediate2, self.g)
        
        self.a0 = aCoefficients[0]
        self.a1 = aCoefficients[1]
        self.a2 = aCoefficients[2]
        self.b0 = bCoefficients[0]
        self.b1 = bCoefficients[1]
        self.b2 = bCoefficients[2]

        rotationAngle_1 = np.degrees(atan(-self.b1/self.a1))
        rotationAngle_2 = np.degrees(atan(self.a2/self.b2))
        self.rotationAngle = .5*(rotationAngle_1 + rotationAngle_2)
        #if self.siteName == 'Eastern Iowa Observatory':
            #self.rotationAngle = .5*(rotationAngle_1 + rotationAngle_2)+180
            #print self.rotationAngle
        self.imageData = np.fliplr(skimage.transform.rotate(self.imageData, self.rotationAngle, order = 3)).astype(float)
        #figure3 = self.imageData
        #plt.imshow(figure3)
        #plt.title('Calibrated Image')
        #plt.gray()
        #plt.show()
        
        #Rotating Zenith Counter-clockwise by rotation angle and flipping it left-right
        zenithI = self.width - int(cos(radians(self.rotationAngle))*(self.zenith[0] - self.width/2) - sin(radians(self.rotationAngle))*(self.zenith[1]-self.height/2) + self.width/2)
        zenithJ = int(sin(radians(self.rotationAngle))*(self.zenith[0] - self.width/2) + cos(radians(self.rotationAngle))*(self.zenith[1]-self.height/2) + self.height/2)
        self.zenith = [zenithI, zenithJ]
        self.setLensFunction()

    def setLensFunction(self):
        #Calculates lens function coefficients for the equation: Angle = a0 + a1.px + a2.px^2 + a3.px^3
        xDiff = [self.i[0] - i for i in self.i]
        yDiff = [self.j[0] - j for j in self.j]
        distanceFromZenith = [np.sqrt(xDiff[k]**2 + yDiff[k]**2) for k in range(len(self.i))]
        angleFromZenith = [self.elevation[0] - self.elevation[k] for k in range(len(self.i))]
        
        firstColumn = np.array(([1]*len(distanceFromZenith)))
        distanceFromZenithMat = np.vstack((firstColumn, distanceFromZenith, [x**2 for x in distanceFromZenith], [x**3 for x in distanceFromZenith])).transpose()
        self.pixToAngleCoefficients = np.dot(linalg.pinv(distanceFromZenithMat), angleFromZenith)[::-1]

        firstColumn = np.array(([1]*len(distanceFromZenith)))
        angleFromZenithMat = np.vstack((firstColumn, angleFromZenith, [x**2 for x in angleFromZenith], [x**3 for x in angleFromZenith])).transpose()
        self.angleToPixCoefficients = np.dot(linalg.pinv(angleFromZenithMat), distanceFromZenith)[::-1]
        
        self.fisheyeRadius = self.getPixelsFromAngle(90)

    def getPixelsFromAngle(self, angle):
        #Input Angle in degrees
        return np.polyval(self.angleToPixCoefficients, angle)

    def mercatorUnwarp(self):
        newImageWidth = 500
        newImageHeight = 500
        finalImage = np.ones([newImageHeight, newImageWidth])*-1

        for j in range(self.imageData.shape[0]):
            for i in range(self.imageData.shape[1]):
                newI = self.newIMatrix[j][i]
                newJ = self.newJMatrix[j][i]
                if (not(isnan(newI)) and not(isnan(newJ)) and not(isnan(self.backgroundCorrection[j, i]))):
                    self.imageData[j, i] = self.imageData[j, i]*self.backgroundCorrection[j, i]
                    finalImage[int(newJ), int(newI)] = self.imageData[j, i]

        (iKnown, jKnown) = np.where(finalImage>=0)
        valuesKnown = np.extract(finalImage>=0,finalImage)
        (iAll, jAll) = np.where(finalImage>=-1)
        interpolatedData = np.reshape(griddata((iKnown, jKnown), valuesKnown, (iAll, jAll), method='cubic', fill_value = -1), [newImageWidth, newImageHeight])
        alphaMask = np.ones([newImageHeight, newImageWidth])*255
        for j in range(interpolatedData.shape[0]):
            for i in range(interpolatedData.shape[1]):
                if interpolatedData[j, i] == -1:
#                    alphaMask[j, i] = 0
                    alphaMask[j, i] = np.nan # to create transparent background

        interpolatedData = (interpolatedData*255/(np.nanmax(interpolatedData))).astype('uint8')
        alphaMask = alphaMask.astype('uint8')
        self.imageData = np.dstack([interpolatedData, alphaMask])
        self.writeMode = 'LA'

    def equalizeHistogram(self):
        #Histogram Equalization to adjust contrast [1%-99%]
        numberBins = 10000 #A good balance between time and space complexity, and well as precision
        flattenedImageData = self.imageData.flatten()
        imageHistogram, bins = np.histogram(flattenedImageData, numberBins)
        imageHistogram = imageHistogram[1:]
        bins = bins[1:]
        cdf = np.cumsum(imageHistogram)
        cdf = cdf[0:9996]
        max_cdf = max(cdf)
        maxIndex = np.argmin(abs(cdf - 0.95*max_cdf))
        minIndex = np.argmin(abs(cdf - 0.05*max_cdf))
        #the following two lines took forever
        #maxIndex = min(range(len(cdf)), key=lambda i: abs(cdf[i]-0.95*max_cdf))
        #minIndex = min(range(len(cdf)), key=lambda i: abs(cdf[i]-0.05*max_cdf))
        
        #print(maxIndex1, minIndex1, maxIndex, minIndex)
    #	maxI = where(cdf <= 0.99*max(cdf))
    #	self.maxIndex = maxI[-1]
    #	minI = where(cdf <= 0.01*max(cdf))
    #	self.minIndex = minI[-1]
        vmax = float(bins[maxIndex])
        vmin = float(bins[minIndex])
        lowValueIndices = flattenedImageData < vmin
        flattenedImageData[lowValueIndices] = vmin
        highValueIndices = flattenedImageData > vmax

        flattenedImageData[highValueIndices] = vmax
        self.imageData = flattenedImageData.reshape(self.imageData.shape)
        #print(self.imageData.shape)
        return self.imageData

    def writePNG(self):
        #Writing for web display - change output directory to something local
        #writeAddress = DIRN + "All Sites Images\\" + self.siteName + ".png" # uncomment for windows
        #print("Saving image: " + self.rawImage[0:8])
        #processeddirprefix = "/media/abhatt/Seagate Backup Plus Drive/workspace/Data/MANGO/InGeO/"
        #processeddirprefix = "/home/ftp/pub/earthcube/provider/asti/MANGOProcessed/"
        #writeAddress = self.rawPath + "Processed/" + self.rawImage[0:8] + ".png"
        #writeAddress = processeddirprefix + self.siteDir + "/" + self.rawFolder + "/"+ self.rawImage[0:8] + ".png"
        writeAddress = self.pfoldersdir + "/" + self.rawImage[0:8] + ".png" 
        #change write addres
        #writeAddress = os.path.dirname(os.getcwd()) + "/Sites/All Sites Images/" + self.siteName + ".png" # uncomment for linux
        #writeAddress = self.rawImageAddress[-12:-4] + ".png"
        #writeAddress = self.rawImageAddress[0:-13] + "/Processed1/" + self.rawImageAddress[-12:-4] + ".png"
        finalImage = PIL.Image.fromarray(self.imageData, self.writeMode)
        #implt =	plt.imshow(finalImage, clim=(0,200))
        #plt.imshow(finalImage)
        #plt.axis('off')
        #implt.set_clim(self.minIndex, self.maxIndex)
        #implt.set_clim(0, 250)
        #plt.savefig(writeAddress, transparent=True, frameon=False, bbox_inches='tight', pad_inches=0)
        finalImage.save(writeAddress, format = 'png')
        #plt.clf() 
        #Writing for archiving to source folder
        #writeAddress = self.rawImageAddress[:-3] + "png"
        #finalImage.save(writeAddress)
        #plt.imshow(finalImage)
        #plt.title('Final Image')
        #plt.show()
        