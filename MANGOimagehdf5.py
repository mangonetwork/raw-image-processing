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
# 2014-10-25  Asti Bhatt
# Adding compatibility for BU software generated .153 files
#######################################################################


from scipy.interpolate import griddata

import numpy as np
import PIL.Image
import os
import copy
import skimage.transform
import pandas as pd
import scipy.linalg as slg
import h5py
from PIL import Image


class MANGOimage:
    def __init__(self, rawimg, config, collective_img_arr):
        # list_of_dirs = ['parent', 'rawData', 'rawSite', 'rawSiteFiles', 'rawImages', 'processedImages']
        self.rawImage = rawimg
        self.config = config
        # self.rawSiteFilesPath = self.config['Data Locations']['rawSiteFiles']
        self.cal_hdf = self.config['Data Locations']['cal_hdf']
        self.allImagesArray = collective_img_arr

        self.contrast = int(config['Specifications']['contrast'])

        # is this nessisary, or do we want to just let the original error appear
        # these statements create a clean interface, but sometimes obtuse what the actual problem is
        try:
            self.loadFITS()
            self.load_files()
            self.process()
        except ValueError:
            raise ValueError('Raw image file cannot be processed.')

    # rename function?
    # FITS are a very specific type of image file used by some (non-MANGO) cameras
    def loadFITS(self):
        data = h5py.File(self.rawImage, 'r')
        self.imageData = data['image']
        self.imageArray = np.array(data['image'])
        self.width = self.imageData.attrs['width']
        self.height = self.imageData.attrs['height']
        self.imageArray1D = self.imageArray.flatten()

    def load_files(self):
        # self.loadCalibrationData()
        self.loadNewIJ()
        self.loadBackgroundCorrection()

    def process(self):
        self.equalizeHistogram()
        # self.showImage()
        # self.setLensFunction()
        self.calibrate()
        # self.showImage()
        self.mercatorUnwrap(self.imageData)
        # self.showImage()
        self.writePNG()
        # self.showImage()

    # function for debugging?
    # fine to keep this, but make a note of it
    def showImage(self):
        img = Image.fromarray(self.imageData)
        img.show()

    # # change with new calibration format
    # def loadCalibrationData(self):
    #     # these values should be in the calibration hdf5 file
    #     calibrationFilename = self.config['Data Locations']['calibrationFile']
    #     calibrationData = pd.read_csv(calibrationFilename, delimiter=',', index_col='Star Name')
    #     self.azimuth = np.array(calibrationData['Azimuth'])
    #     self.elevation = np.array(calibrationData['Elevation'])
    #     self.i = np.array(calibrationData['i Coordinate'])
    #     self.j = np.array(calibrationData['j Coordinate'])
    #     self.zenith = np.array(calibrationData.iloc[0].loc[['i Coordinate', 'j Coordinate']])
    #     self.zenithI = self.zenith[0]
    #     self.zenithJ = self.zenith[1]

    # change with new calibration format
    def loadNewIJ(self):
        # newIFilename = os.path.join(self.rawSiteFilesPath, 'newI.csv')
        # newJFilename = os.path.join(self.rawSiteFilesPath, 'newJ.csv')
        # nif_df = pd.read_csv(newIFilename, delimiter=',', header=None)
        # self.newIMatrix = np.array(nif_df)
        # njf_df = pd.read_csv(newJFilename, delimiter=',', header=None)
        # self.newJMatrix = np.array(njf_df)
        # print('NewI, NewJ')
        # print(self.newIMatrix.shape, self.newJMatrix.shape)
        with h5py.File(self.cal_hdf, 'r') as f:
            self.newIMatrix = f['New I array'][:]
            self.newJMatrix = f['New J array'][:]
        # print(self.newIMatrix.shape, self.newJMatrix.shape)

    # change with new calibration format
    def loadBackgroundCorrection(self):
        # backgroundCorrectionFilename = os.path.join(self.rawSiteFilesPath, 'backgroundCorrection.csv')
        # bg_corr_df = pd.read_csv(backgroundCorrectionFilename, delimiter=',', header=None)
        # self.backgroundCorrection = np.array(bg_corr_df)
        # print('Background Correction')
        # print(self.backgroundCorrection.shape)
        with h5py.File(self.cal_hdf, 'r') as f:
            self.backgroundCorrection = f['Background Correction Array'][:]
        # print(self.backgroundCorrection.shape)

    # We determined skimage does this equivilently, correct?
    # We should use that if so
    def equalizeHistogram(self):
        # Histogram Equalization to adjust contrast [1%-99%]
        # config file
        numberBins = 10000  # A good balance between time and space complexity, and well as precision
        imageHistogram, bins = np.histogram(self.imageArray1D, numberBins)
        imageHistogram = imageHistogram[1:]
        bins = bins[1:]
        cdf = np.cumsum(imageHistogram)

        # spliced to cut off non-image area
        # any way to determine this dynamically?  How periminant is it?
        cdf = cdf[:9996]

        max_cdf = max(cdf)
        contrast = self.contrast
        maxIndex = np.argmin(abs(cdf - contrast/100 * max_cdf))
        minIndex = np.argmin(abs(cdf - (100 - contrast)/100 * max_cdf))
        vmax = float(bins[maxIndex])
        vmin = float(bins[minIndex])
        lowValueIndices = self.imageArray1D < vmin
        self.imageArray1D[lowValueIndices] = vmin
        highValueIndices = self.imageArray1D > vmax
        self.imageArray1D[highValueIndices] = vmax
        self.imageData = self.imageArray1D.reshape(self.imageData.shape)

    # function not currently used - revise later?
    def removeStars(self):
        filteredData = copy.copy(self.imageData).astype(float)
        for i in range(self.width):
            leftPixels = self.imageData[:, max(i - 5, 0):max(i - 2, 1)]
            rightPixels = self.imageData[:, min(i + 3, self.width - 2):min(i + 6, self.width - 1)]
            allPixels = np.append(leftPixels, rightPixels, 1)
            pixelData = self.imageData[:, i]
            stdVal = np.std(allPixels, 1)
            meanVal = np.mean(allPixels, 1)
            for j in range(self.height):
                if pixelData[j] > meanVal[j] + 3 * stdVal[j]:
                    filteredData[j - 1:j + 2, i - 1:i + 2] = -1
                (iKnown, jKnown) = np.where(filteredData >= 0)
            valuesKnown = np.extract(filteredData >= 0, filteredData)
        # find minimum value from the imageData to find threshold
        m = np.empty(self.width)
        for i in range(self.width):
            m[i] = min(self.imageData[:, i])
        threshold = min(m)
        (iAll, jAll) = np.where(filteredData >= threshold)  # min value from imageData
        self.imageData = np.reshape(
            griddata((iKnown, jKnown), valuesKnown, (iAll, jAll), method='linear', fill_value=0), filteredData.shape)

    # # these will be read from the calibration file, not recalculated here, correct?
    # # If we don't currently have these values calcuatated from the calibration procedure, use the linear approximation
    # def setLensFunction(self):
    #     # Calculates lens function coefficients for the equation: Angle = a0 + a1.px + a2.px^2 + a3.px^3
    #     xDiff = self.zenithI - self.i
    #     yDiff = self.zenithJ - self.j
    #     distanceFromZenith = np.sqrt(xDiff ** 2 + yDiff ** 2)
    #     angleFromZenith = self.elevation[0] - self.elevation
    #     firstColumn = np.ones(len(distanceFromZenith))
    #
    #     angleFromZenithMat = np.vstack(
    #         [firstColumn, angleFromZenith, angleFromZenith ** 2, angleFromZenith ** 3]).transpose()
    #     self.angleToPixCoefficients = np.flip(np.dot(np.linalg.pinv(angleFromZenithMat), distanceFromZenith))
    #
    #     self.fisheyeRadius = self.getPixelsFromAngle(90)
    #
    # def getPixelsFromAngle(self, angle):
    #     # Input Angle in degrees
    #     return np.polyval(self.angleToPixCoefficients, angle)

    def calibrate(self):
        # # Spatial Calibration based on star data
        # # self.zenith = np.array([self.i[0], self.j[0]])
        # G_el = 1.0 - (self.getPixelsFromAngle(self.elevation) / self.fisheyeRadius)
        # self.f = G_el * np.sin(np.radians(self.azimuth))
        # self.g = G_el * np.cos(np.radians(self.azimuth))
        # firstColumn = np.ones(len(self.i))
        # oneIJ = np.vstack((firstColumn, self.i, self.j)).transpose()
        # oneIJInverse = slg.pinv(oneIJ)
        # aCoefficients = np.dot(oneIJInverse, self.f)
        # bCoefficients = np.dot(oneIJInverse, self.g)
        #
        # a0 = aCoefficients[0]
        # a1 = aCoefficients[1]
        # a2 = aCoefficients[2]
        # b0 = bCoefficients[0]
        # b1 = bCoefficients[1]
        # b2 = bCoefficients[2]
        #
        # rotationAngle_1 = np.degrees(np.arctan(-b1 / a1))
        # rotationAngle_2 = np.degrees(np.arctan(a2 / b2))
        # self.rotationAngle = .5 * (rotationAngle_1 + rotationAngle_2)
        # # just read rotation angle from calibration file, skip everything above this
        # # Let's correct this just in the calibration file for now
        # if 'EIO' in self.rawSiteFilesPath:
        #     self.rotationAngle = self.rotationAngle + 180

        self.rotationAngle = 10.
        self.imageData = np.fliplr(skimage.transform.rotate(self.imageData,
                                                            self.rotationAngle, order=3)).astype(float)

        # What exactly is zenith and why is this nessisary?

        # # Rotating Zenith Counter-clockwise by rotation angle and flipping it left-right
        # zenithI = self.width - int(np.cos(np.radians(self.rotationAngle)) * (self.zenith[0] - self.width / 2) - np.sin(
        #     np.radians(self.rotationAngle)) * (self.zenith[1] - self.height / 2) + self.width / 2)
        # zenithJ = int(np.sin(np.radians(self.rotationAngle)) * (self.zenith[0] - self.width / 2) + np.cos(
        #     np.radians(self.rotationAngle)) * (self.zenith[1] - self.height / 2) + self.height / 2)
        # self.zenith = [zenithI, zenithJ]
        # #self.setLensFunction()

        # why do we call the lens function after calibration?

    # What does this do??
    # is the effetively where the Fish-Eye dewarping occurs?
    # Should be able to add a NaN mask to just image array
    def mercatorUnwrap(self, ID_array):
        newImageWidth = 500
        newImageHeight = 500
        finalImage = np.ones([newImageHeight, newImageWidth]) * -1

        #it doesn't work for differently sized image arrays

        for j in range(ID_array.shape[0]):
            for i in range(ID_array.shape[1]):
                newI = self.newIMatrix[j][i]
                newJ = self.newJMatrix[j][i]
                if not (np.isnan(newI)) and not (np.isnan(newJ)) and not \
                        (np.isnan(self.backgroundCorrection[j, i])):
                    ID_array[j, i] = ID_array[j, i] * self.backgroundCorrection[j, i]
                    finalImage[int(newJ), int(newI)] = ID_array[j, i]

        (iKnown, jKnown) = np.where(finalImage >= 0)
        valuesKnown = np.extract(finalImage >= 0, finalImage)
        (iAll, jAll) = np.where(finalImage >= -1)
        interpolatedData = np.reshape(
            griddata((iKnown, jKnown), valuesKnown, (iAll, jAll), method='cubic', fill_value=-1),
            [newImageWidth, newImageHeight])
        alphaMask = np.ones([newImageHeight, newImageWidth]) * 255
        for j in range(interpolatedData.shape[0]):
            for i in range(interpolatedData.shape[1]):
                if interpolatedData[j, i] == -1:
                    # alphaMask[j, i] = 0
                    alphaMask[j, i] = np.nan  # to create transparent background

        interpolatedData = (interpolatedData * 255 / (np.nanmax(interpolatedData))).astype('uint8')
        alphaMask = alphaMask.astype('uint8')
        ID_array = np.dstack([interpolatedData, alphaMask])
        self.writeMode = 'LA'
        self.imageArray = ID_array

    def writePNG(self):
        self.imageArray = [np.array(self.imageArray)]
        if len(self.allImagesArray) == 0:
            self.allImagesArray = self.imageArray
        else:
            # This stacking should probably happen in main processing image script
            # MangoImage.py should only be aware of a single image file
            self.allImagesArray = np.append(self.allImagesArray, self.imageArray, axis=0)
