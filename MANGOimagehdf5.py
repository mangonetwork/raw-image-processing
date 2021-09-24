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
import copy
import skimage.transform
import h5py
from PIL import Image


class MANGOimage:
    # inheret Process Image class?
    def __init__(self, rawimg, config, collective_img_arr):

        self.rawImage = rawimg
        self.config = config
        self.cal_hdf = self.config['Data Locations']['cal_hdf']
        self.contrast = int(config['Specifications']['contrast'])

        self.load_image()
        self.load_caldata()
        self.process()

    def load_image(self):
        self.imageData = self.rawImage['image'][:]
        self.width = self.rawImage['image'].attrs['width']
        self.height = self.rawImage['image'].attrs['height']

    def load_caldata(self):
        # load data from calibration array
        with h5py.File(self.cal_hdf, 'r') as f:
            self.newIMatrix = f['New I array'][:]
            self.newJMatrix = f['New J array'][:]
            self.backgroundCorrection = f['Background Correction Array'][:]
            self.rotationAngle = f['Calibration Angle'][()]

    def process(self):
        self.equalizeHistogram()
        # self.showImage()
        # self.setLensFunction()
        self.calibrate()
        # self.showImage()
        self.mercatorUnwrap(self.imageData)
        # self.showImage()
        # self.writePNG()
        # self.showImage()

    # function for debugging?
    # fine to keep this, but make a note of it
    def showImage(self):
        img = Image.fromarray(self.imageData)
        img.show()


    def equalizeHistogram(self):
        # Histogram Equalization to adjust contrast [1%-99%]
        imageArray1D = self.imageData.flatten()
        # config file
        numberBins = 10000  # A good balance between time and space complexity, and well as precision
        imageHistogram, bins = np.histogram(imageArray1D, numberBins)
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
        lowValueIndices = imageArray1D < vmin
        imageArray1D[lowValueIndices] = vmin
        highValueIndices = imageArray1D > vmax
        imageArray1D[highValueIndices] = vmax
        self.imageData = imageArray1D.reshape(self.imageData.shape)

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


    def calibrate(self):
        self.imageData = np.fliplr(skimage.transform.rotate(self.imageData,
                                                            self.rotationAngle, order=3)).astype(float)


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
