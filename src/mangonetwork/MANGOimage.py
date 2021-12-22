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

# script containing image manipulation functions for processing and quicklooks


from scipy.interpolate import griddata

import numpy as np
import copy
import skimage.transform
import h5py
from PIL import Image



def equalizeHistogram(imageData, contrast):
    # This function contains some questionable hard-coded things?
    # Histogram Equalization to adjust contrast [1%-99%]
    imageArray1D = imageData.flatten()
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
    # contrast = self.calParams['contrast']
    maxIndex = np.argmin(abs(cdf - contrast/100 * max_cdf))
    minIndex = np.argmin(abs(cdf - (100 - contrast)/100 * max_cdf))
    vmax = float(bins[maxIndex])
    vmin = float(bins[minIndex])
    lowValueIndices = imageArray1D < vmin
    imageArray1D[lowValueIndices] = vmin
    highValueIndices = imageArray1D > vmax
    imageArray1D[highValueIndices] = vmax
    imageData = imageArray1D.reshape(imageData.shape)

    return imageData


# function not currently used - revise later?
def removeStars(imageData, width, height):
    filteredData = copy.copy(imageData).astype(float)
    for i in range(width):
        leftPixels = imageData[:, max(i - 5, 0):max(i - 2, 1)]
        rightPixels = imageData[:, min(i + 3, width - 2):min(i + 6, width - 1)]
        allPixels = np.append(leftPixels, rightPixels, 1)
        pixelData = imageData[:, i]
        stdVal = np.std(allPixels, 1)
        meanVal = np.mean(allPixels, 1)
        for j in range(height):
            if pixelData[j] > meanVal[j] + 3 * stdVal[j]:
                filteredData[j - 1:j + 2, i - 1:i + 2] = -1
            (iKnown, jKnown) = np.where(filteredData >= 0)
        valuesKnown = np.extract(filteredData >= 0, filteredData)
    # find minimum value from the imageData to find threshold
    m = np.empty(width)
    for i in range(width):
        m[i] = min(imageData[:, i])
    threshold = min(m)
    (iAll, jAll) = np.where(filteredData >= threshold)  # min value from imageData
    imageData = np.reshape(
        griddata((iKnown, jKnown), valuesKnown, (iAll, jAll), method='linear', fill_value=0), filteredData.shape)

    return imageData


def rotateImage(imageData, angle):
    imageData = np.fliplr(skimage.transform.rotate(imageData,
                                                        angle, order=3)).astype(float)
    return imageData


# What does this do??
# Effetively where the Fish-Eye dewarping occurs
# Certain points are mapped to new coordinates and then the rest of the array is filled in with interpolation
# There are probably better ways to do this ...
# Should be able to add a NaN mask to just image array
# imageData = data
# imageArray = data + masking array
def mercatorUnwrap(imageData, newIMatrix, newJMatrix, backgroundCorrection):
    ID_array = imageData
    newImageWidth = 500
    newImageHeight = 500
    finalImage = np.ones([newImageHeight, newImageWidth]) * -1

    #it doesn't work for differently sized image arrays

    for j in range(ID_array.shape[0]):
        for i in range(ID_array.shape[1]):
            newI = newIMatrix[j][i]
            newJ = newJMatrix[j][i]
            if not (np.isnan(newI)) and not (np.isnan(newJ)) and not \
                    (np.isnan(backgroundCorrection[j, i])):
                ID_array[j, i] = ID_array[j, i] * backgroundCorrection[j, i]
                finalImage[int(newJ), int(newI)] = ID_array[j, i]


    # most of this can probably be done with np.where(np.isfinite())
    (iKnown, jKnown) = np.where(finalImage >= 0)
    valuesKnown = np.extract(finalImage >= 0, finalImage)
    (iAll, jAll) = np.where(finalImage >= -1)
    interpolatedData = np.reshape(
        griddata((iKnown, jKnown), valuesKnown, (iAll, jAll), method='cubic', fill_value=-1),
        [newImageWidth, newImageHeight])

    # This only needs to be done once and is slow
    # Move to outer ProcessImage class
    alphaMask = np.ones([newImageHeight, newImageWidth]) * 255
    alphaMask[interpolatedData == -1] = 0
    alphaMask = alphaMask.astype('uint8')
    # for j in range(interpolatedData.shape[0]):
    #     for i in range(interpolatedData.shape[1]):
    #         if interpolatedData[j, i] == -1:
    #             # alphaMask[j, i] = 0
    #             alphaMask[j, i] = np.nan  # to create transparent background


    # interpolatedData = (interpolatedData * 255 / (np.nanmax(interpolatedData))).astype('uint8')
    # alphaMask = alphaMask.astype('uint8')
    # ID_array = np.dstack([interpolatedData, alphaMask])
    # self.writeMode = 'LA'
    # self.imageArray = ID_array
    interpolatedData = (interpolatedData * 255 / (np.nanmax(interpolatedData)))
    # interpolatedData[np.isnan(alphaMask)] = np.nan
    # alphaMask = alphaMask.astype('uint8')
    # ID_array = np.dstack([interpolatedData, alphaMask])
    # self.writeMode = 'LA'
    # self.imageArray = interpolatedData
    imageData = interpolatedData.astype('uint8')

    return imageData, alphaMask


# function for debugging?
# fine to keep this, but make a note of it
def showImage(imageData):
    img = Image.fromarray(imageData.astype('uint8'))
    img.show()
