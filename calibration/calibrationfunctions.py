##############################################################################
#
#       Written by:  Maria Rangel
#       Started:     02-06-2013       
#       Finished:    03-05-2013
#       Description: Contains all the functions used for the calibration pages
#
##############################################################################
from math import atan, asin, atan2
import numpy as np
from flask import *
from pylab import *
import matplotlib.pyplot as plt
import scipy.linalg as slg

mod = Blueprint('users', __name__, url_prefix='')

# SYSTEM VARIABLES
SITES_DIR = 'app/static/Sites/'
CAL_DIR = '/calibration/'
SYS_DIR = '/system/'
SCHEDULE_DIR = '/schedule/'
IMG_DIR = '/images/'
ALL_IMG_DIR = '/static/Sites/All Sites Images/'


# Data verification for site calibration
def lensfunction(zenith, azimuth, elevation, iCoord, jCoord, siteName):
    imageWidth = 695  # Atik 314L+ CCD at 2x2 binning
    imageHeight = 519  # Atik 314L+ CCD at 2x2 binning

    lenOfData = len(elevation)

    # Calculates the distance and angle from zenith
    distance = np.sqrt((zenith[1] - iCoord) ** 2 + (zenith[2] - jCoord) ** 2)
    angle = zenith[0] - elevation

    # Performs Fitting
    firstColumn = np.ones(len(angle))
    angleMat = np.vstack((firstColumn, angle, angle ** 2, angle ** 3)).transpose()
    angleCoefficients = np.dot(np.linalg.pinv(angleMat), distance)
    angleCoefficients = np.flip(angleCoefficients)

    firstColumn = np.ones(len(distance))
    distanceMat = np.vstack((firstColumn, distance, distance ** 2, distance ** 3)).transpose()
    distanceCoefficients = np.dot(np.linalg.pinv(distanceMat), angle)
    distanceCoefficients = np.flip(distanceCoefficients)
    distances = np.arange(260)
    testangles = np.polyval(distanceCoefficients, distances)

    anglefit = np.arange(91)
    pixfit = np.polyval(angleCoefficients, anglefit)

    empirical = np.polyval(angleCoefficients, angle)
    # observed = distance
    varianceDistance = np.std(distance) ** 2

    errorValue = (distance - empirical) ** 2
    residualSumSquares = np.sum(errorValue)

    plt.plot(anglefit, pixfit, 'r-')
    plt.scatter(angle, distance)
    plt.xlim((0, 91))
    plt.ylim(0)
    plt.grid(True)
    plt.ylabel('Distance')
    plt.xlabel('Angle')
    plt.title('Lens Function')
    plt.savefig(SITES_DIR + siteName + CAL_DIR + 'lensfunction.png')
    figureAddress = ('/static/Sites/' + siteName + CAL_DIR + 'lensfunction.png')
    plt.clf()

    # MODIFYING ZENITH
    G_el = 1.0 - (np.polyval(angleCoefficients, elevation) / np.polyval(angleCoefficients, 90))
    f = G_el * np.sin(np.radians(azimuth))
    g = G_el * np.cos(np.radians(azimuth))
    firstColumn = np.ones(len(iCoord))
    oneIJ = np.vstack((firstColumn, iCoord, jCoord)).transpose()
    oneIJInverse = slg.pinv(oneIJ)
    aCoefficients = np.dot(oneIJInverse, f)
    bCoefficients = np.dot(oneIJInverse, g)

    a0 = aCoefficients[0]
    a1 = aCoefficients[1]
    a2 = aCoefficients[2]
    b0 = bCoefficients[0]
    b1 = bCoefficients[1]
    b2 = bCoefficients[2]

    rotationAngle_1 = np.degrees(atan(-b1 / a1))
    rotationAngle_2 = np.degrees(atan(a2 / b2))
    rotationAngle = .5 * (rotationAngle_1 + rotationAngle_2)
    zenithI = int(np.cos(np.radians(rotationAngle)) * (zenith[1] - imageWidth / 2) - np.sin(np.radians(rotationAngle)) *
                  (zenith[2] - imageHeight / 2) + imageWidth / 2)
    zenithJ = int(np.sin(np.radians(rotationAngle)) * (zenith[1] - imageWidth / 2) + np.cos(np.radians(rotationAngle)) *
                  (zenith[2] - imageHeight / 2) + imageHeight / 2)
    deltaI = (zenithI - imageWidth / 2)
    zenithI = int(zenithI - 2 * deltaI)
    zenith = [zenith[0], zenithI, zenithJ]

    # DONE MODIFYING ZENITH

    return figureAddress, angle, distance, individualErrors, residualSumSquares, distanceCoefficients, zenith


# This function needs the zenith information (zenith elevation, zenith i coordinate, zenith j coordinate),
# pixels to Angle Coefficients which were calculated in the lens function, the zenith longitude (site longitude) ,
# zenith latitude (site latitude), and the site name. The function calculates the latitude and longitudes for each
# pixel and stores them in a CSV file. The function also calculates the southwest and northeast coordinates
def pixelsToLatLon(zenith, pixToAngleCoefficients, zenithLatitude, zenithLongitude, siteName):
    viewingAngle = 75
    imageWidth = 695  # Atik 314L+ CCD at 2x2 binning
    imageHeight = 519  # Atik 314L+ CCD at 2x2 binning
    newImageWidth = 500  # Web display images will be 500 pixels wide
    newImageHeight = 500  # Web display images will be 500 pixels tall
    zenithLatitude = float(zenithLatitude)  # converting the zenith latitude to a float
    zenithLongitude = float(zenithLongitude)  # converting the zenith longitude to a float
    zenithLatitude = np.radians(zenithLatitude)  # converting the zenith latitude from degrees to radians
    zenithLongitude = np.radians(zenithLongitude)  # converting the zenith latitude from degrees to radians

    latitudes = np.empty((imageHeight, imageWidth))  # creating an empty array with the same dimensions as the raw image
    longitudes = np.empty(
        (imageHeight, imageWidth))  # creating an empty array with the same dimensions as the raw image

    latitudes[:] = np.nan  # Setting the latitude array to NANs
    longitudes[:] = np.nan  # Setting the longitude array to NANs

    yDistance = imageHeight - zenith[2]
    xDistance = imageWidth - zenith[1]
    distanceFromZenith = np.sqrt(yDistance ** 2 + xDistance ** 2)
    angleFromZenith = np.polyval(pixToAngleCoefficients, distanceFromZenith)
    for j in range(imageHeight):
        for i in range(imageWidth):
            if angleFromZenith <= viewingAngle:
                earthDistance = earthDistanceFromZenithAngle(angleFromZenith)
                bearing = np.pi / 2 + atan2(yDistance, xDistance)
                [latitudes[j, i], longitudes[j, i]] = greatCircleArc(bearing, earthDistance, zenithLatitude,
                                                                     zenithLongitude)

    latBounds = np.array([np.nanmin(latitudes), np.nanmax(latitudes)])
    latRange = latBounds[1] - latBounds[0]
    lonBounds = np.array([np.nanmin(longitudes), np.nanmax(longitudes)])
    lonRange = lonBounds[1] - lonBounds[0]
    scaleBounds = np.log(np.tan(np.radians(latBounds) * 0.5 + np.pi / 4))
    scaleRange = scaleBounds[1] - scaleBounds[0]

    newImatrix = np.empty(
        (imageHeight, imageWidth))  # creating an empty array with the same dimensions as the raw image
    newJmatrix = np.empty(
        (imageHeight, imageWidth))  # creating an empty array with the same dimensions as the raw image
    newImatrix[:] = np.nan  # Setting the newI array to NANs
    newJmatrix[:] = np.nan  # Setting the newJ array to NANs

    for j in range(imageHeight):
        for i in range(imageWidth):
            pixLat = latitudes[j][i]
            pixLon = longitudes[j][i]
            if not np.isnan(pixLon) and not np.isnan(pixLat):
                newI = int((newImageWidth - 1) * (pixLon - lonBounds[0]) / lonRange)
                scaleJ = np.log(np.tan(np.radians(pixLat / 2.0) + np.pi / 4))
                newJ = (newImageHeight - 1) - int((newImageHeight - 1) * (scaleJ - scaleBounds[0]) / scaleRange)
                newImatrix[j, i] = newI
                newJmatrix[j, i] = newJ

    backgroundCorrection = np.empty(
        (imageHeight, imageWidth))  # creating an empty array with the same dimensions as the raw image
    backgroundCorrection[:] = np.nan  # Setting the backgroundCorrection array to NANs

    xDistance = imageWidth - zenith[1]
    yDistance = imageHeight - zenith[2]
    for j in range(imageHeight):
        for i in range(imageWidth):
            yDistance = j - zenith[2]
            xDistance = i - zenith[1]
            distanceFromZenith = sqrt(yDistance ** 2 + xDistance ** 2)
            angleFromZenith = polyval(pixToAngleCoefficients, distanceFromZenith)
            backgroundCorrection[j, i] = getBackgroundCorrection(angleFromZenith)

    # Stores the latitude and longitude path to file location
    newIPath = SITES_DIR + siteName + CAL_DIR + 'newI.csv'
    newJPath = SITES_DIR + siteName + CAL_DIR + 'newJ.csv'
    backgroundCorrectionPath = SITES_DIR + siteName + CAL_DIR + 'backgroundCorrection.csv'

    imageLatLonPath = SITES_DIR + siteName + CAL_DIR + 'imageLatLong.csv'

    # calculates the southwest's latitude and longitude and northeast's latitude and longtiude for image overlay for
    # the map
    southWestLatLon = [np.nanmin(latitudes), np.nanmin(longitudes)]
    northEastLatLon = [np.nanmax(latitudes), np.nanmax(longitudes)]
    imageLatitudesLongitudes = np.vstack(
        (southWestLatLon, northEastLatLon))  # creates an array stacked for writing to a CSV

    # Creates CSV files for Latitudes and Longitudes data for each pixel
    savetxt(backgroundCorrectionPath, backgroundCorrection, delimiter=",")  # writes the backgroungCorrection CSV file
    savetxt(newJPath, newJmatrix, delimiter=",")  # writes the newJ CSV file
    savetxt(newIPath, newImatrix, delimiter=",")  # writes the newI CSV file
    savetxt(imageLatLonPath, imageLatitudesLongitudes,
            delimiter=",")  # writes the image coordinates (sotuhwest and northeast) CSV file


def getBackgroundCorrection(angleFromZenith):
    # Effects corrected: Van Rhijn, Lens Vignetting, M. KUBOTA et al, 2001
    earthRadius = 6371.0  # kilometers
    airglowHeight = 400.0  # kilometers
    vanRhijnFactor = np.sqrt(
        1.0 - ((earthRadius / (earthRadius + airglowHeight)) ** 2 * np.sin(np.radians(angleFromZenith)) ** 2))
    # Effects corrected: Atmospheric Extinction
    a = 0.2  # atmospheric extinction coefficient, M. KUBOTA et al, 2001
    F = 1 / (np.cos(np.radians(angleFromZenith)) + 0.15 * (93.885 - angleFromZenith) ** (-1.253))
    extinctionFactor = 10.0 ** (0.4 * a * F)
    return vanRhijnFactor * extinctionFactor


def greatCircleArc(bearing, earthDistance, zenithLatitude, zenithLongitude):
    earthRadius = 6371.0
    try:
        lat = np.asin(np.sin(zenithLatitude) * np.cos(earthDistance / earthRadius) +
                      np.cos(zenithLatitude) * np.sin(earthDistance / earthRadius) * np.cos(bearing))
        lon = zenithLongitude + atan2(np.sin(bearing) * np.sin(earthDistance / earthRadius) * np.cos(zenithLatitude),
                                      np.cos(earthDistance / earthRadius) - np.sin(zenithLatitude) * np.sin(lat))
        lat = np.degrees(lat)
        lon = np.degrees(lon)
        return lat, lon
    except:
        flash('Failure in function greatCircleArc(...)!')
        # raise('Failure in function greatCircleArc(...)!')


def earthDistanceFromZenithAngle(zenithAngle):
    try:
        earthRadius = 6371.0  # kilometers
        airglowHeight = 400.0  # kilometers
        a = 1
        b = -2 * earthRadius * cos(radians(180 - zenithAngle))
        c = -(airglowHeight ** 2) - (2 * earthRadius * airglowHeight)
        possibleDistances = roots([a, b, c])
        distance = max(possibleDistances)
        earthAngle = arcsin(distance * sin(radians(180 - zenithAngle)) / (earthRadius + airglowHeight))
        earthDistance = earthRadius * earthAngle
        return earthDistance
    except:
        flash('Failure in function earthDistanceFromZenithAngle(...)!')
        # raise('Failure in function earthDistanceFromZenithAngle(...)!')
