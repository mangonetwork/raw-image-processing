##############################################################################
#
#       Written by:  Maria Rangel
#       Started:     02-06-2013       
#       Finished:    03-05-2013
#       Description: Contains all the functions used for the calibration pages
#
##############################################################################
from flask                  import *
from numpy                  import *
from StringIO               import StringIO
from math                   import atan, acos, asin, atan2
import matplotlib           as plt
import matplotlib.pyplot    as plt
from pylab                  import *

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
    imageWidth = 695                            # Atik 314L+ CCD at 2x2 binning
    imageHeight = 519                           # Atik 314L+ CCD at 2x2 binning

    lenOfData = len(elevation)

    # Creating empty arrays   
    distance = []
    angle = []
    individualErrors = []
    residualSumSquares=[]

    
    # Calculates the distance and angle from zenith
    for k in range(0,lenOfData ):
        temp = sqrt(((zenith[1]- iCoord[k])**2)+((zenith[2]- jCoord[k])**2))        # distance formula
        distance.append(temp)
        temp2 = (zenith[0] - elevation[k])     # angle from zenith = subtract the angle from the star from the zenith angle
        angle.append(temp2)
            
    # Performs Fitting
    firstColumn = array(([1]*len(angle)))
    angleMat = vstack((firstColumn, angle, [x**2 for x in angle], [x**3 for x in angle])).transpose()
    angleCoefficients = dot(linalg.pinv(angleMat), distance)
    angleCoefficients = angleCoefficients[::-1] 

    firstColumn = array(([1]*len(angle)))
    distanceMat = vstack((firstColumn, distance, [x**2 for x in distance], [x**3 for x in distance])).transpose()
    distanceCoefficients = dot(linalg.pinv(distanceMat), angle)
    distanceCoefficients = distanceCoefficients[::-1] 
    distances = range(260)
    testangles = polyval(distanceCoefficients, distances)

    anglefit = range(91)
    pixfit = polyval(angleCoefficients, anglefit)


    empirical = polyval(angleCoefficients, angle)
    #observed = distance
    varianceDistance = (std(distance)**2)

    residualSumSquares = 0
    individualErrors = []
    for i in range(len(distance)):
        errorValue = ((distance[i]-empirical[i])**2) #residual sum of squares
        individualErrors = append(individualErrors, errorValue)
        residualSumSquares += errorValue
 
    plot(anglefit, pixfit, 'r-')    
    scatter(angle, distance)
    xlim((0,91))
    ylim(0)
    grid(True)
    ylabel('Distance')
    xlabel('Angle')
    title('Lens Function')
    savefig(SITES_DIR + siteName + CAL_DIR + 'lensfunction.png')
    figureAddress = ('/static/Sites/' + siteName + CAL_DIR + 'lensfunction.png')
    clf()
    
    ###### MODIFYING ZENITH ###########
    G_el = 1.0 - [polyval(angleCoefficients, value) for value in elevation]/polyval(angleCoefficients, 90)
    f = G_el*sin(radians(azimuth))
    g = G_el*cos(radians(azimuth))
    firstColumn = array(([1]*len(iCoord)))
    oneIJ = vstack((firstColumn, iCoord, jCoord)).transpose()
    intermediate0 = dot(oneIJ.transpose(), oneIJ)
    intermediate1 = linalg.pinv(intermediate0)
    intermediate2 = dot(intermediate1, oneIJ.transpose())

    aCoefficients = dot(intermediate2, f)
    bCoefficients = dot(intermediate2, g)
    
    a0 = aCoefficients[0]
    a1 = aCoefficients[1]
    a2 = aCoefficients[2]
    b0 = bCoefficients[0]
    b1 = bCoefficients[1]
    b2 = bCoefficients[2]

    
    rotationAngle_1 = degrees(atan(-b1/a1))
    rotationAngle_2 = degrees(atan(a2/b2))
    rotationAngle = .5*(rotationAngle_1 + rotationAngle_2)
    zenithI = int(cos(radians(rotationAngle))*(zenith[1] - imageWidth/2) - sin(radians(rotationAngle))*(zenith[2]-imageHeight/2) + imageWidth/2)
    zenithJ = int(sin(radians(rotationAngle))*(zenith[1] - imageWidth/2) + cos(radians(rotationAngle))*(zenith[2]-imageHeight/2) + imageHeight/2)
    deltaI = (zenithI-imageWidth/2)
    zenithI = int(zenithI - 2*deltaI)
    zenith = [zenith[0], zenithI, zenithJ]
   
    ######### DONE MODIFYING ZENITH ########


    return figureAddress, angle, distance, individualErrors, residualSumSquares, distanceCoefficients, zenith
    

# This function needs the zenith information (zenith elevation, zenith i coordinate, zenith j coordinate), pixels to Angle Coefficients
# which were calculated in the lens function, the zenith longitude (site longitude) , zenith latitude (site latitude), and the site name. 
# The function calculates the latitude and logitudes for each pixel and stores them in a CSV file. The function also calculates the southwest
# and northeast coordinates 
def pixelsToLatLon(zenith, pixToAngleCoefficients, zenithLatitude, zenithLongitude, siteName):
    viewingAngle = 75
    imageWidth = 695                            # Atik 314L+ CCD at 2x2 binning
    imageHeight = 519                           # Atik 314L+ CCD at 2x2 binning
    newImageWidth = 500                         # Web display images will be 500 pixels wide
    newImageHeight = 500                        # Web display images will be 500 pixels tall
    zenithLatitude = float(zenithLatitude)      # converting the zenith latitude to a float
    zenithLongitude = float(zenithLongitude)    # converting the zenith longitude to a float
    zenithLatitude = radians(zenithLatitude)    # converting the zenith latitude from degrees to radians
    zenithLongitude = radians(zenithLongitude)  # converting the zenith latitude from degrees to radians
    
    latitudes = empty(([imageHeight, imageWidth]))          # creating an empty array with the same dimensions as the raw image
    longitudes = empty(([imageHeight, imageWidth]))         # creating an empty array with the same dimensions as the raw image
   
    latitudes[:] = NAN          # Setting the laitude array to NANs
    longitudes[:] = NAN         # Setting the longitude array to NANs

    for j in range(imageHeight):
        for i in range(imageWidth):           
            yDistance = j - zenith[2]
            xDistance = i - zenith[1]
            distanceFromZenith = sqrt(yDistance**2 + xDistance**2)
            angleFromZenith = polyval(pixToAngleCoefficients, distanceFromZenith)
            if angleFromZenith<=viewingAngle:
                earthDistance = earthDistanceFromZenithAngle(angleFromZenith)
                bearing = pi/2 + atan2(yDistance, xDistance)
                [latitudes[j, i], longitudes[j, i]] = greatCircleArc(bearing, earthDistance, zenithLatitude, zenithLongitude)
    


    latBounds = [nanmin(latitudes), nanmax(latitudes)]
    latRange = latBounds[1]-latBounds[0]
    lonBounds = [nanmin(longitudes), nanmax(longitudes)]
    lonRange = lonBounds[1]-lonBounds[0]
    scaleBounds = log(tan(radians(latBounds)*0.5 + pi/4))
    scaleRange = scaleBounds[1]-scaleBounds[0]

    newImatrix = empty(([imageHeight, imageWidth]))          # creating an empty array with the same dimensions as the raw image
    newJmatrix = empty(([imageHeight, imageWidth]))         # creating an empty array with the same dimensions as the raw image   
    newImatrix[:] = NAN          # Setting the newI array to NANs
    newJmatrix[:] = NAN         # Setting the newJ array to NANs

    for j in range(imageHeight):
            for i in range(imageWidth):
                pixLat = latitudes[j][i]
                pixLon = longitudes[j][i]
                if (not(isnan(pixLon)) and not(isnan(pixLat))):
                    newI = int((newImageWidth-1)*(pixLon - lonBounds[0])/lonRange)
                    scaleJ = log(tan(radians(pixLat/2.0) + pi/4))
                    newJ = (newImageHeight-1) - int((newImageHeight-1)*(scaleJ- scaleBounds[0])/scaleRange)
                    newImatrix[j, i] = newI
                    newJmatrix[j, i] = newJ
    

    backgroundCorrection    = empty(([imageHeight, imageWidth]))    # creating an empty array with the same dimensions as the raw image
    backgroundCorrection[:] = NAN                                   # Setting the backgroundCorrection array to NANs
    
    for j in range(imageHeight):
            for i in range(imageWidth):
                    yDistance = j - zenith[2]
                    xDistance = i - zenith[1]
                    distanceFromZenith = sqrt(yDistance**2 + xDistance**2)
                    angleFromZenith = polyval(pixToAngleCoefficients, distanceFromZenith)
                    backgroundCorrection[j, i] = getBackgroundCorrection(angleFromZenith)

    # Stores the latitude and longitude path to file location
    newIPath = SITES_DIR + siteName + CAL_DIR + 'newI.csv'
    newJPath = SITES_DIR + siteName + CAL_DIR + 'newJ.csv'
    backgroundCorrectionPath = SITES_DIR + siteName + CAL_DIR + 'backgroundCorrection.csv'

    imageLatLonPath = SITES_DIR + siteName + CAL_DIR + 'imageLatLong.csv'

    # calculates the southwest's latitude and longitude and northeast's latitude and longtiude for image overlay for the map
    southWestLatLon = [nanmin(latitudes), nanmin(longitudes)]
    northEastLatLon = [nanmax(latitudes), nanmax(longitudes)]
    imageLatitudesLongitudes = array((vstack((southWestLatLon, northEastLatLon))))  # creates an array stacked for writing to a CSV
    
    # Creates CSV files for Latitudes and Longitudes data for each pixel
    savetxt(backgroundCorrectionPath , backgroundCorrection, delimiter=",")     # writes the backgroungCorrection CSV file
    savetxt(newJPath , newJmatrix, delimiter=",")                               # writes the newJ CSV file
    savetxt(newIPath , newImatrix, delimiter=",")                               # writes the newI CSV file
    savetxt(imageLatLonPath, imageLatitudesLongitudes, delimiter=",")           # writes the image coordinates (sotuhwest and northeast) CSV file

def getBackgroundCorrection(angleFromZenith):
    #Effects corrected: Van Rhijn, Lens Vignetting, M. KUBOTA et al, 2001
    earthRadius = 6371.0      #kilometers
    airglowHeight = 400.0     #kilometers
    vanRhijnFactor = sqrt(1.0 - ((earthRadius/(earthRadius + airglowHeight))**2 * sin(radians(angleFromZenith))**2))
    #Effects corrected: Atmospheric Extinction
    a = 0.2                 #atmospheric extinction coefficient, M. KUBOTA et al, 2001
    F = 1/(cos(radians(angleFromZenith)) + 0.15*(93.885 - angleFromZenith)**(-1.253))
    extinctionFactor = 10.0**(0.4*a*F)
    return vanRhijnFactor*extinctionFactor

def greatCircleArc(bearing, earthDistance, zenithLatitude, zenithLongitude):
    earthRadius = 6371.0
    try:    
        lat = asin(sin(zenithLatitude)*cos(earthDistance/earthRadius) + 
                  cos(zenithLatitude)*sin(earthDistance/earthRadius)*cos(bearing))
        lon = zenithLongitude + atan2(sin(bearing)*sin(earthDistance/earthRadius)*cos(zenithLatitude), 
                         cos(earthDistance/earthRadius)-sin(zenithLatitude)*sin(lat))
        lat = degrees(lat)
        lon = degrees(lon)
        return lat, lon
    except: 
        flash('Failure in function greatCircleArc(...)!')

def earthDistanceFromZenithAngle(zenithAngle):
    try:
        earthRadius = 6371.0      #kilometers
        airglowHeight = 400.0     #kilometers
        a = 1
        b = -2*earthRadius*cos(radians(180 - zenithAngle))
        c = -(airglowHeight**2)-(2*earthRadius*airglowHeight)
        possibleDistances = roots([a, b, c])
        distance = max(possibleDistances)
        earthAngle = arcsin(distance*sin(radians(180 - zenithAngle))/(earthRadius+airglowHeight))
        earthDistance = earthRadius*earthAngle
        return earthDistance
    except:
        flash('Failure in function earthDistanceFromZenithAngle(...)!')
