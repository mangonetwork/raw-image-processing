######################################################################
#
#   Written by:  Maria Rangel 02-06-2013
#   Description: Code performs coordinate mapping verification
#
######################################################################
from numpy    import *
from math     import *




def kilometersPerPixels():
    filename= 'CoordinateMappingVerification.csv'
    data = genfromtxt(filename, delimiter = ',') # reads csv file

    # Known variables
    airglowheight= 250.00                         # km
    
    # Extracting data and converting to array
    starname =  array((data[:,0]))               # extracts the star name from data
    elevation = array((data[:,1]))               # extracts the elevation from data
    icoord =    array((data[:,2]))               # extracts the i coord from data
    jcoord =    array((data[:,3]))               # extracts the star name from data
    zenith =    array((data[10,:]))              # extract zenith from data

    lenofdata = len(data)                        # gives the length of csv file

    # Creating empty arrays   
    distance = []
    theta =  []
    kilometers = []
    kmperpix = []
    
    # Calculates the kilometers per pixels
    for k in range(11,lenofdata):
        # Distance
        dis = sqrt(((zenith[2]- icoord[k])**2)+((zenith[3]- jcoord[k])**2))
        distance.append(dis)
        # Theta
        zenithangle = (90 - elevation[k])
        zenithangle = float(zenithangle)
        zenithangle = radians(zenithangle)
        theta.append(zenithangle)
        print theta
        # Kilometers
        km= 0.0
        km = airglowheight*tan(zenithangle)
        kilometers.append(km)
        # Kilometers per pixel
        pixDis = 0.0
        pixDis= km/dis
        kmperpix.append(pixDis)

    return kmperpix
    

