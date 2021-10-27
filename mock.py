import h5py
import os
import numpy as np
import pandas as pd

def write_to_mock():

    f = h5py.File('mock.hdf5', 'w')
    sitepath = "C:\\Users\\padma\\MANGO SU21\\raw_data\\CFS\\site_files\\calibration"
    newIFilename = os.path.join(sitepath, 'newI.csv')
    nifarray = np.array(pd.read_csv(newIFilename))

    newJFilename = os.path.join(sitepath, 'newJ.csv')
    njfarray = np.array(pd.read_csv(newJFilename))

    backgroundCorrectionFilename = os.path.join(sitepath, 'backgroundCorrection.csv')
    bcfarray = np.array(pd.read_csv(backgroundCorrectionFilename))

    calibrationFile = "C:\\Users\\padma\\MANGO SU21\\raw_data\\CFS\\site_files\\calibration\\Calibration.csv"
    calarray = np.array(pd.read_csv(calibrationFile))

    latitudeFile = "C:\\Users\\padma\\MANGO SU21\\raw_data\\CFS\\site_files\\calibration\\Latitudes.csv"
    latarray = np.array(pd.read_csv(latitudeFile))

    longitudeFile = "C:\\Users\\padma\\MANGO SU21\\raw_data\\CFS\\site_files\\calibration\\Longitudes.csv"
    lonarray = np.array(pd.read_csv(longitudeFile))

    Lat = f.create_dataset('Latitude', data=latarray, compression='gzip', compression_opts=1)
    Lat.attrs['Description'] = 'geodetic latitude of each pixel projected to 250 km'
    Lat.attrs['Size'] = 'Ipixels x Jpixels'
    Lat.attrs['Projection Altitude'] = 250
    Lat.attrs['Unit'] = 'degrees'


    Lon = f.create_dataset('Longitude', data=lonarray, compression='gzip', compression_opts=1)
    Lon.attrs['Description'] = 'geodetic longitude of each pixel projected to 250 km'
    Lon.attrs['Size'] = 'Ipixels x Jpixels'
    Lon.attrs['Projection Altitude'] = 250
    Lon.attrs['Unit'] = 'degrees'

    newI = f.create_dataset('New I array', data=nifarray, compression='gzip', compression_opts=1)
    newI.attrs['Description'] = 'New I array used for masking in mercator unwrapping function in MANGOimage.py'

    newJ = f.create_dataset('New J array', data=njfarray, compression='gzip', compression_opts=1)
    newJ.attrs['Description'] = 'New J array used for masking in mercator unwrapping function in MANGOimage.py'

    backgroundCorrection = f.create_dataset('Background Correction Array', data=bcfarray, compression='gzip', compression_opts=1)
    backgroundCorrection.attrs['Description'] = 'Background correction array used for masking in mercator unwrapping function in MANGOimage.py'

    calibration = f.create_dataset('Calibration Angle', data=[200], compression='gzip', compression_opts=1)
    calibration.attrs['Unit'] = 'degrees (to rotate anticlockwise)'
    f.close()


if __name__ == '__main__':
    write_to_mock()