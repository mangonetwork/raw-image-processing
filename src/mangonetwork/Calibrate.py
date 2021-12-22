# Calibrate.py
import numpy as np
import pandas as pd
import h5py
import os
import argparse
import configparser


class Calibrate:

    def __init__(self, configFile, outputFile):
        self.outputFile = outputFile

        self.calParams = {}

        self.read_config(configFile)
        self.mock()

    def read_config(self, configFile):
        # read in config file
        config = configparser.ConfigParser()
        config.read(configFile)

        self.siteCalFilepath = config['MOCK']['SITE_CAL_FILEPATH']
        # site calibration folders are all located in:
        # ~/Sites/<site name>/calibration
        self.rotationAngle = float(config['MOCK']['ROTATION_ANGLE'])

    def mock(self):
        print('WARNING: THIS IS A MOCK ROUTINE.  Instead of calculating new calibration parameters, it converts the old calibration files to the new format and assigns a pre-designated rotation angle.')

        f = h5py.File(self.outputFile, 'w')
        newIFilename = os.path.join(self.siteCalFilepath, 'newI.csv')
        nifarray = np.array(pd.read_csv(newIFilename, header=None))

        newJFilename = os.path.join(self.siteCalFilepath, 'newJ.csv')
        njfarray = np.array(pd.read_csv(newJFilename, header=None))

        backgroundCorrectionFilename = os.path.join(self.siteCalFilepath, 'backgroundCorrection.csv')
        bcfarray = np.array(pd.read_csv(backgroundCorrectionFilename, header=None))

        calibrationFile = os.path.join(self.siteCalFilepath, 'Calibration.csv')
        calarray = np.array(pd.read_csv(calibrationFile, header=None))

        latitudeFile = os.path.join(self.siteCalFilepath, 'Latitudes.csv')
        latarray = np.array(pd.read_csv(latitudeFile, header=None))

        longitudeFile = os.path.join(self.siteCalFilepath, 'Longitudes.csv')
        lonarray = np.array(pd.read_csv(longitudeFile, header=None))

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

        calibration = f.create_dataset('Calibration Angle', data=[self.rotationAngle], compression='gzip', compression_opts=1)
        calibration.attrs['Unit'] = 'degrees (to rotate anticlockwise)'
        f.close()


def parse_args():
    """
    Handle the command line arguments.
    Returns:
    Output of argparse.ArgumentParser.parse_args.
    """

    parser = argparse.ArgumentParser(description='Accepting config, input files and output file'
                                                 'to process MANGOImage.')
    parser.add_argument('-c', '--config', dest='config', type=str,
                        help='Config file containing data locations and image specs.')
    parser.add_argument('-o', '--output', dest='output', type=str,
                        help='Output file to write processed images to.')
    args = parser.parse_args()

    return args


def main():
    command_line_args = parse_args()
    conf = command_line_args.config
    output = command_line_args.output
    Calibrate(conf, output)


if __name__ == '__main__':
    main()
