import warnings
import configparser

from .MANGOImage import MANGOImage
import argparse
import re
import pandas as pd
# import pymap3d as pm
import datetime as dt
import numpy as np
import h5py

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")

# Need docstrings throughout

# PASSING ARGUMENTS vs class attributes???

class ProcessImage:

    def __init__(self, configFile, inputList, outputFile):
        self.inputList = inputList
        self.outputFile = outputFile

        self.RE = 6371.0
        self.calParams = {}

        self.read_config(configFile)
        self.read_calfile()
        self.process_images()
        self.write_to_hdf5()

    def read_config(self, configFile):
        # read in config file
        config = configparser.ConfigParser()
        config.read(configFile)

        self.siteName = config.get('PROCESSING','SITE_NAME')
        self.calFile = config.get('CALIBRATION','CALIBRATION_FILE')
        self.calParams['contrast'] = config.getint('PROCESSING','CONTRAST')
        # specify in config file
        self.newImax = config.getint('PROCESSING','NEWIMAX')
        self.newJmax = config.getint('PROCESSING','NEWJMAX')
        # self.elevCutoff = config.getfloat('PROCESSING','ELEVCUTOFF')
        # read from calibration file?
        # self.targAlt = config.getfloat('PROCESSING','IMAGEALTITUDE')


    def read_calfile(self):

        with h5py.File(self.calFile, 'r') as f:
            self.calParams['transformedCoords'] = f['transformed_coords'][:]
            self.calParams['atmosphericCorrection'] = f['atmos_corr'][:]
            self.calParams['mask'] = f['mask'][:]
            self.ha = f['altitude'][()]
            self.elevCutoff = f['elevation_cutoff'][()]
            # self.calParams['newIMatrix'] = f['New I array'][:]
            # self.calParams['newJMatrix'] = f['New J array'][:]
            # self.calParams['backgroundCorrection'] = f['Background Correction Array'][:]
            # self.calParams['rotationAngle'] = f['Calibration Angle'][()]
            # self.latitude = f['Latitude'][:]
            # self.longitude = f['Longitude'][:]
        # self.longitude[self.longitude<0.] +=360.


    def process_images(self):
        self.startTime = np.array([])
        self.exposureTime = np.array([])
        self.ccdTemp = np.array([])
        self.rawList = self.inputList
        self.imageArrays = []
        self.get_time_independent_information(self.rawList[0])
        self.create_position_arrays()
        for file in self.rawList:
            print(file)
            with h5py.File(file, 'r') as hdf5_file:
                imageData = hdf5_file['image'][:]
                self.get_time_dependent_information(hdf5_file['image'])
                # picture = MANGOimagehdf5.MANGOimage(hdf5_file['image'], self.calParams)

            image = MANGOImage(imageData)
            image.equalizeHistogram(self.calParams['contrast'])
            # self.showImage()
            # self.setLensFunction()
            # image.rotateImage(self.calParams['rotationAngle'])
            # self.showImage()
            # image.mercatorUnwrap(self.calParams['newIMatrix'], self.calParams['newJMatrix'], self.calParams['backgroundCorrection'])
            image.transformImage(self.calParams['transformedCoords'], self.calParams['atmosphericCorrection'], self.calParams['mask'], self.pixelCoords)
            image.applyMask(self.imageMask)
            # self.showImage()
            # self.writePNG()
            # self.showImage()

            self.imageArrays.append(image.imageData)

        self.imageArrays = np.array(self.imageArrays)


    def get_time_independent_information(self, file):
        with h5py.File(file, 'r') as hdf5_file:
            data = hdf5_file['image']
            self.code = data.attrs['station']
            self.siteLat = data.attrs['latitude']
            self.siteLon = data.attrs['longitude']
            self.bin_x = data.attrs['bin_x']
            self.bin_y = data.attrs['bin_y']
            self.label = data.attrs['label']
            # get this from starttime rather than file name?
            data_split = re.findall(r'-(\d+)-', file)[0]
            self.siteDate = dt.datetime.strptime(data_split, '%Y%m%d').date()

    # consider appending this in loop?
    # Just a style preference
    def get_time_dependent_information(self, img):
        '''
        Obtains the following attributes:
        1. start_time
        2. exposure_time
        3. ccd_temp
        :param file: input hdf5 file
        :return: None
        '''

        self.startTime = np.append(self.startTime, img.attrs['start_time'])
        self.exposureTime = np.append(self.exposureTime, img.attrs['start_time'])
        self.ccdTemp = np.append(self.ccdTemp, img.attrs['start_time'])


    # def create_position_arrays(self):
    #
    #     i_grid, j_grid = np.meshgrid(np.arange(self.newImax),np.arange(self.newJmax))
    #     RL = self.newJmax/2.
    #     x_grid = (self.newImax/2.-i_grid)/RL
    #     y_grid = (self.newJmax/2.-j_grid)/RL
    #     r2_grid = x_grid**2 + y_grid**2
    #     z_grid = np.sqrt(1. - r2_grid)
    #
    #     # Set mask based on elevation angle limit
    #     r_cutoff = np.cos(self.elevCutoff*np.pi/180.)
    #     self.imageMask = r2_grid > r_cutoff**2
    #     x_grid[self.imageMask] = np.nan
    #     y_grid[self.imageMask] = np.nan
    #     z_grid[self.imageMask] = np.nan
    #
    #     self.azimuth = np.arctan2(x_grid, y_grid)*180./np.pi
    #     self.elevation = np.arccos(np.sqrt(x_grid**2 + y_grid**2))*180./np.pi
    #
    #     x, y, z = pm.geodetic2ecef(self.site_lat, self.site_lon, 0.)
    #     vx, vy, vz = pm.enu2uvw(x_grid, y_grid, z_grid, self.site_lat, self.site_lon)
    #
    #     earth = pm.Ellipsoid()
    #     a2 = (earth.semimajor_axis + self.targAlt*1000.)**2
    #     b2 = (earth.semimajor_axis + self.targAlt*1000.)**2
    #     c2 = (earth.semiminor_axis + self.targAlt*1000.)**2
    #
    #     A = vx**2/a2 + vy**2/b2 + vz**2/c2
    #     B = x*vx/a2 + y*vy/b2 + z*vz/c2
    #     C = x**2/a2 + y**2/b2 + z**2/c2 -1
    #
    #     alpha = (np.sqrt(B**2-A*C)-B)/A
    #
    #     self.latitude, self.longitude, alt = pm.ecef2geodetic(x + alpha*vx, y + alpha*vy, z + alpha*vz)

    def create_position_arrays(self):

        transformedX = self.calParams['transformedCoords'][0][~self.calParams['mask']]
        transformedY = self.calParams['transformedCoords'][1][~self.calParams['mask']]
        newX = np.linspace(np.min(transformedX),np.max(transformedX), self.newImax)
        newY = np.linspace(np.min(transformedY),np.max(transformedY), self.newJmax)

        newXGrid, newYGrid = np.meshgrid(newX, newY)
        self.pixelCoords = np.array([newXGrid,newYGrid])

        psi = np.sqrt(newXGrid**2 + newYGrid**2)/(self.RE+self.ha)
        c = np.sqrt(self.RE**2 + (self.RE+self.ha)**2 - 2*self.RE*(self.RE+self.ha)*np.cos(psi))  # Law of Cosines
        b = (self.RE+self.ha) - self.RE*np.cos(psi)
        elv = np.pi/2. - np.arccos(b/c) - psi
        azm = np.arctan2(newXGrid, newYGrid)
        self.imageMask = elv<self.elevCutoff*np.pi/180.

        # Use Haversine equations to find lat/lon of each point
        lat = np.arcsin(np.sin(self.siteLat*np.pi/180.) * np.cos(psi) +
                      np.cos(self.siteLat*np.pi/180.) * np.sin(psi) * np.cos(azm))
        lon = self.siteLon*np.pi/180. + np.arctan2(np.sin(azm) * np.sin(psi) * np.cos(self.siteLat*np.pi/180.),
                                      np.cos(psi) - np.sin(self.siteLat*np.pi/180.) * np.sin(lat))

        # import pymap3d as pm
        # x, y, z = pm.geodetic2ecef(self.siteLat, self.siteLon, 0.)
        # vx, vy, vz = pm.enu2uvw(np.cos(elv)*np.sin(azm), np.cos(elv)*np.cos(azm), np.sin(elv), self.siteLat, self.siteLon)
        #
        # earth = pm.Ellipsoid()
        # a2 = (earth.semimajor_axis + self.ha*1000.)**2
        # b2 = (earth.semimajor_axis + self.ha*1000.)**2
        # c2 = (earth.semiminor_axis + self.ha*1000.)**2
        #
        # A = vx**2/a2 + vy**2/b2 + vz**2/c2
        # B = x*vx/a2 + y*vy/b2 + z*vz/c2
        # C = x**2/a2 + y**2/b2 + z**2/c2 -1
        #
        # alpha = (np.sqrt(B**2-A*C)-B)/A
        #
        # self.latitude, self.longitude, alt = pm.ecef2geodetic(x + alpha*vx, y + alpha*vy, z + alpha*vz)

        self.azimuth = azm*180./np.pi
        self.elevation = elv*180./np.pi
        self.latitude = lat*180./np.pi
        self.longitude = lon*180./np.pi


    def write_to_hdf5(self):

        self.endTime = self.startTime + self.exposureTime
        # tstmp_s = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.startTime])
        # tstmp_e = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in self.endTime])
        tstmp_s = self.startTime
        tstmp_e = self.endTime

        # save hdf5 file
        with h5py.File(self.outputFile, 'w') as f:
            f.create_group('SiteInfo')

            T = f.create_dataset('UnixTime', data=np.array([tstmp_s, tstmp_e]), compression='gzip', compression_opts=1)
            T.attrs['Description'] = 'unix time stamp'
            T.attrs['Unit'] = 'seconds'

            # 250 km should be read in from somewhere - calibration file?
            # this is associated with the locations of the latitude, so that would make sense
            # Consider just copying array with attached attributes from calibration file?
            Lat = f.create_dataset('Latitude', data=self.latitude, compression='gzip', compression_opts=1)
            Lat.attrs['Description'] = 'geodetic latitude of each pixel'
            Lat.attrs['Size'] = 'Ipixels x Jpixels'
            Lat.attrs['Projection Altitude'] = self.ha
            Lat.attrs['Unit'] = 'degrees'

            Lon = f.create_dataset('Longitude', data=self.longitude, compression='gzip', compression_opts=1)
            Lon.attrs['Description'] = 'geodetic longitude of each pixel'
            Lon.attrs['Size'] = 'Ipixels x Jpixels'
            Lon.attrs['Projection Altitude'] = self.ha
            Lon.attrs['Unit'] = 'degrees'

            az = f.create_dataset('Azimuth', data=self.azimuth, compression='gzip', compression_opts=1)
            az.attrs['Description'] = 'azimuth of each pixel'
            az.attrs['Size'] = 'Ipixels x Jpixels'
            az.attrs['Unit'] = 'degrees'

            el = f.create_dataset('Elevation', data=self.elevation, compression='gzip', compression_opts=1)
            el.attrs['Description'] = 'elevation of each pixel'
            el.attrs['Size'] = 'Ipixels x Jpixels'
            el.attrs['Unit'] = 'degrees'

            el = f.create_dataset('PixelCoordinates', data=self.pixelCoords, compression='gzip', compression_opts=1)
            el.attrs['Description'] = 'coordinates of each pixel in grid at the airglow altitude'
            el.attrs['Size'] = '2 (X,Y) x Ipixels x Jpixels'
            el.attrs['Unit'] = 'km'

            I = f.create_dataset('ImageData', data=self.imageArrays, compression='gzip', compression_opts=1)
            I.attrs['Description'] = 'pixel values for images'
            I.attrs['Site Abbreviation'] = self.code
            I.attrs['Image label'] = self.label
            # I.attrs['x/y binning'] = np.array([self.bin_x, self.bin_y])

            M = f.create_dataset('Mask', data=self.imageMask, compression='gzip', compression_opts=1)
            M.attrs['Description'] = 'image mask'

            CCD = f.create_dataset('CCDTemperature', data=self.ccdTemp)
            CCD.attrs['Description'] = 'Temperature of CCD'
            CCD.attrs['Unit'] = 'degrees C'

            N = f.create_dataset('SiteInfo/Name', data=self.siteName)
            N.attrs['Description'] = 'site name'
            C = f.create_dataset('SiteInfo/Code', data=self.code)
            C.attrs['Description'] = 'one letter site abbreviation/code'
            L = f.create_dataset('SiteInfo/Coordinates', data=[self.siteLat, self.siteLon])
            L.attrs['Description'] = 'geodetic coordinates of site; [lat, lon]'
            L.attrs['Unit'] = 'degrees'


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
    parser.add_argument('-i', '--input', dest='inputs', nargs='+', type=str,
                        help='A list of .hdf5 files to process and store in output file.')
    parser.add_argument('-o', '--output', dest='output', type=str,
                        help='Output file to write processed images to.')
    args = parser.parse_args()

    return args


def main():
    command_line_args = parse_args()
    conf = command_line_args.config
    inputs = command_line_args.inputs
    output = command_line_args.output
    ProcessImage(conf, inputs, output)

if __name__ == '__main__':
    main()
