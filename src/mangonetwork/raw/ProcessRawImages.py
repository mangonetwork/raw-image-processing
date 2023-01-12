import warnings
import configparser
import logging

from .MANGOImage import MANGOImage
import argparse
import re
import pandas as pd
# import pymap3d as pm
import datetime as dt
import numpy as np
import h5py
import sys
import os
from scipy.interpolate import griddata

if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources

warnings.filterwarnings("ignore", message="Reloaded modules: MANGOimage")

# Need docstrings throughout

# PASSING ARGUMENTS vs class attributes???

class ProcessImage:

    def __init__(self, configFile, inputList, outputFile):
        self.inputList = inputList
        self.outputFile = outputFile

        self.RE = 6371.0
        self.calParams = {}
        # Pull these from file?
        self.Imax = 695
        self.Jmax = 519

        self.read_config(configFile)
        # self.read_calfile()

        self.get_time_independent_information(self.inputList[0])
        self.create_transform_grids()
        self.create_position_arrays()

        # import matplotlib.pyplot as plt
        # plt.scatter(self.transXGrid[::10,::10], self.transYGrid[::10,::10], c=self.elevGrid[::10,::10]*180./np.pi, vmin=14., vmax=16.)
        # plt.scatter(self.newXGrid[::50,::50], self.newYGrid[::50,::50], c=self.elevation[::50,::50], vmin=14., vmax=16., edgecolors='magenta')
        # # plt.scatter(x_grid[::50,::50], y_grid[::50,::50])
        # # plt.scatter(self.x0, self.y0)
        # plt.scatter([-self.rl, self.rl, 0., 0.], [0., 0., -self.rl, self.rl])
        # plt.show()



        self.atmosCorr = self.atmospheric_corrections()

        self.process_images()
        self.write_to_hdf5()

    def read_config(self, config):
        # # read in config file
        # config = configparser.ConfigParser()
        # config.read(configFile)

        self.siteName = config.get('PROCESSING','SITE_NAME')
        # self.calFile = config.get('CALIBRATION','CALIBRATION_FILE')
        self.calParams['contrast'] = config.getint('PROCESSING','CONTRAST')
        # specify in config file
        self.newImax = config.getint('PROCESSING','NEWIMAX')
        self.newJmax = config.getint('PROCESSING','NEWJMAX')
        self.elevCutoff = config.getfloat('PROCESSING','ELEVCUTOFF')
        # read from calibration file?
        # self.targAlt = config.getfloat('PROCESSING','IMAGEALTITUDE')
        # self.Imax = config.getint('CALIBRATION', 'IMAX')
        # self.Jmax = config.getint('CALIBRATION', 'JMAX')
        self.ha = config.getfloat('PROCESSING', 'ALTITUDE')

        self.x0 = config.getfloat('CALIBRATION_PARAMS', 'X0')
        self.y0 = config.getfloat('CALIBRATION_PARAMS', 'Y0')
        self.rl = config.getfloat('CALIBRATION_PARAMS', 'RL')
        self.theta = config.getfloat('CALIBRATION_PARAMS', 'THETA')
        self.A = config.getfloat('CALIBRATION_PARAMS', 'A')
        self.B = config.getfloat('CALIBRATION_PARAMS', 'B')
        self.C = config.getfloat('CALIBRATION_PARAMS', 'C')
        self.D = config.getfloat('CALIBRATION_PARAMS', 'D')



    # def read_calfile(self):
    #
    #     with h5py.File(self.calFile, 'r') as f:
    #         self.calParams['transformedCoords'] = f['transformed_coords'][:]
    #         self.calParams['atmosphericCorrection'] = f['atmos_corr'][:]
    #         self.calParams['mask'] = f['mask'][:]
    #         self.ha = f['altitude'][()]
    #         self.elevCutoff = f['elevation_cutoff'][()]
    #         # self.calParams['newIMatrix'] = f['New I array'][:]
    #         # self.calParams['newJMatrix'] = f['New J array'][:]
    #         # self.calParams['backgroundCorrection'] = f['Background Correction Array'][:]
    #         # self.calParams['rotationAngle'] = f['Calibration Angle'][()]
    #         # self.latitude = f['Latitude'][:]
    #         # self.longitude = f['Longitude'][:]
    #     # self.longitude[self.longitude<0.] +=360.

    # def normalize_pixel_coords(self, i, j):
    #     RL = self.Jmax/2.
    #     x = (self.Imax/2.-i)/RL
    #     y = (self.Jmax/2.-j)/RL
    #     return x, y

    # move these functions to processing class
    def create_transform_grids(self):

        # Create a grid and find the transformed coordinates of each grid point
        x_grid, y_grid = np.meshgrid(np.arange(self.Imax),np.arange(self.Jmax))
        # x_grid, y_grid = self.normalize_pixel_coords(i_grid, j_grid)    # Don't need this - happens internally


        self.transXGrid, self.transYGrid, self.elevGrid = self.transform(x_grid, y_grid, self.x0, self.y0, self.rl, self.theta, self.A, self.B, self.C, self.D, return_elev=True)

        # transformedX = self.transXGrid[elevGrid>self.elevCutoff*np.pi/180.]
        # transformedY = self.transYGrid[elevGrid>self.elevCutoff*np.pi/180.]
        # newX = np.linspace(np.min(transformedX),np.max(transformedX), self.newImax)
        # newY = np.linspace(np.min(transformedY),np.max(transformedY), self.newJmax)

        # find boundary based on where the lens function equals the cutoff angle
        # rb = np.roots([self.D, self.C, self.B, self.A-self.elevCutoff*np.pi/180.])
        # el_fac = rb[np.isreal(rb)][0].real
        # print(el_fac*self.rl)
        b = np.arcsin(self.RE/(self.RE+self.ha)*np.sin(self.elevCutoff*np.pi/180.+np.pi/2))  # Law of Sines
        psi = np.pi/2-self.elevCutoff*np.pi/180.-b
        d = (self.RE+self.ha)*psi


        newX = np.linspace(-d, d, self.newImax)
        newY = np.linspace(-d, d, self.newJmax)

        self.newXGrid, self.newYGrid = np.meshgrid(newX, newY)
        # self.pixelCoords = np.array([newXGrid,newYGrid])


        # # create mask - possibly redundant?
        # self.mask = np.zeros((self.Jmax, self.Imax), dtype=bool)
        # self.mask[lam_grid<self.elevCutoff*np.pi/180.] = True
        #
        # # atmospheric corrections
        # self.atmosphericCorrection = self.atmospheric_corrections(lam_grid)
        # self.atmosphericCorrection[self.mask] = np.nan

    def transform(self, x, y, x0, y0, rl, theta, A, B, C, D, unwarp=True, return_elev=False):
        x1 = (x - x0)/rl
        y1 = (y - y0)/rl

        t = theta*np.pi/180.
        x2 = np.cos(t)*x1 - np.sin(t)*y1
        y2 = np.sin(t)*x1 + np.cos(t)*y1
        # x2 = np.cos(t)*x1 + np.sin(t)*y1
        # y2 = -np.sin(t)*x1 + np.cos(t)*y1

        r = np.sqrt(x2**2+y2**2)
        lam = A + B*r + C*r**2 + D*r**3

        # import matplotlib.pyplot as plt
        # plt.scatter(r, lam*180./np.pi)
        # plt.grid()
        # plt.show()

        if unwarp:
            # unwarping
            b = np.arcsin(self.RE/(self.RE+self.ha)*np.sin(lam+np.pi/2))  # Law of Sines
            psi = np.pi/2-lam-b
            d = (self.RE+self.ha)*psi
        else:
            d = np.cos(lam)

        x3 = d*x2/r
        y3 = d*y2/r

        if return_elev:
            return x3, y3, lam
        else:
            return x3, y3


    def create_position_arrays(self):

        psi = np.sqrt(self.newXGrid**2 + self.newYGrid**2)/(self.RE+self.ha)
        c = np.sqrt(self.RE**2 + (self.RE+self.ha)**2 - 2*self.RE*(self.RE+self.ha)*np.cos(psi))  # Law of Cosines
        b = (self.RE+self.ha) - self.RE*np.cos(psi)
        elv = np.pi/2. - np.arccos(b/c) - psi
        azm = np.arctan2(self.newXGrid, self.newYGrid)
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

        # import matplotlib.pyplot as plt
        # plt.scatter(np.sqrt(self.newXGrid**2 + self.newYGrid**2), elv*180./np.pi)
        # plt.grid()
        # plt.show()



        self.azimuth = azm*180./np.pi
        self.elevation = elv*180./np.pi
        self.latitude = lat*180./np.pi
        self.longitude = lon*180./np.pi


    def atmospheric_corrections(self):
        # Atmospheric corrections are taken from Kubota et al., 2001
        # Kubota, M., Fukunishi, H. & Okano, S. Characteristics of medium- and
        #   large-scale TIDs over Japan derived from OI 630-nm nightglow observation.
        #   Earth Planet Sp 53, 741–751 (2001). https://doi.org/10.1186/BF03352402

        # calculate zenith angle
        za = np.pi/2-self.elevation*np.pi/180.

        # Kubota et al., 2001; eqn. 6
        vanRhijnFactor = np.sqrt(1.0 - np.sin(za)**2*(self.RE/(self.RE+self.ha))**2)

        # Kubota et al., 2001; eqn. 7,8
        a = 0.2
        F = 1. / (np.cos(za) + 0.15 * (93.885 - za*180./np.pi) ** (-1.253))
        extinctionFactor = 10.0 ** (0.4 * a * F)

        return vanRhijnFactor * extinctionFactor



    def process_images(self):
        self.startTime = np.array([])
        self.exposureTime = np.array([])
        self.ccdTemp = np.array([])
        self.rawList = self.inputList
        self.imageArrays = []
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
            # image.transformImage(np.array([self.xt_grid, self.yt_grid]), self.atmosphericCorrection, self.mask, self.pixelCoords)

            # Flip image first
            # image.invertImage()

            # img_corr = self.imageData
            # xt_grid = transformedCoords[0]
            # yt_grid = transformedCoords[1]
            #
            # x_grid = newCoords[0]
            # y_grid = newCoords[1]

            # import matplotlib.pyplot as plt
            # # plt.pcolormesh(image.imageData)
            # # plt.scatter(self.x0,self.y0)
            # plt.scatter(self.transXGrid, self.transYGrid, c=image.imageData, s=2)
            # # plt.scatter(self.newXGrid[::50,::50], self.newYGrid[::50,::50])
            # # plt.scatter(x_grid[::50,::50], y_grid[::50,::50])
            # # plt.scatter(self.x0, self.y0)
            # plt.show()


            # interpolatedData = griddata((xt_grid[~mask], yt_grid[~mask]), img_corr[~mask], (x_grid, y_grid), fill_value=0)
            newImage = griddata((self.transXGrid.flatten(), self.transYGrid.flatten()), image.imageData.flatten(), (self.newXGrid, self.newYGrid), fill_value=0)

            # The rest of these functions can possibly be moved outside the loop
            # Atmospheric correction
            newImage = newImage*self.atmosCorr

            # Apply mask outside elevation cutoff
            newImage[self.elevation<self.elevCutoff] = 0.

            # Renormalize each image and convert to int
            newImage = (newImage * 255 / np.nanmax(newImage)).astype('uint8')
            # self.imageData = interpolatedData.astype('uint8')

            # image.transformImage(np.array([self.transXGrid, self.yt_grid]), self.pixelCoords)
            # # Apply atmospheric corrections
            # image.applyMask(self.imageMask)
            # # self.showImage()
            # # self.writePNG()
            # # self.showImage()

            self.imageArrays.append(newImage)

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
        self.exposureTime = np.append(self.exposureTime, img.attrs['exposure_time'])
        self.ccdTemp = np.append(self.ccdTemp, img.attrs['ccd_temp'])


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

            el = f.create_dataset('PixelCoordinates', data=np.array([self.newXGrid,self.newYGrid]), compression='gzip', compression_opts=1)
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


# def parse_args():
#     """
#     Handle the command line arguments.
#     Returns:
#     Output of argparse.ArgumentParser.parse_args.
#     """
#
#     parser = argparse.ArgumentParser(description='Accepting config, input files and output file'
#                                                  'to process MANGOImage.')
#     parser.add_argument('-c', '--config', dest='config', type=str,
#                         help='Config file containing data locations and image specs.')
#     parser.add_argument('-i', '--input', dest='inputs', nargs='+', type=str,
#                         help='A list of .hdf5 files to process and store in output file.')
#     parser.add_argument('-o', '--output', dest='output', type=str,
#                         help='Output file to write processed images to.')
#     args = parser.parse_args()
#
#     return args
#
#
# def main():
#     command_line_args = parse_args()
#     conf = command_line_args.config
#     inputs = command_line_args.inputs
#     output = command_line_args.output
#     ProcessImage(conf, inputs, output)
#
# if __name__ == '__main__':
#     main()
#

def parse_args():

    parser = argparse.ArgumentParser(description='Create a quick look movie')

    parser.add_argument('-c', '--config', metavar='FILE',
                        help='Alternate configuration file')
    parser.add_argument('-f', '--filelist', metavar='FILE',
                        help='A file with a list of .hdf5 file names')
    parser.add_argument('-o', '--output', default='mango.hdf5',
                        help='Output filename (default is mango.hdf5)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')

    parser.add_argument('inputfiles', nargs='*')

    return parser.parse_args()

def find_config(filename):

    with h5py.File(filename) as h5:
        station = h5['image'].attrs['station']
        instrument = h5['image'].attrs['instrument']

    name = '%s-%s.ini' % (station, instrument)
    print(name)

    logging.debug('Using package configuration file: %s' % name)

    return resources.files('mangonetwork.raw.data').joinpath(name).read_text()

def main():

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.filelist:
        inputs = [line.strip() for line in open(args.filelist)]
        inputs = [line for line in inputs if line]
    else:
        inputs = args.inputfiles

    if not inputs:
        logging.error('No input files listed')
        sys.exit(1)

    logging.debug('Processing %d files' % len(inputs))

    if args.config:
        logging.debug('Alternate configuration file: %s' % args.config)
        if not os.path.exists(args.config):
            logging.error('Config file not found')
            sys.exit(1)
        contents = open(args.config).read()
    else:
        contents = find_config(inputs[0])

    config = configparser.ConfigParser()
    config.read_string(contents)

    logging.debug('Configuration file: %s' % args.config)

    ProcessImage(config, inputs, args.output)

    sys.exit(0)

if __name__ == '__main__':
    main()
