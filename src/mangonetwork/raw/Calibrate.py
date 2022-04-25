# Calibrate.py
import numpy as np
import pandas as pd
import h5py
import os
import argparse
import configparser


class Calibrate:

    def __init__(self, configFile):

        self.calParams = {}

        self.read_config(configFile)

        self.star_data = np.loadtxt(self.starCalFile, usecols=(1,2,3,4))

        self.RE = 6371.0

        self.find_calibration_params()
        self.create_calibration_arrays()
        self.save_calibration_file()


    def read_config(self, configFile):
        # read in config file
        config = configparser.ConfigParser()
        config.read(configFile)

        # self.siteCalFilepath = config['MOCK']['SITE_CAL_FILEPATH']
        # site calibration folders are all located in:
        # ~/Sites/<site name>/calibration
        # self.rotationAngle = float(config['MOCK']['ROTATION_ANGLE'])
        self.outputFile = config.get('CALIBRATION','CALIBRATION_FILE')
        self.starCalFile = config.get('CALIBRATION', 'STAR_CAL_FILE')
        self.ha = config.getfloat('CALIBRATION', 'ALTITUDE')
        self.elevCutoff = config.getfloat('CALIBRATION', 'ELEV_CUTOFF')
        self.Imax = config.getint('CALIBRATION','IMAX')
        self.Jmax = config.getint('CALIBRATION','JMAX')
        # self.newImax = config.getint('CALIBRATION','NEWIMAX')
        # self.newJmax = config.getint('CALIBRATION','NEWJMAX')
        # self.siteLat = config.getfloat('CALIBRATION','SITELAT')
        # self.siteLon = config.getfloat('CALIBRATION','SITELON')

    def normalize_pixel_coords(self, i, j):
        RL = self.Jmax/2.
        x = (self.Imax/2.-i)/RL
        y = (self.Jmax/2.-j)/RL
        return x, y

    def find_calibration_params(self):

        star_el = self.star_data[:,0]*np.pi/180.
        star_az = self.star_data[:,1]*np.pi/180.
        i = self.star_data[:,2]
        j = self.star_data[:,3]
        x, y = self.normalize_pixel_coords(i, j)
        xp = np.cos(star_el)*np.sin(star_az)
        yp = np.cos(star_el)*np.cos(star_az)


        from scipy.optimize import least_squares

        init_params = self.initial_params(x, y, xp, yp)
        params = least_squares(self.residuals, init_params, args=(x, y, xp, yp))
        self.x0, self.y0, self.theta, self.A, self.B, self.C, self.D = params.x


    def transform(self, x, y, x0, y0, theta, A, B, C, D, unwarp=True, return_elev=False):
        x1 = x - x0
        y1 = y - y0

        t = theta*np.pi/180.
        x2 = np.cos(t)*x1 - np.sin(t)*y1
        y2 = np.sin(t)*x1 + np.cos(t)*y1

        r = np.sqrt(x2**2+y2**2)
        lam = A + B*r + C*r**2 + D*r**3

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

    def residuals(self, params, x, y, xp, yp):
        x0, y0, theta, A, B, C, D = params
        xt, yt = self.transform(x, y, x0, y0, theta, A, B, C, D, unwarp=False)
        res = np.sqrt((xp-xt)**2 + (yp-yt)**2)
        return res

    def initial_params(self, x, y, xp, yp):
        from scipy.spatial.transform import Rotation

        # Get these from config file?
        # TBH, unless a camera is oriented REALLY wierdly, these initial parameters should be fine
        x0, y0 = [0., 0.]
        A, B, C, D = [np.pi/2, -np.pi/2, 0., 0.]

        # calculate transformation with initial tranlation and lens function params but no rotation
        xu, yu = self.transform(x, y, x0, y0, 0., A, B, C, D, unwarp=False)

        # Find rotation matrix such that the vectors to the star locations roughly match
        Pu = np.array([xu, yu, np.zeros(len(x))]).T
        Pp = np.array([xp, yp, np.zeros(len(x))]).T
        R, _ = Rotation.align_vectors(Pp, Pu)
        # Find euler angles of rotation matrix and select "z" rotation as an approximate theta
        theta = R.as_euler('xyz', degrees=True)[2]

        return [x0, y0, theta, A, B, C, D]


    def create_calibration_arrays(self):

        # Create a grid and find the transformed coordinates of each grid point
        i_grid, j_grid = np.meshgrid(np.arange(self.Imax),np.arange(self.Jmax))
        x_grid, y_grid = self.normalize_pixel_coords(i_grid, j_grid)

        self.xt_grid, self.yt_grid, lam_grid = self.transform(x_grid, y_grid, self.x0, self.y0, self.theta, self.A, self.B, self.C, self.D, return_elev=True)

        # create mask
        self.mask = np.zeros((self.Jmax, self.Imax), dtype=bool)
        self.mask[lam_grid<self.elevCutoff*np.pi/180.] = True

        # atmospheric corrections
        self.atmosphericCorrection = self.atmospheric_corrections(lam_grid)
        self.atmosphericCorrection[self.mask] = np.nan


    def atmospheric_corrections(self, elev):
        # Atmospheric corrections are taken from Kubota et al., 2001
        # Kubota, M., Fukunishi, H. & Okano, S. Characteristics of medium- and
        #   large-scale TIDs over Japan derived from OI 630-nm nightglow observation.
        #   Earth Planet Sp 53, 741–751 (2001). https://doi.org/10.1186/BF03352402

        # calculate zenith angle
        za = np.pi/2-elev

        # Kubota et al., 2001; eqn. 6
        vanRhijnFactor = np.sqrt(1.0 - np.sin(za)**2*(self.RE/(self.RE+self.ha))**2)

        # Kubota et al., 2001; eqn. 7,8
        a = 0.2
        F = 1. / (np.cos(za) + 0.15 * (93.885 - za*180./np.pi) ** (-1.253))
        extinctionFactor = 10.0 ** (0.4 * a * F)

        return vanRhijnFactor * extinctionFactor



    # def create_new_position_arrays():
    #
    #     newX = np.linspace(np.min(self.xt_grid[~self.mask]),np.max(self.xt_grid[~self.mask]), self.newImax)
    #     newY = np.linspace(np.min(self.yt_grid[~self.mask]),np.max(self.yt_grid[~self.mask]), self.newJmax)
    #
    #     newXGrid, newYGrid = np.meshgrid(newX, newY)
    #
    #     psi = np.sqrt(newXGrid**2 + newYGrid**2)/(self.RE+self.ha)
    #     c = np.sqrt(self.RE**2 + self.ha**2 - 2*self.RE*(self.RE+self.ha)*np.cos(psi))  # Law of Cosines
    #     elv = np.arcsin((self.RE+self.ha)/c*np.sin(psi)) - np.pi/2.     # Law of Sines
    #     azm = np.arctan2(newXGrid, newYGrid)
    #
    #     lat = np.arcsin(np.sin(self.siteLat*np.pi/180.) * np.cos(psi) +
    #                   np.cos(self.siteLat*np.pi/180.) * np.sin(psi) * np.cos(azm))
    #     lon = self.siteLon*np.pi/180. + np.arctan2(np.sin(azm) * np.sin(psi) * np.cos(self.siteLon*np.pi/180.),
    #                                   np.cos(psi) - np.sin(self.siteLat*np.pi/180.) * np.sin(lat))
    #
    #     self.latitude = lat*180./np.pi
    #     self.longitude = lon*180./np.pi




    def save_calibration_file(self):

        with h5py.File(self.outputFile, 'w') as f:

            newI = f.create_dataset('transformed_coords', data=np.array([self.xt_grid, self.yt_grid]), compression='gzip', compression_opts=1)
            newI.attrs['Description'] = 'Transformed x, y coordinates of each point in the image grid'

            backgroundCorrection = f.create_dataset('atmos_corr', data=self.atmosphericCorrection, compression='gzip', compression_opts=1)
            backgroundCorrection.attrs['Description'] = 'Background atmospheric correction factors to apply to the image pre-regridding'

            mask = f.create_dataset('mask', data=self.mask, compression='gzip', compression_opts=1)
            backgroundCorrection.attrs['Description'] = 'Mask of portion of array that does not include the ASI FoV'

            # transformation/calibration parameters
            lensfun = f.create_dataset('lens_function', data=np.array([self.A, self.B, self.C, self.D]))
            lensfun.attrs['Description'] = 'Lens function coefficients of the form [A,B,C,D] (A+Br+Cr^2+Dr^3)'

            trans = f.create_dataset('translation', data=np.array([self.x0, self.y0]))
            trans.attrs['Description'] = 'Linear translation [x0, y0]'

            rotate = f.create_dataset('rotation_angle', data=self.theta)
            rotate.attrs['Description'] = 'Rotation angle'

            alt = f.create_dataset('altitude', data=self.ha)
            alt.attrs['Description'] = 'Assumed airglow altituded used for atmospheric correction (km)'

            elev = f.create_dataset('elevation_cutoff', data=self.elevCutoff)
            elev.attrs['Description'] = 'Lowest allowed elevation angle'
    # def mock(self):
    #     print('WARNING: THIS IS A MOCK ROUTINE.  Instead of calculating new calibration parameters, it converts the old calibration files to the new format and assigns a pre-designated rotation angle.')
    #
    #     f = h5py.File(self.outputFile, 'w')
    #     newIFilename = os.path.join(self.siteCalFilepath, 'newI.csv')
    #     nifarray = np.array(pd.read_csv(newIFilename, header=None))
    #
    #     newJFilename = os.path.join(self.siteCalFilepath, 'newJ.csv')
    #     njfarray = np.array(pd.read_csv(newJFilename, header=None))
    #
    #     backgroundCorrectionFilename = os.path.join(self.siteCalFilepath, 'backgroundCorrection.csv')
    #     bcfarray = np.array(pd.read_csv(backgroundCorrectionFilename, header=None))
    #
    #     calibrationFile = os.path.join(self.siteCalFilepath, 'Calibration.csv')
    #     calarray = np.array(pd.read_csv(calibrationFile, header=None))
    #
    #     latitudeFile = os.path.join(self.siteCalFilepath, 'Latitudes.csv')
    #     latarray = np.array(pd.read_csv(latitudeFile, header=None))
    #
    #     longitudeFile = os.path.join(self.siteCalFilepath, 'Longitudes.csv')
    #     lonarray = np.array(pd.read_csv(longitudeFile, header=None))
    #
    #     Lat = f.create_dataset('Latitude', data=latarray, compression='gzip', compression_opts=1)
    #     Lat.attrs['Description'] = 'geodetic latitude of each pixel projected to 250 km'
    #     Lat.attrs['Size'] = 'Ipixels x Jpixels'
    #     Lat.attrs['Projection Altitude'] = 250
    #     Lat.attrs['Unit'] = 'degrees'
    #
    #
    #     Lon = f.create_dataset('Longitude', data=lonarray, compression='gzip', compression_opts=1)
    #     Lon.attrs['Description'] = 'geodetic longitude of each pixel projected to 250 km'
    #     Lon.attrs['Size'] = 'Ipixels x Jpixels'
    #     Lon.attrs['Projection Altitude'] = 250
    #     Lon.attrs['Unit'] = 'degrees'
    #
    #     newI = f.create_dataset('New I array', data=nifarray, compression='gzip', compression_opts=1)
    #     newI.attrs['Description'] = 'New I array used for masking in mercator unwrapping function in MANGOimage.py'
    #
    #     newJ = f.create_dataset('New J array', data=njfarray, compression='gzip', compression_opts=1)
    #     newJ.attrs['Description'] = 'New J array used for masking in mercator unwrapping function in MANGOimage.py'
    #
    #     backgroundCorrection = f.create_dataset('Background Correction Array', data=bcfarray, compression='gzip', compression_opts=1)
    #     backgroundCorrection.attrs['Description'] = 'Background correction array used for masking in mercator unwrapping function in MANGOimage.py'
    #
    #     calibration = f.create_dataset('Calibration Angle', data=[self.rotationAngle], compression='gzip', compression_opts=1)
    #     calibration.attrs['Unit'] = 'degrees (to rotate anticlockwise)'
    #     f.close()


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
    args = parser.parse_args()

    return args


def main():
    command_line_args = parse_args()
    conf = command_line_args.config
    Calibrate(conf)


if __name__ == '__main__':
    main()
