#!/usr/bin/env python
"""Calibration"""

##########################################################################
#
#   Calibration
#
#   2022-xx-xx  Leslie Lamarche and Asti Bhatt
#               Initial implementation
#
#   2023-03-08  Todd Valentic
#               Make PEP8 compliant
#
##########################################################################

import argparse
import configparser
import io
import logging
import os
import sys

import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation

if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources


class Calibrate:
    """Camera calibration"""

    def __init__(self, config_file, starcal_file, output):
        self.find_calibration_params(starcal_file)
        self.save_calibration_params(output, config_file)

    def find_calibration_params(self, starcal_file):
        """Load calibration parameters from starcal file"""

        # read in data from starcal file
        star_az, star_el, x, y = np.loadtxt(
            io.StringIO(starcal_file), usecols=(1, 2, 3, 4), unpack=True
        )

        # true x,y positions of stars
        xp = np.cos(star_el * np.pi / 180.0) * np.sin(star_az * np.pi / 180.0)
        yp = np.cos(star_el * np.pi / 180.0) * np.cos(star_az * np.pi / 180.0)

        init_params = self.initial_params(x, y, xp, yp)
        params = least_squares(self.residuals, init_params, args=(x, y, xp, yp))
        self.x0, self.y0, self.rl, self.theta, self.C, self.D = params.x

        # NOTE: A and B are fully constrained when fitting for rl
        self.A = np.pi / 2.0
        self.B = -(np.pi / 2.0 + self.C + self.D)

        ## DEBUG: To confirm star locations match after transformation
        #xt, yt = self.transform(x, y, self.x0, self.y0, self.rl,
        #                          self.theta, self.A, self.B, self.C, self.D)
        #import matplotlib.pyplot as plt
        #plt.scatter(xp, yp)
        #plt.scatter(xt, yt)
        #plt.show()

    # pylint: disable=too-many-arguments, too-many-locals

    def transform(self, x, y, x0, y0, rl, theta, A, B, C, D):
        """Transformation"""

        x1 = (x - x0) / rl
        y1 = (y - y0) / rl

        t = theta * np.pi / 180.0
        x2 = np.cos(t) * x1 - np.sin(t) * y1
        y2 = np.sin(t) * x1 + np.cos(t) * y1

        r = np.sqrt(x2**2 + y2**2)
        lam = A + B * r + C * r**2 + D * r**3
        d = np.cos(lam)

        x3 = d * x2 / r
        y3 = d * y2 / r

        return x3, y3

    def residuals(self, params, x, y, xp, yp):
        """Residuals"""

        x0, y0, rl, theta, C, D = params
        A = np.pi / 2.0
        B = -(np.pi / 2.0 + C + D)
        xt, yt = self.transform(x, y, x0, y0, rl, theta, A, B, C, D)
        res = np.sqrt((xp - xt) ** 2 + (yp - yt) ** 2)
        return res

    def initial_params(self, x, y, xp, yp):
        """Initial parameters"""

        # Use center of image and half of y distance for x0, y0, and rl
        x0, y0, rl = [347.5, 259.5, 259.5]

        # appriximate lense function with line
        A, B, C, D = [np.pi / 2, -np.pi / 2, 0.0, 0.0]

        # calculate transformation with initial tranlation and lens function params but no rotation
        xu, yu = self.transform(x, y, x0, y0, rl, 0.0, A, B, C, D)

        # Find rotation matrix such that the vectors to the star locations roughly match
        Pu = np.array([xu, yu, np.zeros(len(x))]).T
        Pp = np.array([xp, yp, np.zeros(len(x))]).T
        R, _ = Rotation.align_vectors(Pp, Pu)

        # Find euler angles of rotation matrix and select "z" rotation as an approximate theta
        theta = R.as_euler("xyz", degrees=True)[2]

        return [x0, y0, rl, theta, C, D]

    def save_calibration_params(self, output, config_file):
        """Save results"""

        config = configparser.ConfigParser()
        config.read_string(config_file)

        config.set("CALIBRATION_PARAMS", "X0", str(self.x0))
        config.set("CALIBRATION_PARAMS", "Y0", str(self.y0))
        config.set("CALIBRATION_PARAMS", "RL", str(self.rl))
        config.set("CALIBRATION_PARAMS", "THETA", str(self.theta))
        config.set("CALIBRATION_PARAMS", "A", str(self.A))
        config.set("CALIBRATION_PARAMS", "B", str(self.B))
        config.set("CALIBRATION_PARAMS", "C", str(self.C))
        config.set("CALIBRATION_PARAMS", "D", str(self.D))

        with open(output, "w", encoding="utf-8") as cf:
            config.write(cf)


# -------------------------------------------------------------------------
# Main application
# -------------------------------------------------------------------------


def parse_args():
    """Command line options"""

    parser = argparse.ArgumentParser(description="Calculate camera calibration")

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")

    parser.add_argument(
        "-c", "--config", metavar="FILE", help="Alternate configuration file"
    )
    parser.add_argument(
        "-sc", "--starcal", metavar="FILE", help="Alternate starcal file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="mango-config.ini",
        help="Output configuration filename (default is mango-config.ini)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    return parser.parse_args()


def find_config(station, instrument):
    """Find configuration file from pacakge data"""

    config_file = f"{station}-{instrument}.ini"

    logging.debug("Using package configuration file: %s", config_file)

    return resources.files('mangonetwork.raw.data').joinpath(config_file).read_text()


def find_starcal(station, instrument):
    """FInd starcal file in pacakge data"""

    starcal_file = f"starcal-{station}-{instrument}.txt"

    logging.debug("Using package starcal file: %s", starcal_file)

    return resources.files("mangonetwork.raw.data").joinpath(starcal_file).read_text()


def main():
    """Main application"""

    args = parse_args()

    fmt = "[%(asctime)s] %(levelname)s %(message)s"

    if args.verbose:
        logging.basicConfig(format=fmt, level=logging.DEBUG)
    else:
        logging.basicConfig(format=fmt, level=logging.INFO)

    if args.config:
        logging.debug("Alternate configuration file: %s", args.config)
        if not os.path.exists(args.config):
            logging.error("Config file not found")
            sys.exit(1)
        with open(args.config, encoding="utf-8") as f:
            config_contents = f.read()
    else:
        config_contents = find_config(args.station, args.instrument)

    if args.starcal:
        logging.debug("Alternate starcal file: %s", args.starcal)
        if not os.path.exists(args.starcal):
            logging.error("StarCal file not found")
            sys.exit(1)
        with open(args.starcal, encoding="utf-8") as f:
            starcal_contents = f.read()
    else:
        starcal_contents = find_starcal(args.station, args.instrument)

    Calibrate(config_contents, starcal_contents, args.output)

    sys.exit(0)


if __name__ == "__main__":
    main()
