import sys

import numpy as np
import PIL
import os
import h5py
import argparse
from datetime import datetime

if sys.version_info < (3,9):
    import importlib_resources as resources
else:
    from importlib import resources

import pandas as pd

class BITtoHDF5:

    def __init__(self, date, cf, bitfile, hdf5file):
        # inputs are config file, bitfile you want to convert, name of output hdf5 file
        self.date = date
        self.configFile = cf
        self.inputFile = bitfile
        self.outputFile = hdf5file

        # extract site, time

    def convert_to_PNG(self):
        f = open(self.inputFile, "rb")
        a = np.fromfile(f, dtype='int16')

        # is this format correct? we're re-writing as int32
        self.imageData1D = a[64:360769].astype('int32')
        self.imageData = self.imageData1D.reshape([519, 695])
        self.width = self.imageData.shape[0]
        self.height = self.imageData.shape[1]
        writeAddress = os.path.join(os.path.splitext(self.outputFile)[0], '.png')
        finalImage = PIL.Image.fromarray(self.imageData)

        finalImage.save(writeAddress, format='png')

    # the old sites are bridger, cfs, hat creek observatory, madison - is bridger brg?
    def write_to_HDF5(self):

        time_str = self.inputFile[1:-4]
        try:
            start_obj = datetime.strptime(self.date + time_str, '%Y%m%d%H%M%S')
        except ValueError:
            start_obj = datetime.strptime(self.date + time_str, '%b%d%Y%H%M%S')

        site_info = pd.read_csv(resources.files('mangonetwork.raw.data').joinpath('SiteInformation.csv'))
        site_code_letter = self.inputFile[0]
        row = site_info[site_info['Site Code']==site_code_letter]

        self.metadata = dict(
            # I think the time is already in UTC - check whether I have to make any conversions.
            start_time=start_obj.timestamp(),
            station=row['Site Abbreviation'].item(),
            latitude=row['Center Latitude'].item(),
            longitude=row['Center Longitude'].item(),

            serialnum="",
            device_name='Atik 414ex', # confirm
            label='Red Line',

            #i think we were only using greenline cameras
            instrument='redline',
            exposure_time="",
            x="",
            y="",
            width= 695,
            height= 519,
            bytes_per_pixel="", #recheck in
            bin_x="", # how did we bin this?
            bin_y="",
            ccd_temp="",
            set_point="",
            image_bytes=""
        )

        UnitsCatalog = dict(
            version='',
            start_time='Unix timetamp (UTC)',
            station='',
            latitude='degrees N',
            longitude='degrees E',
            serialnum='',
            device_name='',
            label='',
            instrument='',
            exposure_time='seconds',
            x='pixels',
            y='pixels',
            width='pixels',
            height='pixels',
            bytes_per_pixel='bytes',
            bin_x='pixels',
            bin_y='pixels',
            ccd_temp='degrees C',
            set_point='degrees C',
            image_bytes='bytes'
        )

        # artemis
        with h5py.File(self.outputFile, 'w') as output:
            image = output.create_dataset('image',
                                          data=self.imageData,
                                          compression='gzip')

            image.attrs.update(self.metadata)

            units = output.create_group('units')
            units.attrs.update(UnitsCatalog)

        return self.outputFile



def parse_args():
    """
    Handle the command line arguments.
    Returns:
    Output of argparse.ArgumentParser.parse_args.
    """

    parser = argparse.ArgumentParser(description='Accepting config, input BIT file and output HDF5 filename'
                                                 'to convert BIT to HDF5.')

    parser.add_argument('-d', '--date', dest = 'date', type=str,
                        help = 'Date in MMDDYYYY format.', required = True)
    parser.add_argument('-c', '--config', dest='config', type=str,
                        help='Config file containing ___.', required = True)
    parser.add_argument('-if', '--inputFile', dest='inputFile', nargs=1, type=str,
                        help='Input BIT file to be converted.')
    parser.add_argument('-il', '--inputList', dest='inputList', nargs = '+', type=str,
                        help='List of BIT files to be converted.')
    parser.add_argument('-o', '--output', dest='output', type=str,
                        help='Output HDF5 filename to write to.')
    args = parser.parse_args()

    return args


def main():
    command_line_args = parse_args()
    date = command_line_args.date
    conf = command_line_args.config
    output = command_line_args.output
    #what should the output be? I think this changes everything
    inputFile = command_line_args.inputFile

    # i've added a list component but idk how we would customize the output names
    # in that case, we might have to remove that argument
    BITtoHDF5(date, conf, inputFile, output)

if __name__ == '__main__':
    main()


