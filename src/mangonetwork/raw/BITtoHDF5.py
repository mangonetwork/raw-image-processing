import numpy as np
import PIL
import os
import h5py
import argparse

class BITtoHDF5:

    def __init__(self, cf, bitfile, hdf5file):
        # inputs are config file, bitfile you want to convert, name of output hdf5 file
        self.configFile = cf
        self.inputFile = bitfile
        self.outputFile = hdf5file
        #self.letter_to_site = {'C': 'Capitol Reef Field Station', 'B': 'Bridger', 'H'}
        self.letter_to_code = {'C': 'cfs', 'm': 'mad', 'h': 'hco'}


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

        self.metadata = dict(
            start_time=self.inputFile[1:-4],
            station=self.letter_to_code[self.inputFile[0]],

            #maybe store lon and lat in config file? it's already there in the csv though.
            # TODO
            latitude="",
            longitude="",

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
            image_bytes="",
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
    conf = command_line_args.config
    output = command_line_args.output
    #what should the output be? I think this changes everything
    inputFile = command_line_args.inputFile

    # i've added a list component but idk how we would customize the output names
    # in that case, we might have to remove that argument
    BITtoHDF5(conf, inputFile, output)

if __name__ == '__main__':
    main()


