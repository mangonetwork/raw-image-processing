# data-processing
This repository contains the scripts needed for low-level MANGO data processing.

## Installation
This package is organized to the PEP 517/518 standards and can be installed with the latest version of pip.  
```
pip install git+https://github.com/mangonetwork/raw-image-processing.git
```
To install for development, clone then pip install with the `-e` flag.
```
git clone https://github.com/mangonetwork/raw-image-processing.git
cd raw-image-processing
pip install -e .
```

## Usage
This package provides command line programs for generating quicklook moves and processing raw images into level 1 data files (usually done nightly).  Both these programs require a configuration file that specifies metadata and some processing parameters.  The processing code requires the configuration file contain the `CALIBRATON_PARAMS` section, which is generated with the [image-calibration](https://github.com/mangonetwork/image-calibration) module.

### mango-process-raw-images
This produces a single processed data file from the provided raw image files where the output data are rotated and unwarpped.  
```
mango-process-raw-images <list of input raw files>
```
This program can take a long time to run.  Use the `-n` flag followed by an integer to customize the number of parallel processes that are run. Conversely, the `-s` flag will force the program to process images sequentially without using the multiprocessing module.  This can be useful for some debugging operations.


### mango-quicklook-movies
This produces a quicklook movie from provided raw image files.  The movie will be rotated and equalized, but not unwarpped.  It requires a configuration file which specifies at minimum a rotation angle (`manual_theta`).  If the configuration file has been calibrated, it will instead use the calibrated rotation angle.  
```
mango-quicklook-movies <list of input raw files>
```



Both programs have the options of using the following additional flags.

`-f` Instead of listing input files as command line arguments, they can be listed in a text file and that text file can be provided as input.  To generate a list of hdf5 files in a particular directory in the text file `filelist`:
```
ls path/to/data/directory > filelist
```
To generate a list of files for a single night if your directory strucure mimics that used in the MANGO database:
```
ls /path/to/local/mango/archive/<site>/<instrument>/raw/<year>/<doy>/??/*.hdf5 > filelist
```

`-c` Specify the config file to use rather than the default.

`-o` Specify the name of the output file instead of using the default.


## Collecting Data
This package operates on data posted in the [MANGO Database](https://data.mangonetwork.org/data/transport/mango/archive/).  Individual data files can be browsed and downloaded in a browser, but it may be useful to bulk download the days you're interested in with [wget](https://www.gnu.org/software/wget/).

Download a singe day's worth of data mimicing the directory structure of the online database:
```
wget -r -np -nH --cut-dir 8 -R "index.html*" https://data.mangonetwork.org/data/transport/mango/archive/<site>/<instrument>/raw/<year>/<doy>/
```
To only download hdf5 files, add the `-R "*.png"` flag.  To only download png files, add the `-R "*.hdf5"` flag.  To download to a specific local location, use the `-P` flag.
```
wget -r -np -nH --cut-dir 4 -R "index.html*" -P /path/to/local/archive https://data.mangonetwork.org/data/transport/mango/archive/<site>/<instrument>/raw/<year>/<doy>/
```
