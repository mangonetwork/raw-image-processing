# data-processing
This repository contains the scripts needed for low-level MANGO data processing.

## Installation
This package is organized to the PEP 517/518 standards and can be installed with the latest version of pip.  To install for development:

1. Clone repository
2. Enter top directory
3. `pip install -e .`

## Usage
This package provides access to a number of command line programs for process raw MANGO data and generating standard data products:

- mango-calibrate
- mango-process-raw-images
- mango-keograms
- mango-quicklook-movies

Usage of each of these programs is described below.

### mango-calibrate
This produces the calibration file needed process data for each camera.  It requires a configuration file and the desired name of the output calibration file as command line options.  

**NOTE**: The calibration procedure has not been finalized, so the current version of this simply copies the old calibration arrays into the new, single-file format.  As developing the calibration routine, make sure it works with the raw data processing routines such that all parameters required to process the data are included in the calibration file.

### mango-process-raw-images
This produces a single processed data file from the provided raw image files where the output data are calibrated, equalized, rotated, and unwarpped.  It requires a configuration file which lists a calibration file (created by `mango-calibrate`) and several other parameters, a list of input raw files, and the name of the output file.  The list of input files can make use of wildcard symbols.

```
$ mango-process-raw-images -c <config file name> -i <list of input raw files> -o <output file name>
```
