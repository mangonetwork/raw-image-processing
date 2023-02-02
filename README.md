# data-processing
This repository contains the scripts needed for low-level MANGO data processing.

## Installation
This package is organized to the PEP 517/518 standards and can be installed with the latest version of pip.  To install for development:

1. Clone repository
2. Enter top directory
3. `pip install -e .`

## Usage
This package provides access to a number of command line programs for process raw MANGO data and generating standard data products:

- mango-process-raw-images
- mango-calibrate
- mango-starcal
- mango-quicklook-movies
- mango-keograms

Usage of each of these programs is described below.

### mango-starcal
Validates a starcal file by overplotting the selected star points on the image.  This is provided mostly as a convenience and as a placeholder for future more sophisticated processes that could generate starcal files automatically.  The process requires the three-character station code and instrument (redline or greenline) as input.  NOTE: This process will likely not work on computers other than the one that was originally used to create the starcal file.

```
$ mango-starcal <station> <instrument> -sc <starcal file name>
```

### mango-calibrate
This adds the calibration parameter section to the provided config file using a specified a starcal file.  If a config file (`-c`) or a starcal file (`-sc`) are not provided, the default from the package data are used.  This program should only need to be run once for each camera.

```
$ mango-calibrate <station> <instrument> -c <config file name> -sc <starcal file name>
```

### mango-process-raw-images
This produces a single processed data file from the provided raw image files where the output data are calibrated, equalized, rotated, and unwarpped.  It requires a configuration file which lists the calibration parameters (created by `mango-calibrate`), metadata, and several other parameters specifying how the images are to be processed.  If a config file is not provide (`-c`), the default from the package data is used. This process also requires a file containing a list of input raw files (`-f`), and the name of the output file (`-o`).

```
$ mango-process-raw-images -c <config file name> -i <list of input raw files> -o <output file name>
```

### mango-quicklook-movies
This produces a quicklook movie from provided raw image files.  The movie will be rotated and equalized, but not unwarpped.  It requires a configuration file which specifies at minimum a rotation angle (`manual_theta`).  If the configuration file has been calibrated, it will instead use the calibrated rotation angle.  If a config file is not provide (`-c`), the default from the package data is used. This process also requires a file containing a list of input raw files (`-f`), and the name of the output movie file (`-o`).

```
$ mango-quicklook-movies -c <config file name> -i <list of input raw files> -o <output file name>
```

### mango-keograms
Placeholder for future keogram program.

## Configuration Files

This package utilizes two types of configuration files: config and starcal. Default configuration files are included as part of the package data under `/raw-image-processing/src/mangonetwork/raw/data`.

### starcal
The starcal files contain lists of stars and their coordinates in an image that can be used to calculate the calibration parameters by `mango-calibrate`.  The header of each file should include the full filename of the raw image used to identify star positions.  Columns in the file include the star name, real azimuth, real elevation, x position in image, and y position in image.  Presently, these files are created manually.  The `mango-starcal` program is helpful to create these files, but an "empty" file must first be created that specifies the raw image filename header

### config
The config file contains all information needed to run `mango-process-raw-images` and `mango-quicklook-movies`.  The `CALIBRATION_PARAMS` section is filled automatically by `mango-calibrate`.  The `PROCESSING` section contains site metadata and specific parameters needed to process images. The `QUICKLOOK` section contains a rotation angle entered manually so quicklook movies can be created before the camera is calibrated.
