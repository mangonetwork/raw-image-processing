# ==================================================================
# create_hdf5_files.py
# Convert individual MANGO *.png files to a single *.h5 file for a
#   given night with latitude/longitude arrays and site information
# Created by: L. Lamarche - 2019-3-14
# Based on MANGO_HDF_convert.py
# Updated:
#
# Notes:
#   - hard-coded file paths need to be updated
#   - requires glob - may not be Windows compatible?
#
# Requires config.ini, which specifies data file locations
# ; config.ini
#  
# [DEFAULT]
# SITEFILE = ???
# DATADIR = ???
# LATLONDIR = ???
#
# To Run:
#  > python create_hdf5_files.py [sites]
#  - the optional sites argument is a list of codes for sites you would
#     like to process
#  - if no sites are given, default is to process all sites
# e.g.
#  > python create_hdf5_files.py
#    will create *.h5 files for all sites in SiteInformation.csv
#  > python create_hdf5_files.py M R B
#    will create *.h5 fiels for Madison, Rainwater Observatory, and Bridger
# ==================================================================

from PIL import Image
import numpy as np
import glob
import regex as re
import h5py
import csv
import datetime as dt
import os
import configparser
import sys
import pandas as pd


# I modified the file so that it only works for one site

def hdf5_file_info(sitefile, processedImageDir, latlondir, site, dateFolder=None):
    # create site list from the site file and user input
    siteData = pd.read_csv(sitefile)
    site_list = siteData[siteData['Site Abbreviation'] == site]

    code = site_list['Site Abbreviation'].item()
    site_name = site_list['Site Name'].item()
    site_lon = site_list['Center Longitude'].item()
    site_lat = site_list['Center Latitude'].item()

    # read lat/lon from where ever Latitude.csv and Longitude.csv are for that site
    latDir = os.path.join(latlondir, 'Latitudes.csv')
    lonDir = os.path.join(latlondir, 'Longitudes.csv')
    try:
        latitude = np.array(pd.read_csv(latDir, dtype=float, delimiter=','))
        longitude = np.array(pd.read_csv(lonDir, dtype=float, delimiter=','))
        longitude[longitude < 0] += 360.
    except IOError:
        print('Could not process {}!'.format(site_name))
    # from MANGO_HDF_convert.py
    # sitepath = os.path.join(processedImageDir, site_name)
    sitepath = processedImageDir
    if not os.path.exists(sitepath):
        print('Could not find data directory for {}!'.format(site_name))

    def single_day_file(folderName):
        try:
            numbersOnly = re.findall(r'[0-9]+', folderName)[0]
            date = dt.datetime.strptime(numbersOnly, '%Y%m%d')
        except ValueError:
            date = dt.datetime.strptime(folderName, '%b%d%y')
        except:
            print('Date of folder is in an unexpected format.')

        datapath = os.path.join(processedImageDir, folderName)
        pimages = [os.path.basename(image) for image in glob.glob(os.path.join(datapath,'*.png'))]
        pimages.sort()
        images = []
        times = []
        filename = os.path.join(datapath, folderName + '.h5')
        if not os.path.isfile(filename):

            for pimageaddress in pimages:
                pimageaddress = os.path.join(datapath, pimageaddress)
                img = Image.open(pimageaddress)
                arrimg = np.array(img)
                images.append(arrimg[:, :, 0])
                timestring = pimageaddress[-11:-5]
                t = dt.datetime.strptime(timestring, '%H%M%S').replace(year=date.year, month=date.month,
                                                                       day=date.day)
                times.append(t)

            tstmp = np.array([(t - dt.datetime.utcfromtimestamp(0)).total_seconds() for t in times])
            images = np.array(images)

            # save hdf5 file
            f = h5py.File(filename, 'w')
            f.create_group('SiteInfo')
            dtype = h5py.special_dtype(vlen=str)
            N = f.create_dataset('SiteInfo/Name', dtype=dtype, data=site_name)
            N.attrs['Description'] = 'site name'
            C = f.create_dataset('SiteInfo/Code', data=code)
            C.attrs['Description'] = 'one letter site abbreviation/code'
            L = f.create_dataset('SiteInfo/Lon', data=site_lon)
            L.attrs['Description'] = 'geodetic longitude of site'
            L.attrs['Unit'] = 'degrees'
            L = f.create_dataset('SiteInfo/Lat', data=site_lat)
            L.attrs['Description'] = 'geodetic latitude of site'
            L.attrs['Unit'] = 'degrees'
            T = f.create_dataset('Time', data=tstmp, compression='gzip', compression_opts=1)
            T.attrs['Description'] = 'unix time stamp'
            T.attrs['Unit'] = 'seconds'
            T.attrs['Size'] = 'Nrecords'
            I = f.create_dataset('ImageData', data=images, compression='gzip', compression_opts=1)
            I.attrs['Description'] = 'pixel values for images'
            I.attrs['Size'] = 'Nrecords x Ipixels x Jpixels'
            Lon = f.create_dataset('Longitude', data=longitude, compression='gzip', compression_opts=1)
            Lon.attrs['Description'] = 'geodetic longitude of each pixel projected to 250 km'
            Lon.attrs['Size'] = 'Ipixels x Jpixels'
            Lon.attrs['Unit'] = 'degrees'
            Lat = f.create_dataset('Latitude', data=latitude, compression='gzip', compression_opts=1)
            Lat.attrs['Description'] = 'geodetic latitude of each pixel projected to 250 km'
            Lat.attrs['Size'] = 'Ipixels x Jpixels'
            Lat.attrs['Unit'] = 'degrees'

            f.close()

    if dateFolder is None:
        dates_in_processedImagesDir = next(os.walk(sitepath))[1]
        for procFolder in dates_in_processedImagesDir:
            single_day_file(procFolder)
    single_day_file(dateFolder)
