#!/usr/bin/env python
"""Image Processing Core Engine for MANGO"""

##########################################################################
#
#   Image Processing Core Engine for MANGO
#   (Midlatitude All-sky-imager Network for Geophysical Observation)
#
#   2013-01-23  Rohit Mundra
#               Initial implementation
#
#   2013-02-05  Rohit Mundra
#               Robust Star Removal implemented
#               Histogram Equalization implemented
#
#   2013-02-22  Rohit Mundra
#               Changing unwarping to map to Normal Mercator Projection
#               Added alpha mask and PIL PNG saving
#
#   2013-02-22  Fabrice Kengne
#               Adding compatibility for FITS read
#
#   2014-10-25  Asti Bhatt
#               Adding compatibility for BU software generated .153 files
#
#   2023-03-07  Todd Valentic
#               Make PEP8 compliant
#               Remove unused code that had been commented out
#
#
#   Note: Most functionality in this file has been moved to imageops.py
#     This file can probably be depricated.
##########################################################################

import copy

import numpy as np
from PIL import Image
from scipy.interpolate import griddata
import skimage.transform

class MANGOImage:
    """Image manipulation functions"""


    def __init__(self, image_data):
        self.image_data = image_data

    def equalize(self, contrast, num_bins=10000):
        """Histogram Equalization to adjust contrast [1%-99%]"""

        # This function contains some questionable hard-coded things?

        image_array_1d = self.image_data.flatten()

        image_histogram, bins = np.histogram(image_array_1d, num_bins)
        image_histogram = image_histogram[1:]
        bins = bins[1:]
        cdf = np.cumsum(image_histogram)

        # spliced to cut off non-image area
        # any way to determine this dynamically?  How periminant is it?
        cdf = cdf[:9996]

        max_cdf = max(cdf)
        max_index = np.argmin(abs(cdf - contrast/100 * max_cdf))
        min_index = np.argmin(abs(cdf - (100 - contrast)/100 * max_cdf))
        vmax = float(bins[max_index])
        vmin = float(bins[min_index])
        low_value_indices = image_array_1d < vmin
        image_array_1d[low_value_indices] = vmin
        high_value_indices = image_array_1d > vmax
        image_array_1d[high_value_indices] = vmax
        self.image_data = image_array_1d.reshape(self.image_data.shape)

    def remove_stars(self, width, height):
        """Remove stars"""

        # function not currently used - revise later?

        filtered_data = copy.copy(self.image_data).astype(float)
        for i in range(width):
            left_pixels = self.image_data[:, max(i - 5, 0):max(i - 2, 1)]
            right_pixels = self.image_data[:, min(i + 3, width - 2):min(i + 6, width - 1)]
            all_pixels = np.append(left_pixels, right_pixels, 1)
            pixel_data = self.image_data[:, i]
            stdval = np.std(all_pixels, 1)
            mean_val = np.mean(all_pixels, 1)
            for j in range(height):
                if pixel_data[j] > mean_val[j] + 3 * stdval[j]:
                    filtered_data[j - 1:j + 2, i - 1:i + 2] = -1
                (i_known, j_known) = np.where(filtered_data >= 0)
            values_known = np.extract(filtered_data >= 0, filtered_data)

        # find minimum value from the image_data to find threshold
        m = np.empty(width)
        for i in range(width):
            m[i] = min(self.image_data[:, i])
        threshold = min(m)
        (i_all, j_all) = np.where(filtered_data >= threshold)  # min value from image_data
        self.image_data = np.reshape(
            griddata((i_known, j_known), values_known, (i_all, j_all), method='linear', fill_value=0), filtered_data.shape)

    def invert(self):
        """Flip image"""

        # Makes the under/above inversion more explicit/easier to understand

        self.image_data = np.flipud(self.image_data)

    def rotate(self, angle):
        """Rotate image"""

        # MUST flip image before rotating
        # Images are captured from below, but visualized from above in most standard formats
        self.image_data = skimage.transform.rotate(self.image_data, angle, order=3).astype(float)

    def transform_image(self, transport_coords, new_coords):
        """Interpolate the transformed coordinates to new (regular grid) coordinates."""

        # Interpolate the transformed coordinates to new (regular grid) coordinates.
        # This takes care of all traslation/rotation/unwarping required to create a
        # fully calibrated image. Also apply atmospheric corrections.

        img_corr = self.image_data
        xt_grid = transport_coords[0]
        yt_grid = transport_coords[1]

        x_grid = new_coords[0]
        y_grid = new_coords[1]

        interpolated_data = griddata((xt_grid.flatten(), yt_grid.flatten()), img_corr.flatten(), (x_grid, y_grid), fill_value=0)

        interpolated_data = (interpolated_data * 255 / (np.nanmax(interpolated_data)))
        self.image_data = interpolated_data.astype('uint8')


    def apply_mask(self, mask, fill_value=0):
        """mask the image by filling the mask region with a set value"""
        self.image_data[mask] = int(fill_value)

    def show_image(self):
        """function for debugging?"""
        img = Image.fromarray(self.image_data.astype('uint8'))
        img.show()
