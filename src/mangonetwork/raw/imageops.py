#!/usr/bin/env python
"""Image Processing Core Engine for MANGO"""

##########################################################################
#
#   Image Processing Functions
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
#               Refactor from class into individual methods
#
##########################################################################

import copy

import numpy as np
import skimage.transform


def equalize(image, contrast, num_bins=10000):
    """Histogram Equalization to adjust contrast [1%-99%]"""

    image_array_1d = image.flatten()

    image_histogram, bins = np.histogram(image_array_1d, num_bins)
    image_histogram = image_histogram[1:]
    bins = bins[1:]
    cdf = np.cumsum(image_histogram)

    # spliced to cut off non-image area
    # any way to determine this dynamically?  How periminant is it?
    cdf = cdf[:9996]

    max_cdf = max(cdf)
    max_index = np.argmin(abs(cdf - contrast / 100 * max_cdf))
    min_index = np.argmin(abs(cdf - (100 - contrast) / 100 * max_cdf))
    vmax = float(bins[max_index])
    vmin = float(bins[min_index])
    low_value_indices = image_array_1d < vmin
    image_array_1d[low_value_indices] = vmin
    high_value_indices = image_array_1d > vmax
    image_array_1d[high_value_indices] = vmax

    return image_array_1d.reshape(image.shape)


def invert(image):
    """Flip image"""

    # Makes the under/above inversion more explicit/easier to understand

    return np.flipud(image)


def rotate(image, angle):
    """Rotate image"""

    # MUST flip image before rotating
    # Images are captured from below, but visualized from above in most standard formats

    return skimage.transform.rotate(image, angle, order=3).astype(float)


def apply_mask(image, mask, fill_value=0):
    """mask the image by filling the mask region with a set value"""

    new_image = copy.copy(image)

    new_image[mask] = int(fill_value)

    return new_image
