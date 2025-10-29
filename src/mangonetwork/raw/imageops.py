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

def background_removal(image, x0, y0, rl):
    import matplotlib.pyplot as plt
    from scipy.stats import norm

    # Offset from edge of image and size of region for determining the background in each corner of image
    offx = int(0.015*image.shape[1])
    sizex = int(0.05*image.shape[1])
    offy = int(0.015*image.shape[0])
    sizey = int(0.05*image.shape[0])
    # This is the best we can do for uncalibrated
    # For calibrated, may be able to make a better estimate from the mask at the el=0 level

    xgrid, ygrid = np.meshgrid(np.arange(image.shape[1]), np.arange(image.shape[0]))
    mask = np.sqrt((xgrid-x0)**2 + (ygrid-y0)**2) > rl+10.

    masked_image = np.ma.masked_array(image, ~mask)
    c = plt.imshow(masked_image, vmin=0., vmax=600.)
    plt.axhline(offy)
    plt.axhline(offy+sizey)
    plt.axvline(offx)
    plt.axvline(offx+sizex)
    plt.colorbar(c)
    #plt.imshow(mask, cmap='Greys')
    plt.show()

    m1 = image[offy:offy+sizey, offx:offx+sizex].flatten()
    m2 = image[-(offy+sizey):-offy, -(offx+sizex):-offx].flatten()
    m3 = image[offy:offy+sizey, -(offx+sizex):-offx].flatten()
    m4 = image[-(offy+sizey):-offy, offx:offx+sizex].flatten()
    
    plt.hist(m1, bins=np.arange(0.,600.,10.), histtype='step', density=True)
    plt.hist(m2, bins=np.arange(0.,600.,10.), histtype='step', density=True)
    plt.hist(m3, bins=np.arange(0.,600.,10.), histtype='step', density=True)
    plt.hist(m4, bins=np.arange(0.,600.,10.), histtype='step', density=True)
    print(image.size, image[mask].size)
    mu = np.mean(image[mask])
    sig = np.std(image[mask])
    print(mu, sig)
    print(np.mean(m1), np.std(m1))
    print(np.mean(m2), np.std(m2))
    print(np.mean(m3), np.std(m3))
    print(np.mean(m4), np.std(m4))
    plt.hist(image[mask], bins=np.arange(0., 600., 10.), histtype='step', density=True)
    #plt.hist(image[mask], histtype='step', density=True)
    x = np.arange(0., 600., 10.)
    plt.plot(x, norm.pdf(x, loc=mu, scale=sig))
    plt.plot(x, norm.pdf(x, loc=np.mean(m1), scale=np.std(m1)))
    plt.show()

    # Calculate means in the four corners
    m1 = image[offy:offy+sizey, offx:offx+sizex].mean()
    m2 = image[-(offy+sizey):-offy, -(offx+sizex):-offx].mean()
    m3 = image[offy:offy+sizey, -(offx+sizex):-offx].mean()
    m4 = image[-(offy+sizey):-offy, offx:offx+sizex].mean()
    m = np.mean([m1,m2,m3,m4])
    # Subtract mean and set a 0 floor
    image = image - m
    image[image < 0] = 0

    return image

def estimate_background(image):

    # Offset from edge of image and size of region for determining the background in each corner of image
    offx = int(0.015*image.shape[1])
    sizex = int(0.05*image.shape[1])
    offy = int(0.015*image.shape[0])
    sizey = int(0.05*image.shape[0])

    # Calculate means in the four corners
    m1 = image[offy:offy+sizey, offx:offx+sizex].mean()
    m2 = image[-(offy+sizey):-offy, -(offx+sizex):-offx].mean()
    m3 = image[offy:offy+sizey, -(offx+sizex):-offx].mean()
    m4 = image[-(offy+sizey):-offy, offx:offx+sizex].mean()
    m = np.mean([m1,m2,m3,m4])
    ## Subtract mean and set a 0 floor
    #image = image - m
    #image[image < 0] = 0

    return m


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
