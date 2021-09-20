import os
import ProcessImage
import time
import create_hdf5_files


if __name__ == '__main__':
    siteName = 'EIO'
    parentFolder = os.path.dirname(os.getcwd())
    rawFolder = "raw_data"
    rawPath = os.path.join(parentFolder, rawFolder)
    siteFolder = os.path.join(rawPath, siteName)
    list_of_dates = os.listdir(siteFolder)
    list_of_dates.remove('site_files')
    start_time_org = time.time()
    for i in list_of_dates:
        start_time = time.time()
        var = ProcessImage.ProcessImage(siteName, i)
        print('For image : ' + i + ' time taken = ', time.time() - start_time)
    print('For all images time taken = ', time.time() - start_time_org)
