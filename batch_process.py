import os
import ProcessImage


def date_formatter(originalName):
    dates = {
        'Jan': '01',
        'Feb': '02',
        'Mar': '03',
        'Apr': '04',
        'May': '05',
        'Jun': '06',
        'Jul': '07',
        'Aug': '08',
        'Sep': '09',
        'Oct': '10',
        'Nov': '11',
        'Dec': '12'
    }
    month = dates[originalName[:3]]
    day = originalName[3:5]
    year = '20' + originalName[5:]
    converted_str = year + month + day
    return converted_str


if __name__ == '__main__':
    siteName = 'EIO'
    parentFolder = os.path.dirname(os.getcwd())
    rawFolder = "raw_data"
    rawPath = os.path.join(parentFolder, rawFolder)
    siteFolder = os.path.join(rawPath, siteName)
    list_of_dates = os.listdir(siteFolder)
    for i in list_of_dates[:1]:
        var = ProcessImage.ProcessImage(siteName, i)
