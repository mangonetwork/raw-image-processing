# process_single_day.py


import MANGOimage

# Read the site information file to get 'PixArray'
siteInfofname = "SiteInformation.csv"
f = open(siteInfofname,'r')
siteText = f.read()
s = StringIO(siteText)
siteData = numpy.genfromtxt(s, dtype="|S", delimiter = ',', autostrip=True)

rows = siteData[:,0]

PixArray = []
for i in range(1,len(rows)):
    SiteInfo = siteData[i,:]
    PixArray.append(SiteInfo)

# directory where unprocessed data are
#rawdirprefix= "/media/abhatt/Seagate Backup Plus Drive/workspace/Data/MANGO/"
rawdirprefix= "/Users/e30737/Desktop/Data/MANGO/"
siteDir = "CRFS"
#rawdirprefix = "/home/abhatt/workspace/WPI-SRI MQP - MANGO SYSTEM1/MANGO2/"

# directory where processed data will be
#processeddirprefix = "/media/abhatt/Seagate Backup Plus Drive/workspace/Data/MANGO/InGeO/"
# processeddirprefix = "/home/ftp/pub/earthcube/provider/asti/MANGOProcessed/"
#processeddirprefix = rawdirprefix + siteDir + "/"
#keogramdirprefix = "/media/abhatt/Seagate Backup Plus Drive/workspace/Data/MANGO/Keograms/"
# keogramdirprefix = "/home/ftp/pub/earthcube/provider/asti/Keograms/"


# make directory based on folder being processed in the appropriate location
pfoldersdir = processeddirprefix + siteDir + "/" + rawFolder
diff = len(rawList)
if os.path.isdir(pfoldersdir):
	pList = glob.glob1(pfoldersdir, '*.*')
	pList.sort()
	diff = len(rawList)-len(pList)
else:
	os.makedirs(pfoldersdir)
#	print rawList[0]
# process individual images
for rawImage in rawList[len(rawList)-diff:len(rawList)]:
        #print rawImage
    	MANGOimage(siteName, siteDir, rawFolder, rawPath, rawImage, PixArray)
