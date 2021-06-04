##########################################################
#	Wrapper to batch process MANGO warped images 
#	and create unwarped images
#	This program needs to be in the same directory
#	as a successfully running MANGOImage.py 
#	Asti Bhatt Dec 5, 2014
##########################################################

import glob
import os
from StringIO import StringIO
from numpy import *
from PIL import Image
import numpy 
import datetime
from datetime import timedelta
import time
from time import mktime
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab
from matplotlib import dates
import MANGOimage
from MANGOimage import *



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
rawdirprefix= "/home/ftp/pub/earthcube/provider/asti/MANGO/"
siteDir = "Iowa"
#rawdirprefix = "/home/abhatt/workspace/WPI-SRI MQP - MANGO SYSTEM1/MANGO2/"

# directory where processed data will be 
#processeddirprefix = "/media/abhatt/Seagate Backup Plus Drive/workspace/Data/MANGO/InGeO/"
processeddirprefix = "/home/ftp/pub/earthcube/provider/asti/MANGOProcessed/"
#processeddirprefix = rawdirprefix + siteDir + "/"
#keogramdirprefix = "/media/abhatt/Seagate Backup Plus Drive/workspace/Data/MANGO/Keograms/"
keogramdirprefix = "/home/ftp/pub/earthcube/provider/asti/Keograms/"

# Process each folder and files within it
for rawFolder in next(os.walk(rawdirprefix + siteDir))[1]:
	rawPath = rawdirprefix + siteDir + "/" + rawFolder + "/"
	rawList = glob.glob1(rawPath, siteDir[0]+ '*.*')
	rawList.sort()
	print("Processing raw folder: " + rawFolder) 

# get Site name based on folder being processed
	if rawFolder[0][0] == 'H':
		siteName = 'Hat Creek Observatory'
	else:
		 if siteDir[0] == 'C':
			siteName = 'Capitol Reef Field Station'
		 else:
			if siteDir[0] == 'I':
				siteName = 'Eastern Iowa Observatory'
			else:
				if siteDir[0] == 'B':
					siteName = 'Bridger'

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
    
    #Make a keogram from the images just made
    	kfoldersdir = keogramdirprefix + siteDir + "/" + rawFolder
	if not os.path.isdir(kfoldersdir):
		os.makedirs(kfoldersdir)
    
		pimagelist = glob.glob1(pfoldersdir, siteDir[0]+'*.*')
		pimagelist.sort()
		# extract time from the filenames in the folder
		x = [time.strptime(rawFolder[-7:len(rawFolder)] + name[1:7], "%b%d%y%H%M%S") for name in pimagelist]  
		dt = [datetime.datetime.fromtimestamp(mktime(y)) for y in x]
		# create a 24-hour datetime array
		t1 = (datetime.datetime.fromtimestamp(mktime(time.strptime(rawFolder[-7:len(rawFolder)]+"000000","%b%d%y%H%M%S"))))
		#t2 = date2num(t1 + timedelta(days=1))	
		t2 = (datetime.datetime.fromtimestamp(mktime(time.strptime(rawFolder[-7:len(rawFolder)]+"235959","%b%d%y%H%M%S"))))
		d1 = dates.date2num(t1)
		d2 = dates.date2num(t2)
		# find the begin and end points for the actual data to fit in the 24-hour keogram 
		cad = dt[1]-dt[0]
		totalk = (t2-t1).total_seconds()/cad.total_seconds()
		begink = (dt[0]-t1).total_seconds()/cad.total_seconds()
		endk = (t2-dt[-1]).total_seconds()/cad.total_seconds()
		keoall = 999*numpy.ones([int(totalk),500])
		# create keogram from available data
		keo = numpy.empty([len(pimagelist),500])
		k=0;
		for img in pimagelist:
        		img1 = Image.open(pfoldersdir + "/" + img)
        		arr1 = numpy.array(img1)
	        	keo[k,:] = arr1[:,250,0]
        		k=k+1
        	# add the data to the empty array created above at right indices
		keoall[int(begink):int(begink)+len(pimagelist)] = keo   
    		# plot the keogram
    		fig = pylab.figure(num=None, figsize=(35, 2), dpi=100, facecolor='w', edgecolor='w')
		fig, ax = pylab.subplots()
		ax.imshow(keoall.transpose(), extent = [d1, d2, 0, 500], aspect ='auto',cmap='gray', clim=(0,250))
		ax.xaxis_date()
		ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
		date_formatter = dates.DateFormatter('%H%M')
		ax.xaxis.set_major_formatter(date_formatter)
		fig.savefig(kfoldersdir+"/"+rawFolder[-7:len(rawFolder)]+'_keogram.png',bbox_inches='tight',facecolor=fig.get_facecolor(), edgecolor='none') 
		print "Keogram saved"
    		
    		
    #create stacked images for every hour
    		
		# extract hours from the datetime object created above
	    	hrs = [(dt[k].hour) for k in range(len(dt))]
		# get counts for number of images in that hour
		b, counts = numpy.unique(hrs, return_counts=True)
		imc = 0
		# make stacked images for each hour
		for l in range(len(b)):
			# create a figure size proportional to number of subplots the hour will need
        		hdim = counts[l]*2+2
	        	fig1 = pylab.figure(num=None, figsize=(hdim, 2), dpi=100, facecolor='k', edgecolor='k')
        		for k in range(counts[l]):
        		    img2 = Image.open(pfoldersdir + "/" + pimagelist[imc])
	        	    fig = fig1.add_subplot(1,counts[l],k+1)
        		    fig.imshow(img2)
			# Turn off ticks from x and y axis
        		    fig.tick_params(axis='x', which='both', bottom='off', top='off', colors = 'white', labelbottom='off') 
	        	    fig.tick_params(axis='y', which='both', left='off', right='off', labelleft='off') 
			# Turn off frame
		            fig.set_frame_on(False) 
        		    fig.set_xlabel(pimagelist[imc][1:5], color='w')
	        	    imc = imc+1
        		    fig1.savefig(kfoldersdir+"/"+rawFolder[-7:len(rawFolder)]+'_hour'+str(b[l])+'.png',bbox_inches='tight',facecolor=fig1.get_facecolor(), edgecolor='none')
        	print "Stack images created"	
	pylab.clf() #clear all the plots made
    		
    
           

