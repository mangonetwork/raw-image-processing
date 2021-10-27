# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 12:56:24 2021

@author: padma
"""

#create keograms
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