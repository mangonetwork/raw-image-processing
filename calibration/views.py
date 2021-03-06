######################################################################
#
#       Written by:  Maria Rangel
#       Started:     02-06-2013       
#       Finished:    
#       Description: MANGO Website
#
######################################################################
from flask                          import Blueprint, request, render_template, flash, g, session, redirect, url_for
from flask                          import *
from werkzeug                       import secure_filename, check_password_hash, generate_password_hash
from app.users.decorators           import requires_login, adminrequires_login
from app                            import db
from app.users.models               import User
#from app.users.forms                import RegisterForm, LoginForm, ChangePassword, AdminRegisterForm, AdminLoginForm
from app.users.calibrationfunctions import lensfunction, pixelsToLatLon, greatCircleArc, earthDistanceFromZenithAngle
from numpy                          import *
from StringIO                       import StringIO
from math                           import atan, acos, asin, atan2
import matplotlib                   as plt
import matplotlib.pyplot            as plt
from pylab                          import *
from datetime                       import *
from calendar                       import *
from glob                           import glob
#from flask.ext.babel                 import Babel
import gc
import os
import csv
import sqlite3
import glob

mod = Blueprint('users', __name__, url_prefix='')

# SYSTEM VARIABLES
SITES_DIR = 'app/static/Sites/'
CAL_DIR = '/calibration/'
SYS_DIR = '/system/'
SCHEDULE_DIR = '/schedule/'
IMG_DIR = '/images/'
ALL_IMG_DIR = '/static/Sites/All Sites Images/'
#ALL_IMG_DIR = '/home/abhatt/workspace/WPI-SRI MQP - MANGO SYSTEM1/SERVER/myproject/MANGO/app/static/Sites/All Sites Images/'

#############################################################
#
#       Functions that aid the functionality of MANGO
#
#############################################################


# Opens the addsite CSV file and creates a dictionary for sites to be used
# in multiple locations throughout the MANGO pages
def sitesInNetwork():
    f = open('app/static/Sites/addsite.csv', 'r')
    text = f.read()
    s = StringIO(text)
    data = genfromtxt(s, dtype=None, delimiter = ',')
    site=[]

    try:
        for i in range(1,len(data)):
            temp= dict(id=data[i,0], name=data[i,1].replace('"',''),lat=data[i,2],lon=data[i,3], elevation=data[i,4], deleted=data[i,5])
            site.append(temp)
    except:
        site = "Empty"

    f.close()
    
    return site


# This function reads the schedule File and extracts the information in the CSV file and returns a dictionary with all the inforamtion
def scheduleFileDataExtraction(siteName):
    scheduleDataExtraction = []                          # Creating an empty array for dictionary usage for returning variables

    # Schedule File data extraction for strings for the schedule data extraction
    scheduleFilePath = SITES_DIR + siteName + SCHEDULE_DIR +  'Schedule.csv'   # path of location for the schedule
    scheduleFileUser = open( scheduleFilePath, 'r')      # opens the schedule CSV file
    text1 = scheduleFileUser.read()                      # reads the CSV file 
    p = StringIO(text1)                                  # reads the strings within the CSV file
    data2 = genfromtxt(p, dtype=None, delimiter = ',')   # extracts the data from CSV file 
    scheduleFileUser.close()                             # closes the CSV file       
    userData = data2[:,1]

    # Assigning variables
    adminName = ''.join(map(str, data2[1][1]))  # extracts the administrator's name from the schedule
    date = ''.join(map(str, data2[2][1]))       # extracts the date from the schedule
    time = ''.join(map(str, data2[3][1]))       # extracts the time from the schedule

    # Schedule File data extraction for intengers for the schedule data extraction
    scheduleFile = open(SITES_DIR + siteName + SCHEDULE_DIR +  'Schedule.csv' , 'r')   # opens the schedule CSV file
    text1 = scheduleFile.read()                 # reads the CSV file 
    p = StringIO(text1)                         # reads the strings within the CSV file
    data1 = genfromtxt(p, delimiter = ',')      # extracts the data from CSV file 
    scheduleFile.close()                        # closes the CSV file       

    scheduleData = data1[:,1]                   # extracting the information extracted by the CSV reader

    exposureTime = int(scheduleData[6])         # extracting the Stop Time from the data obtained from CSV file 
    exposureSeconds = exposureTime%100
    exposureMinutes= int(exposureTime/100)

    restPeriodTime = int(scheduleData[7])       # extract the Rest Period Time from the data obtained from the CSV file
    restPeriodTimeSeconds = restPeriodTime      # renamed to match the rest of the format of the form in HTML

    binning = int(scheduleData[8])              # extract the Binnign from the data obtained from the CSV file
    binX = binning
    binY = binX                                 # X binnging and Y binning will be equal to each other

    targetTemperature = int(scheduleData[10])   # extract the Target Temperature from the data obtained from the CSV file

    preCoolingDuration = int(scheduleData[12])  # extract the Target Temperature from the data obtained from the CSV file
    
    if preCoolingDuration == 0:                 # setting Pre-Cooling Duration 
        preCooling = 0                          # turns OFF Pre-Cooling
        preCoolingMinutes = preCoolingDuration  # renamed to match the rest of the format of the form in HTML

    if preCoolingDuration > 0: 
        preCooling = 1                          # turns ON Pre-Cooling
        preCoolingMinutes = preCoolingDuration  # renamed to match the rest of the format of the form in HTML
   
   # Creating a dictionary with all the variables extracted from the CSV for returning in the function
    scheduleDataExtraction = dict(adminName=adminName, date=date, time=time, exposureSeconds=exposureSeconds, exposureMinutes=exposureMinutes, restPeriodTimeSeconds=restPeriodTimeSeconds, binX=binX, binY=binY, targetTemperature=targetTemperature, 
        preCoolingMinutes=preCoolingMinutes, preCoolingDuration=preCoolingDuration)
    
    return scheduleDataExtraction

# This function extracts all the data from the CSV necessary for the graphs in the site monitoring page
def graphData(siteName, filename):
     # Opening file to obtain Floats
    FilePath = SITES_DIR + siteName + SYS_DIR + filename      # File Path
    
    fileOpen = open(FilePath,  'r')                           # Open File and in read mode
    text = fileOpen.read()
    s = StringIO(text)
    data = genfromtxt(s, dtype= float, delimiter = ',')       # extracts the data in float mode
    yAxis = array(([]))
    yAxis = data[:,1]
    fileOpen.close()
    
    numOfElements = len(data)                                 # extracts the amount of elements in the data set
    
    # Opening CSV file to open Floats
    fileOpen = open(FilePath,  'r')                                   # Open File and in read mode
    text = fileOpen.read()
    s = StringIO(text)
    data2 = genfromtxt(s, dtype=str, delimiter = ',')                   # extracts the data in strings format
    xAxis = map(str, data2[:,0])
    fileOpen.close()

    return xAxis, yAxis, numOfElements

# This function modifies the data to be passed out the HTML and java script. It modifies the time by making a timedate object. 
def modifyData(day, time, elements, yAxis):
    modifiedYAxis = []
    modifiedTime = []
    for i in range(elements):
        temp = time[i].split('.')            # splits seconds from miliseconds 
        date_object = datetime.strptime(temp[0], '%Y-%m-%d %H:%M:%S')
        currentTime = datetime.today()
        if currentTime - timedelta(days=day) < date_object:
            modifiedTime.append(date_object)
            modifiedYAxis.append(yAxis[i])
    numOfElements = len(modifiedTime)
    return modifiedTime, modifiedYAxis, numOfElements

# This function modifies the data from the abiemt time to be passed out the HTML and java script. It modifies the time by making a timedate object. 
def modifyAmbient(day, time, elements, yAxis):
    modifiedYAxis = []
    modifiedTime = []
    for i in range(elements): 
        date_object = datetime.strptime(time[i], '%m/%d/%Y %I:%M:%S %p')
        currentTime = datetime.today()
        if currentTime - timedelta(days=day) < date_object:
            modifiedTime.append(date_object)
            modifiedYAxis.append(yAxis[i])
    numOfElements = len(modifiedTime)
    return modifiedTime, modifiedYAxis, numOfElements

#############################################################
#
#               Public Webpage Control
#
#############################################################
# pulls the users profile from the database
@mod.before_request
def before_request():
    g.role = None
    g.user = None
    if 'user_id' in session:
        g.user = User.query.get(session['user_id'])


#Home page Map- this page displays images on the map
@mod.route('/', methods = ['GET', 'POST'])
def home():
    # Creating empty arrays 
    sitesAvailable = []                     # site available array
    imageAddress = []                       # storing the image address
    imageInformation =[]                    # storing image information
    site = sitesInNetwork()                 # calls the function to find the sites in the network using the addsite.csv file
    
    for sites in site: 
        if sites['deleted'] == '0':                 # Deleted column must equal to '0' to be identified on the network
            sitesAvailable.append(sites['name'])

    amountOfSites = len(sitesAvailable)             # The amount of sites currently "active" on the network
    print sitesAvailable
    
    for i in range(amountOfSites):
            siteName = sitesAvailable[i]  
	    print siteName      
            imageAddress = ALL_IMG_DIR + siteName + '.png'
#	    imageAddress = glob.glob(ALL_IMG_DIR + siteName[0] + '*.png')[0]
	    
            
            filename = SITES_DIR + siteName + CAL_DIR + 'imageLatLong.csv' 
            
            f = open(filename,'r')
            text = f.read()
            s = StringIO(text)
            data = genfromtxt(s, dtype=float, delimiter = ',')

            southWestLatLon = data[0,:]
            northEastLatLon = data[1,:]
            #print southWestLatLon,northEastLatLon
            tmpAddress = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))+imageAddress
#	    tmpAddress = imageAddress
            
            modifiedTime = datetime.fromtimestamp(os.path.getmtime(tmpAddress))

            currentTime = datetime.now()

            if (currentTime - modifiedTime) < timedelta(minutes=10):
                temp = dict(siteName=siteName, imageAddress=imageAddress, southWestLatLon=southWestLatLon, northEastLatLon=northEastLatLon)
                imageInformation.append(temp)
            else:
                temp = dict(siteName=siteName, imageAddress=imageAddress, southWestLatLon=southWestLatLon, northEastLatLon=northEastLatLon)
                imageInformation.append(temp)
                pass

    return render_template('users/Public/home.html', imageInformation=imageInformation, amountOfSites=amountOfSites)

#Info Public page
@mod.route('/info/')
def info():
    return render_template('users/Public/info.html')

#Info Public page
@mod.route('/ionosphere/')
def ionosphere():
    return render_template('users/Public/ionosphere.html')

#Space Phenomena page
@mod.route('/spacephenomena/')
def spacephenomena():
    return render_template('users/Public/spacephenomena.html')

# Register Public page
@mod.route('/register/', methods = ['GET', 'POST'])
def register():
    form = RegisterForm(request.form)
    if form.validate_on_submit():
        # create an user instance this information is not in the datbase
        user = User(form.firstname.data, form.lastname.data, form.affiliation.data, form.username.data, form.email.data, generate_password_hash(form.password.data), role=1)        
        
        try: 
            username_exists =db.session.query(User.username).filter_by(username=form.username.data).one()
            if username_exists:
                flash('Username is already taken. Please select a different username.')
                return render_template("users/Public/register.html", form=form)
        except:
            pass
        
        try:
            email_exists = db.session.query(User.email).filter_by(email=form.email.data).one()
            if email_exists:
                    flash('Email not unique. Email already registered to the site.')
                    return render_template("users/Public/register.html", form=form)
        except:
            pass
     


        db.session.add(user)
        db.session.commit()

        session['user_id'] = user.id
        session['user_role'] = user.role
        session['user_email'] = user.email

        flash('User Registered')
        return redirect(url_for("users.home"))

    return render_template('users/Public/register.html', form=form)

# Login Public Page   
@mod.route('/login/', methods = ['GET', 'POST'])
def login():
    form = LoginForm(request.form)
    # make sure all data is valid 
    # Does not validate password is right
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        # we use werzeug to validate user's password
        if user and check_password_hash(user.password, form.password.data):
            # the session can't be modified as it's signed, 
            # it's a safe place to store the user id
            session['user_id'] = user.id
            session['user_role'] = user.role
            session['user_email'] = user.email

            flash('Welcome,')
            flash('%s' % user.username)
            return redirect(url_for('users.home'))
        flash('Wrong email or password', 'error-message')

    
    return render_template('users/Public/login.html', form=form)

#Public Logout page
@mod.route('/logout/')
def logout():
    session.pop('user_id',None)
    flash("Logged out.")
    return render_template('users/Public/logout.html')

#Logout page
@mod.route('/comment/')
def comment():
    return render_template('users/Public/comment.html')


#############################################################
#
#           Administrators Webpage Control
#
#############################################################
# MANGO System Main Page
@mod.route('/MANGO/', methods = ['GET', 'POST'])
def MANGO():
    if "submit" in request.form:
        return redirect(url_for('users.adminlogin'))

    return render_template('users/Private/MANGO.html')

#Admin Login Page
@mod.route('/adminlogin/', methods = ['GET', 'POST'])
def adminlogin():
    form = LoginForm(request.form)
     # make sure all data is valid 
     # Does not validate password is right
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        # we use werzeug to validate user's password
        if user and check_password_hash(user.password, form.password.data):
            # the session can't be modified as it's signed, 
            # it's a safe place to store the user id
            session['user_id'] = user.id
            session['user_role'] = user.role
            session['user_email'] = user.email

            flash('Welcome,')
            flash('%s' % user.username)
            return redirect(url_for('users.main'))
        flash('Wrong email or password', 'error-message')

    return render_template('users/Private/adminlogin.html', form=form)

#Admin Logout page
@mod.route('/logout/')
def adminlogout():
    session.pop('user_id',None)
    redirect(url_for('users.MANGO'))
    
    return render_template('users/Private/adminlogout.html')

# Register Private Page
@mod.route('/adminregister/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def adminregister():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc'))
    
    form = RegisterForm(request.form)
    if form.validate_on_submit():
        # create an user instance this information is not in the datbase
        user = User(form.firstname.data, form.lastname.data, form.affiliation.data, form.username.data, form.email.data, generate_password_hash(form.password.data), role=0)        
        
        try: 
            username_exists =db.session.query(User.username).filter_by(username=form.username.data).one()
            if username_exists:
                flash('Username is already taken. Please select a different username.')
                return render_template("users/Public/register.html", form=form)
        except:
            pass
        
        try:
            email_exists = db.session.query(User.email).filter_by(email=form.email.data).one()
            if email_exists:
                    flash('Email not unique. Email already registered to the site.')
                    return render_template("users/Public/register.html", form=form)
        except:
            pass
     
        db.session.add(user)
        db.session.commit()

        session['user_id'] = user.id
        session['user_role'] = user.role
        session['user_email'] = user.email

        flash('User Registered')

        return redirect(url_for("users.main"))

    return render_template('users/Private/adminregister.html', form=form, site=site)

#Admin Home page
@mod.route('/account/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def account():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc'))
    
    if "changePassword" in request.form:
        return redirect(url_for("users.changepassword"))

    if "register" in request.form:
        return redirect(url_for("users.adminregister"))


    return render_template('users/Private/account.html',site=site)

#Change Password page
@mod.route('/changepassword/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def changepassword():
    site = sitesInNetwork()
    for sites in site:
        if 'mysite'+sites['id'] in request.form:
            session['current_id'] = sites['id']
            session['current_site'] = sites['name']
            return redirect(url_for('users.site_desc'))

    form = ChangePassword(request.form)
    if form.validate_on_submit():
        user = User.query.filter_by(id=session['user_id']).first()
        if check_password_hash(user.password, form.current_password.data): 
            user.password = generate_password_hash(form.new_password.data)
            db.session.commit()
            flash('Password successfully changed.')
        else:
            flash('Incorrect password')

    return render_template('users/Private/changepassword.html',site=site, form=form)


#############################################################
#
#           MANGO System Webpage Control
#
#############################################################

#Admin Home page
@mod.route('/main/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def main():
    
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    
    user = User.query.filter_by(id=session['user_id']).first()
    firstName = user.firstname
   
    return render_template('users/Private/main.html',site=site, firstName=firstName)


#############################################################
#
#                       Sites Webpage Control
#
#############################################################

#Site Monitoring System page
@mod.route('/site_desc/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def site_desc():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

    siteName = session['current_site']          # Variable used to extract the site name 

    # Ambient Temperature Graph empty variables 
    ambientTemperature = []     # ambient temperature
    ambientTime = []            # ambien time
    numOfAmbientElements = []   
    modifiedAmbientTime = []
    ambientTime, ambientTemperature, numOfAmbientElements = graphData(siteName, 'ambient.csv')     # Extarcting data from the CSV files for ambient temp

    # CCD Temperature Graph
    ccdTemperature = []
    ccdTime = []
    numOfCCDElements = []
    modifiedCCDTime =[]
    modifiedCCDTemp = []
    ccdTime, ccdTemperature, numOfCCDElements = graphData(siteName, 'ccd.csv')                      # Extarcting data from the CSV files for ccd temp
    

     # Disk Space Graph
    diskSpace = []
    diskTime = []
    numOfDiskElements = []
    modifiedDiskTime =[]
    diskTime, diskSpace, numOfDiskElements = graphData(siteName, 'disk.csv')                        # Extarcting data from the CSV files for disk space

    # Battery
    batteryLevel = []
    batteryTime = []
    numOfBatteryElements = []
    modifiedBatteryTime =[]
    batteryTime, batteryLevel, numOfBatteryElements = graphData(siteName, 'battery.csv')            # Extarcting data from the CSV files for battery level


    # day in the modifydata() and modifyambient() functions is set to 1 and this function modifies the time extracted from the CSV files to make them datetime objects
    modifiedAmbientTime, modifiedAmbientTemp, numOfAmbientElements = modifyAmbient(1,ambientTime,numOfAmbientElements,ambientTemperature)
    modifiedCCDTime, modifiedCCDTemp, numOfCCDElements = modifyData (1,ccdTime,numOfCCDElements,ccdTemperature)
    modifiedBatteryTime, modifiedBatteryLevel, numOfBatteryElements = modifyData (1,batteryTime,numOfBatteryElements,batteryLevel)
    modifiedDiskTime, modifiedDiskSpace, numOfDiskElements = modifyData (1,diskTime,numOfDiskElements,diskSpace)
    

    # If the user wishes to view the last 24 hours worth of data. 
    if "1Day" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

        siteName = session['current_site']          # Variable used to extract the site name 

        # Extarcting the necessary information from the CSV files for ambient temp, ccd temp, disk space, and battery level
        ambientTime, ambientTemperature, numOfAmbientElements = graphData(siteName, 'ambient.csv')
        ccdTime, ccdTemperature, numOfCCDElements = graphData(siteName, 'ccd.csv')
        diskTime, diskSpace, numOfDiskElements = graphData(siteName, 'disk.csv')
        batteryTime, batteryLevel, numOfBatteryElements = graphData(siteName, 'battery.csv')


        # day in the modifydata() and modifyambient() functions is set to 1 since 24 hours = 1 day
        modifiedAmbientTime, modifiedAmbientTemp, numOfAmbientElements = modifyAmbient(1,ambientTime,numOfAmbientElements,ambientTemperature)
        modifiedCCDTime, modifiedCCDTemp, numOfCCDElements = modifyData (1,ccdTime,numOfCCDElements,ccdTemperature)
        modifiedBatteryTime, modifiedBatteryLevel, numOfBatteryElements = modifyData (1,batteryTime,numOfBatteryElements,batteryLevel)
        modifiedDiskTime, modifiedDiskSpace, numOfDiskElements = modifyData (1,diskTime,numOfDiskElements,diskSpace)
    

        return render_template('users/Private/site_desc.html', site=site, mysite = session['current_site'], siteName=siteName, modifiedAmbientTime=modifiedAmbientTime, modifiedAmbientTemp=modifiedAmbientTemp, numOfAmbientElements=numOfAmbientElements,
            modifiedCCDTime=modifiedCCDTime, modifiedCCDTemp=modifiedCCDTemp, numOfCCDElements=numOfCCDElements, modifiedBatteryTime=modifiedBatteryTime, modifiedBatteryLevel=modifiedBatteryLevel, numOfBatteryElements=numOfBatteryElements, modifiedDiskTime=modifiedDiskTime, 
            modifiedDiskSpace=modifiedDiskSpace, numOfDiskElements=numOfDiskElements)

    # If the user wishes to view the last 5 days worth of data. 
    if "5Day" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

        siteName = session['current_site']          # Variable used to extract the site name 

        # Extracting the necessary information from the CSV files for ambient temp, ccd temp, disk space, and battery level
        ambientTime, ambientTemperature, numOfAmbientElements = graphData(siteName, 'ambient.csv')
        ccdTime, ccdTemperature, numOfCCDElements = graphData(siteName, 'ccd.csv')
        diskTime, diskSpace, numOfDiskElements = graphData(siteName, 'disk.csv')
        batteryTime, batteryLevel, numOfBatteryElements = graphData(siteName, 'battery.csv')

        # day in the modifydata() and modifyambient() functions is set to 5 since the user wishes to view 5 days worth of data
        modifiedAmbientTime, modifiedAmbientTemp, numOfAmbientElements = modifyAmbient(5,ambientTime,numOfAmbientElements,ambientTemperature)
        modifiedCCDTime, modifiedCCDTemp, numOfCCDElements = modifyData (5,ccdTime,numOfCCDElements,ccdTemperature)
        modifiedBatteryTime, modifiedBatteryLevel, numOfBatteryElements = modifyData (5,batteryTime,numOfBatteryElements,batteryLevel)
        modifiedDiskTime, modifiedDiskSpace, numOfDiskElements = modifyData (5,diskTime,numOfDiskElements,diskSpace)
    

        return render_template('users/Private/site_desc.html', site=site, mysite = session['current_site'], siteName=siteName, modifiedAmbientTime=modifiedAmbientTime, modifiedAmbientTemp=modifiedAmbientTemp, numOfAmbientElements=numOfAmbientElements,
            modifiedCCDTime=modifiedCCDTime, modifiedCCDTemp=modifiedCCDTemp, numOfCCDElements=numOfCCDElements, modifiedBatteryTime=modifiedBatteryTime, modifiedBatteryLevel=modifiedBatteryLevel, numOfBatteryElements=numOfBatteryElements, modifiedDiskTime=modifiedDiskTime, 
            modifiedDiskSpace=modifiedDiskSpace, numOfDiskElements=numOfDiskElements)

    # If the user wishes to view the last month worth of data. 
    if "1Month" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

        siteName = session['current_site']          # Variable used to extract the site name 

        # Extracting the necessary information from the CSV files for ambient temp, ccd temp, disk space, and battery level
        ambientTime, ambientTemperature, numOfAmbientElements = graphData(siteName, 'ambient.csv')
        ccdTime, ccdTemperature, numOfCCDElements = graphData(siteName, 'ccd.csv')
        diskTime, diskSpace, numOfDiskElements = graphData(siteName, 'disk.csv')
        batteryTime, batteryLevel, numOfBatteryElements = graphData(siteName, 'battery.csv')

        currentTime = datetime.today()                                      # obtains the current time of the system
        daysInMonth = monthrange(currentTime.year, currentTime.month)       # obtains the amount of days in the month using the current time year and month


        # day in the modifydata() and modifyambient() functions is set to days in month depending on the day of the month since the user wishes to view data from
        # the last month from today's date
        modifiedAmbientTime, modifiedAmbientTemp, numOfAmbientElements = modifyAmbient(daysInMonth[1],ambientTime,numOfAmbientElements,ambientTemperature)
        modifiedCCDTime, modifiedCCDTemp, numOfCCDElements = modifyData (daysInMonth[1],ccdTime,numOfCCDElements,ccdTemperature)
        modifiedBatteryTime, modifiedBatteryLevel, numOfBatteryElements = modifyData (daysInMonth[1],batteryTime,numOfBatteryElements,batteryLevel)
        modifiedDiskTime, modifiedDiskSpace, numOfDiskElements = modifyData (daysInMonth[1],diskTime,numOfDiskElements,diskSpace)
    

        return render_template('users/Private/site_desc.html', site=site, mysite = session['current_site'], siteName=siteName, modifiedAmbientTime=modifiedAmbientTime, modifiedAmbientTemp=modifiedAmbientTemp, numOfAmbientElements=numOfAmbientElements,
            modifiedCCDTime=modifiedCCDTime, modifiedCCDTemp=modifiedCCDTemp, numOfCCDElements=numOfCCDElements, modifiedBatteryTime=modifiedBatteryTime, modifiedBatteryLevel=modifiedBatteryLevel, numOfBatteryElements=numOfBatteryElements, modifiedDiskTime=modifiedDiskTime, 
            modifiedDiskSpace=modifiedDiskSpace, numOfDiskElements=numOfDiskElements)

    # If the user wishes to view the last six months worth of data. 
    if "6Month" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        siteName = session['current_site']          # Variable used to extract the site name 
        
        # Extracting the necessary information from the CSV files for ambient temp, ccd temp, disk space, and battery level
        ambientTime, ambientTemperature, numOfAmbientElements = graphData(siteName, 'ambient.csv')
        ccdTime, ccdTemperature, numOfCCDElements = graphData(siteName, 'ccd.csv')
        diskTime, diskSpace, numOfDiskElements = graphData(siteName, 'disk.csv')
        batteryTime, batteryLevel, numOfBatteryElements = graphData(siteName, 'battery.csv')

        currentTime = datetime.today()      # obtains the current time of the system
        totalDays = 0                       # intiating a count variable
       
        # A loop to obtain the total amount of days from today's current time.  
        for i in range(6): 
            month = currentTime.month - i           # sets month to the current month minus the variable of the loop 
            year = currentTime.year                 # extracts the current year
            if month <= 0:                          # if the loop reaches january and it needs to go to decemeber by reseting to 12 and subtracting from there
                month = 12 + month
                year -= 1
            daysInMonth = monthrange(year, month)  
            totalDays+= daysInMonth[1]                  #Adds all days from the loop

        # day in the modifydata() and modifyambient() functions is set to "totalDays" depending on the day of the month since the user wishes to view data from
        # the lastsix months from today's date
        modifiedAmbientTime, modifiedAmbientTemp, numOfAmbientElements = modifyAmbient(totalDays,ambientTime,numOfAmbientElements,ambientTemperature)
        modifiedCCDTime, modifiedCCDTemp, numOfCCDElements = modifyData (totalDays,ccdTime,numOfCCDElements,ccdTemperature)
        modifiedBatteryTime, modifiedBatteryLevel, numOfBatteryElements = modifyData (totalDays,batteryTime,numOfBatteryElements,batteryLevel)
        modifiedDiskTime, modifiedDiskSpace, numOfDiskElements = modifyData (totalDays,diskTime,numOfDiskElements,diskSpace)
    

        return render_template('users/Private/site_desc.html', site=site, mysite = session['current_site'], siteName=siteName, modifiedAmbientTime=modifiedAmbientTime, modifiedAmbientTemp=modifiedAmbientTemp, numOfAmbientElements=numOfAmbientElements,
            modifiedCCDTime=modifiedCCDTime, modifiedCCDTemp=modifiedCCDTemp, numOfCCDElements=numOfCCDElements, modifiedBatteryTime=modifiedBatteryTime, modifiedBatteryLevel=modifiedBatteryLevel, numOfBatteryElements=numOfBatteryElements, modifiedDiskTime=modifiedDiskTime, 
            modifiedDiskSpace=modifiedDiskSpace, numOfDiskElements=numOfDiskElements)

    # If the user wishes to view the last year worth of data. 
    if "1Year" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

        siteName = session['current_site']          # Variable used to extract the site name 

        # Extracting the necessary information from the CSV files for ambient temp, ccd temp, disk space, and battery level
        ambientTime, ambientTemperature, numOfAmbientElements = graphData(siteName, 'ambient.csv')
        ccdTime, ccdTemperature, numOfCCDElements = graphData(siteName, 'ccd.csv')
        diskTime, diskSpace, numOfDiskElements = graphData(siteName, 'disk.csv')
        batteryTime, batteryLevel, numOfBatteryElements = graphData(siteName, 'battery.csv')

        currentTime = datetime.today()      # obtains the current time of the system
        totalDays = 0                       # intiating a count variable
        
         # A loop to obtain the total amount of days from today's current time.  
        for i in range(12):                         # range of 12 because thats the amount of months in a year
            month = currentTime.month - i           # sets month to the current month minus the variable of the loop 
            year = currentTime.year                 # extracts the current year
            if month <= 0:                          # if the loop reaches january and it needs to go to decemeber by reseting to 12 and subtracting from there
                month = 12 + month
                year -= 1
            daysInMonth = monthrange(year, month)  
            totalDays+= daysInMonth[1]                  #Adds all days from the loop
        
        # day in the modifydata() and modifyambient() functions is set to "totalDays" depending on the day of the month since the user wishes to view data from
        # the twelve months from today's date
        modifiedAmbientTime, modifiedAmbientTemp, numOfAmbientElements = modifyAmbient(totalDays,ambientTime,numOfAmbientElements,ambientTemperature)
        modifiedCCDTime, modifiedCCDTemp, numOfCCDElements = modifyData (totalDays,ccdTime,numOfCCDElements,ccdTemperature)
        modifiedBatteryTime, modifiedBatteryLevel, numOfBatteryElements = modifyData (totalDays,batteryTime,numOfBatteryElements,batteryLevel)
        modifiedDiskTime, modifiedDiskSpace, numOfDiskElements = modifyData (totalDays,diskTime,numOfDiskElements,diskSpace)
    

        return render_template('users/Private/site_desc.html', site=site, mysite = session['current_site'], siteName=siteName, modifiedAmbientTime=modifiedAmbientTime, modifiedAmbientTemp=modifiedAmbientTemp, numOfAmbientElements=numOfAmbientElements,
            modifiedCCDTime=modifiedCCDTime, modifiedCCDTemp=modifiedCCDTemp, numOfCCDElements=numOfCCDElements, modifiedBatteryTime=modifiedBatteryTime, modifiedBatteryLevel=modifiedBatteryLevel, numOfBatteryElements=numOfBatteryElements, modifiedDiskTime=modifiedDiskTime, 
            modifiedDiskSpace=modifiedDiskSpace, numOfDiskElements=numOfDiskElements)


    return render_template('users/Private/site_desc.html', site=site, mysite = session['current_site'], siteName=siteName, modifiedAmbientTime=modifiedAmbientTime, modifiedAmbientTemp=modifiedAmbientTemp, numOfAmbientElements=numOfAmbientElements,
        modifiedCCDTime=modifiedCCDTime, modifiedCCDTemp=modifiedCCDTemp, numOfCCDElements=numOfCCDElements, modifiedBatteryTime=modifiedBatteryTime, modifiedBatteryLevel=modifiedBatteryLevel, numOfBatteryElements=numOfBatteryElements, modifiedDiskTime=modifiedDiskTime, 
        modifiedDiskSpace=modifiedDiskSpace, numOfDiskElements=numOfDiskElements)


#############################################################
#
#                   Modifying Sites
#
#############################################################

#Add a Site
@mod.route('/sites/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def sites():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

    if "addsite" in request.form:
        return redirect(url_for("users.addsite"))

    if "modifysite" in request.form:
        session['form'] = request.form
        return redirect(url_for("users.modifysite"))

    if "removesite" in request.form:
        session['form'] = request.form
        return redirect(url_for("users.removesite"))


    return render_template('users/Private/sites.html', site=site)

#Add a Site
@mod.route('/addsite/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def addsite():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

    addSiteSuccess = 0          # creating a variable for success of addition of a new site
    verifyLocation = 0          # creating a variable for latitude and longitude verification

    if "verify" in request.form:
        siteName = request.form['siteName']           # extracts the site name from the form
        siteLatitude = request.form['latitude']       # extracts the latitude from the form
        siteLongitude = request.form['longitude']     # extracts the longitude from the form
        siteElevation = request.form['elevation']     # extracts the elevation from the form
        verifyLocation = 1
    
        return render_template('users/Private/addsite.html', site=site, verifyLocation=verifyLocation, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)

    if "submit" in request.form:
        siteName = request.form['siteName']           # extracts the site name from the form
        siteLatitude = request.form['latitude']       # extracts the latitude from the form
        siteLongitude = request.form['longitude']     # extracts the longitude from the form
        siteElevation = request.form['elevation']     # extracts the elevation from the form

        f = open('app/static/Sites/addsite.csv', 'ab')
        g= csv.DictWriter(f, ['ID','NAME','LAT','LON', 'ELEVATION','DELETED'])
        
        try:
            # if the csv file is not empty it will write the new site and add 1 to the ID number
            g.writerow({'ID': str(int(site[len(site)-1]['id']) +1), 'NAME':request.form['siteName'], 'LAT':request.form['latitude'], 'LON':request.form['longitude'],'ELEVATION': request.form['elevation'], "DELETED":0 })
        except:
            # if the csv file is empty then it will assign 1 to the ID
            g.writerow({'ID':1, 'NAME':request.form['siteName'], 'LAT':request.form['latitude'], 'LON':request.form['longitude'], 'ELEVATION': request.form['elevation'], "DELETED":0})
        f.close()    

        #create directory from newly added site
        try:
            os.mkdir(SITES_DIR + siteName)                #creates site directory
            os.mkdir(SITES_DIR + siteName + CAL_DIR)      #creates calibration directory
            os.mkdir(SITES_DIR + siteName + SYS_DIR)      #creates system directory
            os.mkdir(SITES_DIR + siteName + SCHEDULE_DIR) #creates schedule directory
            os.mkdir(SITES_DIR + siteName + IMG_DIR)      #creates schedule directory

            addSiteSuccess = 1
        except:
            flash("Could not make directory.")
        
        return render_template('users/Private/addsite.html',site=site, addSiteSuccess=addSiteSuccess, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)

    return render_template('users/Private/addsite.html',site=site, addSiteSuccess=addSiteSuccess)

#Modify Site
@mod.route('/modifysite/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def modifysite():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

    # keeps the information from the calibration main page
    form = session['form']                   # extracts the site information from the addsite.csv file
    for sites in site:
        if int(sites['id']) == int(form['siteID']):
            session['mySite'] = sites

    siteName = session['mySite']['name']        
    siteLatitude = session['mySite']['lat']
    siteLongitude = session['mySite']['lon']
    siteElevation = session['mySite']['elevation']
    siteID = session['mySite']['id']  
    modifySuccess = 0
    # Creates empty arrays    
    imageName = []                                     # array to store the image name from the image information form
    date = []                                          # array to store thed date from the image information form
    utcTime = []                                       # array to store the UTC Time from when the image was taken
    adminName = []                                     # array to store the administrator who is modifying the form

   

    verifyLocation = 0

    if "verify" in request.form:
        siteName = request.form['siteName']           # extracts the site name from the form
        siteLatitude = request.form['latitude']       # extracts the latitude from the form
        siteLongitude = request.form['longitude']     # extracts the longitude from the form
        siteElevation = request.form['elevation']     # extracts the elevation from the form
        verifyLocation = 1
    
        return render_template('users/Private/modifysite.html', modifySuccess=modifySuccess, site=site, verifyLocation=verifyLocation, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)

    if "submit" in request.form:
        siteName = request.form['siteName']           # extracts the site name from the form
        siteLatitude = float(request.form['latitude'])       # extracts the latitude from the form
        siteLongitude = float(request.form['longitude'])     # extracts the longitude from the form
        siteElevation = float(request.form['elevation'])     # extracts the elevation from the form

        f = open('app/static/Sites/addsite.csv', 'wb')
        g = csv.writer(f)
        g.writerow(['ID','NAME','LAT','LON', 'ELEVATION', 'DELETED'])
        f.close()
        
        f = open('app/static/Sites/addsite.csv', 'ab')
        g= csv.DictWriter(f, ['ID','NAME','LAT','LON','ELEVATION', 'DELETED'])
        for sites in site:
            if siteID == sites['id']:
                sites['name'] = siteName
                sites['lat']  = siteLatitude          # extracts the latitude from the form
                sites['lon']  = siteLongitude         # extracts the longitude from the form
                sites['elevation']  = siteElevation   # extracts the elevation from the form
            g.writerow({'ID':sites['id'], 'NAME':sites['name'], 'LAT':sites['lat'], 'LON':sites['lon'],'ELEVATION':sites['elevation'], 'DELETED':sites['deleted']})
        f.close()

        existingName = SITES_DIR + session['mySite']['name']
        flash(existingName)
        newName = SITES_DIR + siteName
        flash(newName)
        
                
        #os.rename(existingName, newName)
     
        modifySuccess = 1
        return render_template('users/Private/modifysite.html',site=site, modifySuccess=modifySuccess, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)


    return render_template('users/Private/modifysite.html',site=site, modifySuccess=modifySuccess, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)


#Remove a Site
@mod.route('/removesite/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def removesite():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page



    form = session['form']                   # extracts the site information from the addsite.csv file
    
    for sites in site:
        if int(sites['id']) == int(form['siteID']):
            session['mySite'] = sites
    removeSuccess = 0           # initializing variable to store if site has been sucessfully removed

    siteName = session['mySite']['name']        
    siteLatitude = session['mySite']['lat']
    siteLongitude = session['mySite']['lon']
    siteElevation = session['mySite']['elevation']
    siteID = session['mySite']['id'] 

    if "submit" in request.form:
        f = open('app/static/Sites/addsite.csv', 'wb')
        g = csv.writer(f)
        g.writerow(['ID','NAME','LAT','LON', 'ELEVATION', 'DELETED'])
        f.close()
        f = open('app/static/Sites/addsite.csv', 'ab')
        g= csv.DictWriter(f, ['ID','NAME','LAT','LON','ELEVATION', 'DELETED'])
        
        for sites in site:
            if siteID == sites['id']:
                sites['deleted'] = 1
            g.writerow({'ID':sites['id'], 'NAME':sites['name'], 'LAT':sites['lat'], 'LON':sites['lon'],'ELEVATION':sites['elevation'], 'DELETED':sites['deleted']})
        f.close()
        removeSuccess = 1
              
    return render_template('users/Private/removesite.html',site=site, removeSuccess=removeSuccess, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)

#############################################################
#
#             Modifying Calibration for Sites
#
#############################################################

#Calibration Main Page
@mod.route('/calibration/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def calibration():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page


    if "addcalibration" in request.form:
        session['form'] = request.form
        return redirect(url_for("users.addcalibration"))

    if "modifycalibration" in request.form:
        session['form'] = request.form
        return redirect(url_for("users.modifycalibration"))

    
    return render_template('users/Private/calibration.html', site=site)

#Add New Calibration Information Page
@mod.route('/addcalibration/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def addcalibration():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

    # keeps the information from the calibration main page
    form = session['form']
    site = sitesInNetwork()                    # extracts the site information from the addsite.csv file
    for sites in site:
        if int(sites['id']) == int(form['siteID']):
            session['mySite'] = sites

    # Initializing variables
    num_stars = 12
    verify_success = 0
    siteName = session['mySite']['name']        
    siteLatitude = session['mySite']['lat']
    siteLongitude = session['mySite']['lon']
    siteElevation = session['mySite']['elevation']

    # Initializing empty arrays    
    stars = []                                   # array to store all stars
    starName = []                                # array to store all star names of all stars
    azimuth = []                                 # array to store all azimuths of all stars
    elevation = []                               # array to store all elevations of all stars
    icoord = []                                  # array to store all i Coordinates of each star 
    jcoord = []                                  # array to store all i Coordinates of each star

    figureAddress =[]
    
    starName = request.form.getlist('starName')                 # extracts the star name from the form
    azimuth = map(float,request.form.getlist('azimuth'))        # extracts the star azimuth from the form
    elevation = map(float,request.form.getlist('elevation'))    # extracts the star elevation from the form
    iCoord = map(int,request.form.getlist('iCoord'))            # extracts the star i Coordinate from the form
    jCoord = map(int,request.form.getlist('jCoord'))            # extracts the star j Coordinate from the form

    # If image information is submitted
    if "calImageInfo" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        try:
            siteName = session['mySite']['name']            # extracts the site name from the session
            siteLatitude = session['mySite']['lat']         # extracts the latitude from the session
            siteLongitude = session['mySite']['lon']        # extracts the longitude from the session
            siteElevation = session['mySite']['elevation']
            imageName =  request.form['imageName']          # extracts the image name from the image information form
            date = request.form['imageDate']                # extracts the image name from the image information form
            utcTime = request.form['utcTime']
            adminName = request.form['adminName']


            # To be able to write information into the CSV in columns instead of rows
            
            # column 1 includes all the headings
            column1 = ['Site Name', 'Site Latitude', 'Site Longitude', 'Site Elevation', 'Date', 'UTC Time', ' Image Name','Administrator Name']        
            # coulmn 2 includes all the information acquired by the user
            column2 = [siteName, siteLatitude, siteLongitude,0, date, utcTime, imageName, adminName]

            # Opening the calibration file that includes all of the image information in write mode
            calibrationInfoFile = open(SITES_DIR + siteName + CAL_DIR + 'CalibrationInfo.csv' , 'wb')
            g = csv.writer(calibrationInfoFile);
            
            # Opening the calibration file that includes all of the image information in apending mode
        
            for i in range(len(column1)):
                g.writerow([column1[i],  column2[i]])

            calibrationInfoFile.close()
            flash("Submisson Sucessful.")
            

        except:
           flash("Submisson FAIL!") 
        
        return render_template('users/Private/addcalibration.html', site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, 
            iCoord=iCoord, jCoord=jCoord, siteName=siteName, imageName=imageName, date=date, utcTime=utcTime, adminName=adminName, verify_success=verify_success, siteLongitude= siteLongitude, siteElevation=siteElevation)

    # If data verification is clicked
    if "verify" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        starName = request.form.getlist('starName')
        azimuth = map(float,request.form.getlist('azimuth'))
        elevation = map(float,request.form.getlist('elevation'))
        iCoord = map(int,request.form.getlist('iCoord'))
        jCoord = map(int,request.form.getlist('jCoord'))

        stars = [starName, azimuth, elevation, iCoord, jCoord]
        zenith = [ elevation[0], iCoord[0], jCoord[0]]
        
        figureAddress, angle, distance, individualErrors, residualSumSquares, distanceCoefficients, zenith = lensfunction(zenith, azimuth, elevation, iCoord, jCoord, siteName)

        flash("Data checking has succeeded.")
        verify_success = 1

        return  render_template('users/Private/addcalibration.html',site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, 
            iCoord=iCoord, jCoord=jCoord, siteName=siteName, verify_success=verify_success,
            figureAddress=figureAddress, angle=angle, distance=distance, individualErrors=individualErrors, residualSumSquares=residualSumSquares, siteLatitude=siteLatitude,
            siteLongitude= siteLongitude, siteElevation=siteElevation)
        

    if "submit" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        try:
            # Extracting site information
            siteName = session['mySite']['name']                        # extracts the site name from the session
            siteLatitude = session['mySite']['lat']                     # extracts the site latitude from the session
            siteLongitude = session['mySite']['lon']                    # extracts the site longitude from the session

            # Star Information from the form
            starName = request.form.getlist('starName')                 # extracts the star name from the form
            azimuth = map(float,request.form.getlist('azimuth'))        # extracts the star azimuth from the form
            elevation = map(float,request.form.getlist('elevation'))    # extracts the star elevation from the form
            iCoord = map(int,request.form.getlist('iCoord'))            # extracts the star i Coordinate from the form
            jCoord = map(int,request.form.getlist('jCoord'))            # extracts the star j Coordinate from the form
            
            zenith = [elevation[0], iCoord[0], jCoord[0]]
            
            figureAddress, angle, distance, individualErrors, residualSumSquares, distanceCoefficients, zenith = lensfunction(zenith, azimuth, elevation, iCoord, jCoord, siteName)
            #pixelsToLatLon(zenith, distanceCoefficients, siteLatitude, siteLongitude)
            generator = pixelsToLatLon(zenith, distanceCoefficients, siteLatitude, siteLongitude, siteName)
            
            # Writing the calibration star files 
            calibrationStarFile = open(SITES_DIR + siteName + CAL_DIR + 'Calibration.csv' , 'wb')
            g = csv.writer(calibrationStarFile);
            g.writerow(['Star Name','Azimuth','Elevation','i Coordinate', 'j Coordinate'])
            calibrationStarFile.close()

            calibrationStarFile = open(SITES_DIR + siteName + CAL_DIR + 'Calibration.csv' , 'ab')
            g = csv.DictWriter(calibrationStarFile,['Star Name','Azimuth','Elevation','i Coordinate', 'j Coordinate'])

            for i in range(len(elevation)):
                g.writerow({'Star Name': starName[i], 'Azimuth':azimuth[i], 'Elevation': elevation[i], 'i Coordinate':iCoord[i], 'j Coordinate': jCoord[i]})

            calibrationStarFile.close()
            flash("Star data submisson sucessful.")

        except:
           flash("Submisson FAIL")

        return render_template('users/Private/addcalibration.html', site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, 
            iCoord=iCoord, jCoord=jCoord, siteName=siteName, imageName=imageName, date=date, utcTime=utcTime, adminName=adminName, verify_success=verify_success,
            figureAddress=figureAddress, angle=angle, distance=distance, individualErrors=individualErrors, residualSumSquares=residualSumSquares, siteLatitude=siteLatitude,
            siteLongitude= siteLongitude, siteElevation=siteElevation)

    return render_template('users/Private/addcalibration.html',site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, iCoord=iCoord, jCoord=jCoord, 
        verify_success=verify_success,figureAddress=figureAddress, siteLatitude=siteLatitude, siteLongitude= siteLongitude, siteName= siteName, siteElevation=siteElevation)

#Add modify calibration information 
@mod.route('/modifycalibration/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def modifycalibration():
    form = session['form']
    
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

    # Current Site in Session
    for sites in site:
        if int(sites['id']) == int(form['siteID']):
            session['mySite'] = sites

    siteName = session['mySite']['name']
    siteLatitude = session['mySite']['lat']
    siteLongitude = session['mySite']['lon']
    siteElevation = session['mySite']['elevation']
    num_stars = 12
    verify_success = 0
    image_success = 0
    submit_success = 0
   
    # Calibration image information extraction
    calibrationInfoFile = open(SITES_DIR + siteName + CAL_DIR + 'CalibrationInfo.csv' , 'r')
    text1 = calibrationInfoFile.read()
    p = StringIO(text1)
    data1 = genfromtxt(p, dtype=None, delimiter = ',')    # extracts the data from CSV file 
    calibrationInfoFile.close()

    # Creates empty arrays    
    imageName = []                                     # array to store the image name from the image information form
    date = []                                          # array to store thed date from the image information form
    utcTime = []                                       # array to store the UTC Time from when the image was taken
    adminName = []                                     # array to store the administrator who is modifying the form

    # Assigning variables
    imageName = ''.join(map(str, data1[6][1]))         # extracts the image name from data1  
    date = ''.join(map(str, data1[4][1]))              # extracts the date from data1 
    utcTime = ''.join(map(str, data1[5][1]))
    adminName = ''.join(map(str, data1[7][1]))


    # Calibration star information extraction
    calibrationStarFile = open(SITES_DIR + siteName + CAL_DIR + 'Calibration.csv' , 'r')
    text = calibrationStarFile.read()
    s = StringIO(text)
    data = genfromtxt(s, delimiter = ',')
    calibrationStarFile.close()

    # Creates empty arrays    
    stars = []                                   # array to store all stars
    starName = []                                # array to store all star names of all stars
    azimuth = []                                 # array to store all azimuths of all stars
    elevation = []                               # array to store all elevations of all stars
    icoord = []                                  # array to store all i Coordinates of each star 
    jcoord = []                                  # array to store all i Coordinates of each star
    figureAddress =[]
    angle = []                                   # array to store angle from zenith obtained from lense function
    distance = []
    individualErrors = [] 
    residualSumSquares = 0
    distanceCoefficients = []
    latitudes = []
    longitudes = []
    
    lengthData = data.shape[0]                  # gives the size of the array
    

    azimuth = array((data[1:lengthData,1]))     # extracts the azimuth from data
    elevation = array((data[1:lengthData,2]))   # extracts the elevation from data
    iCoord = array((data[1:lengthData,3]))      # extracts the i Coordinate from data
    jCoord = array((data[1:lengthData,4]))      # extracts the j Coordinate from data
    iCoord = map(int, iCoord)
    jCoord = map(int, jCoord)

    # To be able to obtain the star names has strings
    s = StringIO(text)
    data = genfromtxt(s, dtype=None, delimiter = ',')
    for i in range(1,lengthData):
        starInfo = data[i]
        starName = append(starName, starInfo[0]) 


    # For data submission for image information calibration when "calImageInfo" clicked
    if "calImageInfo" in request.form:
        site = sitesInNetwork()
    
        # For Navigation Bar
        for sites in site:
            if 'mysite'+sites['id'] in request.form:
                session['current_id'] = sites['id']
                session['current_site'] = sites['name']
                return redirect(url_for('users.site_desc'))

        try:
            siteName = session['mySite']['name']
            siteLatitude = session['mySite']['lat']
            siteLongitude = session['mySite']['lon']
            siteElevation = session['mySite']['elevation']
            imageName =  request.form['imageName']
            date = request.form['imageDate']
            utcTime = request.form['utcTime']
            adminName = request.form['adminName']

            # To be able to write information into the CSV in columns instead of rows
            
            # column 1 includes all the headings
            column1 = ['Site Name', 'Site Latitude', 'Site Longitude', 'Site Elevation', 'Date', 'UTC Time', ' Image Name','Administrator Name']        
            # coulmn 2 includes all the information acquired by the user
            column2 = [siteName, siteLatitude, siteLongitude,0, date, utcTime, imageName, adminName]

            # Opening the calibration file that includes all of the image information in write mode
            calibrationInfoFile = open(SITES_DIR + siteName + CAL_DIR + 'CalibrationInfo.csv' , 'wb')
            g = csv.writer(calibrationInfoFile);
            
            for i in range(len(column1)):
                g.writerow([column1[i],  column2[i]])

            calibrationInfoFile.close()
            flash("Submisson Sucessful.")
          
        except:
           flash("Submisson FAIL - Maybe file is open?") 

        return render_template('users/Private/modifycalibration.html', site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, 
            iCoord=iCoord, jCoord=jCoord, siteName=siteName, imageName=imageName, date=date, utcTime=utcTime, adminName=adminName, verify_success=verify_success,
            figureAddress=figureAddress, angle=angle, distance=distance, individualErrors=individualErrors, residualSumSquares=residualSumSquares, siteLatitude=siteLatitude,
            siteLongitude= siteLongitude, siteElevation=siteElevation)


    # For data verification when "Verify" button has been clicked
    if "verify" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        starName = request.form.getlist('starName')
        azimuth = map(float,request.form.getlist('azimuth'))
        elevation = map(float,request.form.getlist('elevation'))
        iCoord = map(int,request.form.getlist('iCoord'))
        jCoord = map(int,request.form.getlist('jCoord'))
        
    
        stars = [starName, azimuth, elevation, iCoord, jCoord]
        zenith = [ elevation[0], iCoord[0], jCoord[0]]
        
        figureAddress, angle, distance, individualErrors, residualSumSquares, distanceCoefficients, zenith = lensfunction(zenith, azimuth, elevation, iCoord, jCoord, siteName)
        flash("Data checking has succeeded.")
        verify_success = 1

        return render_template('users/Private/modifycalibration.html', site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, 
        iCoord=iCoord, jCoord=jCoord, siteName=siteName, imageName=imageName, date=date, utcTime=utcTime, adminName=adminName, verify_success=verify_success,
        figureAddress=figureAddress, angle=angle, distance=distance, individualErrors=individualErrors, residualSumSquares=residualSumSquares, siteLatitude=siteLatitude,
        siteLongitude= siteLongitude, siteElevation=siteElevation)

    # For data verification when "Submit" button has been clicked
    if "submit" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
        
        try:
            
            siteName = session['mySite']['name']
            siteLatitude = session['mySite']['lat']
            siteLongitude = session['mySite']['lon']

            # obtaining  necessary fields for csv writting       
            starName = request.form.getlist('starName')
            azimuth = map(float,request.form.getlist('azimuth'))
            elevation = map(float,request.form.getlist('elevation'))
            iCoord = map(int,request.form.getlist('iCoord'))
            jCoord = map(int,request.form.getlist('jCoord'))

            zenith = [elevation[0], iCoord[0], jCoord[0]]
            
            figureAddress, angle, distance, individualErrors, residualSumSquares, distanceCoefficients, zenith = lensfunction(zenith, azimuth, elevation, iCoord, jCoord, siteName)
            #pixelsToLatLon(zenith, distanceCoefficients, siteLatitude, siteLongitude)
            generator = pixelsToLatLon(zenith, distanceCoefficients, siteLatitude, siteLongitude, siteName)
                        
            # Writing the CSV file for star information
            calibrationStarFile = open(SITES_DIR + siteName + CAL_DIR + 'Calibration.csv' , 'wb')
            g = csv.writer(calibrationStarFile);
            g.writerow(['Star Name','Azimuth','Elevation','i Coordinate', 'j Coordinate'])
            calibrationStarFile.close()

            calibrationStarFile = open(SITES_DIR + siteName + CAL_DIR + 'Calibration.csv' , 'ab')
            g = csv.DictWriter(calibrationStarFile,['Star Name','Azimuth','Elevation','i Coordinate', 'j Coordinate'])

            for i in range(len(elevation)):
                g.writerow({'Star Name': starName[i], 'Azimuth':azimuth[i], 'Elevation': elevation[i], 'i Coordinate':iCoord[i], 'j Coordinate': jCoord[i]})

            calibrationStarFile.close()
            flash("Submisson Sucessful.")


        except:
           flash("Submisson FAIL - Maybe file is open?")  

        return render_template('users/Private/modifycalibration.html', site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, iCoord=iCoord, jCoord=jCoord, siteName=siteName, imageName=imageName, date=date, utcTime=utcTime, adminName=adminName, verify_success=verify_success, figureAddress=figureAddress, angle=angle, distance=distance,   
            siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation)

    return render_template('users/Private/modifycalibration.html', site=site, num_stars=num_stars, starName=starName, azimuth=azimuth, elevation=elevation, 
        iCoord=iCoord, jCoord=jCoord, siteName=siteName, imageName=imageName, date=date, utcTime=utcTime, adminName=adminName, verify_success=verify_success,
        figureAddress=figureAddress, angle=angle, distance=distance, individualErrors=individualErrors, residualSumSquares=residualSumSquares, siteLatitude=siteLatitude,
        siteLongitude= siteLongitude, siteElevation=siteElevation, distanceCoefficients=distanceCoefficients)


#############################################################
#
#             Modifying Schedule for Sites
#
#############################################################
#Sheduling Main Page
@mod.route('/schedule/', methods = ['GET', 'POST'])
#@adminrequires_login     # Requires adminitrator to be logged in to access page
def schedule():
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page
    

    if "addschedule" in request.form:
        session['form'] = request.form
        return redirect(url_for("users.addschedule"))

    if "modifyschedule" in request.form:
        session['form'] = request.form
        return redirect(url_for("users.modifyschedule"))
        
    return render_template('users/Private/schedule.html', site=site)

#Add New Schedule Information Main Page
@mod.route('/addschedule/', methods = ['GET', 'POST'])
@adminrequires_login     # Requires adminitrator to be logged in to access page
def addschedule():
    #Site information
    form = session['form']
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

    # Extracting Site Information
    siteName = session['mySite']['name']
    siteLatitude = session['mySite']['lat']
    siteLongitude = session['mySite']['lon']
    siteElevation = session['mySite']['elevation']

    # Form variables initialization
    secondsSize = 60            # Amount of seconds in a minute
    minutesSize = 60            # Amount of minutes in an hour
    hoursSize = 24              # Amount of hours in a day
    sizeBinX = 5
    targetTempSize = 21
    preCoolingSize = array(([0, 15, 30, 45, 60]))
    restPeriodSize = array(([15, 30, 45, 60]))
    addNewSchedule = 0
    adminName = []
    date = []

    # Initializing variables

    if "submit" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        # Obtaining information
        adminName = request.form['adminName']
        date = request.form['date']
        time = request.form['time']

        # Giving nonmeaning variables a value
        startTime = 000000
        stopTime = 235959
        
        # Obtaining the "exposure time" from the form
        exposureMinutes = int(request.form['exposureMinutes'])
        exposureSeconds = int(request.form['exposureSeconds'])
        
        # Renaming variable for better track keeping
        exposureTime = exposureSeconds + exposureMinutes*100

        # Obtaining the "Rest Period Time" from the form
        restPeriodTimeSeconds = int(request.form['restPeriod'])
        # Renaming variable for better track keeping
        restPeriodTime = restPeriodTimeSeconds
        
        # Obtaining the "Binning X" from the form in the X direction
        
        binning = int(request.form['binX'])
        binX = binning
        binY = binX
        
        # Obtaining the "Target Temperature" from the form
        targetTemperature = int(request.form['targetTemp'])
        
        # Obtaining the "Pre-Cooling" from the form
        preCoolingMinutes = int(request.form['preCooling'])

        # Setting Pre-Cooling Duration 
        if preCoolingMinutes == 0:
            preCoolingDuration = preCoolingMinutes
            preCooling = 0                                      # Turns OFF Pre-Cooling

        if preCoolingMinutes > 0: 
            preCoolingDuration = preCoolingMinutes
            preCooling = 1                                      # Turns ON Pre-Cooling

        column1 = ['Site', 'Admin', 'Date', 'Time', 'Start Time', 'Stop Time', 'Exposure Time', 'Rest Period', 'BinX', 'BinY', 'Target Temperature', 'Pre-Cooling', 'Pre-Cooling Duration']
        column2 = [siteName, adminName, date, time, startTime, stopTime, exposureTime, restPeriodTime, binX, binY, targetTemperature, preCooling, preCoolingDuration ]
        
        try:
            # Opening the schedule file that includes all of the informaiton
            scheduleFile = open(SITES_DIR + siteName + SCHEDULE_DIR + 'Schedule.csv' , 'wb')
            g = csv.writer(scheduleFile);

        
            for i in range(len(column1)):
                g.writerow([column1[i],  column2[i]])

            scheduleFile.close()
            addNewSchedule = 1
            
        except:
            flash("Submission Fail")
         
        return render_template('users/Private/addschedule.html', site=site, addNewSchedule=addNewSchedule, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation, secondsSize=secondsSize, hoursSize=hoursSize, 
            minutesSize=minutesSize, sizeBinX=sizeBinX, targetTempSize=targetTempSize, preCoolingSize=preCoolingSize, restPeriodSize=restPeriodSize, exposureMinutes=exposureMinutes, exposureSeconds=exposureSeconds,
            targetTemperature=targetTemperature, preCoolingMinutes=preCoolingMinutes, binning=binning, 
            restPeriodTimeSeconds=restPeriodTimeSeconds, date=date, time=time, adminName=adminName)
    
    return render_template('users/Private/addschedule.html',site=site, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation, secondsSize=secondsSize, hoursSize=hoursSize, 
        minutesSize=minutesSize, sizeBinX=sizeBinX, targetTempSize=targetTempSize, preCoolingSize=preCoolingSize, addNewSchedule=addNewSchedule, restPeriodSize=restPeriodSize )

#Add new calibration information Main Page
@mod.route('/addscheduleconfirmation/', methods = ['GET', 'POST'])

def addscheduleconfirmation():

    return render_template('users/Private/addscheduleconfirmation.html')


#Add new calibration information Main Page
@mod.route('/modifyschedule/', methods = ['GET', 'POST'])
@adminrequires_login     # Requires adminitrator to be logged in to access page
def modifyschedule():
    #Site information
    form = session['form']
    
    # For the navigation bar in all the pages
    site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
    for sites in site:
        if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
            session['current_id'] = sites['id']         # the site id from the current session
            session['current_site'] = sites['name']     # the site name from the current session
            return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

    # Site information
    site = sitesInNetwork()
    for sites in site:
        if int(sites['id']) == int(form['siteID']):
            session['mySite'] = sites
    
    # Extracting Site Information
    siteName = session['mySite']['name']            # extracts the site name from the session
    siteLatitude = session['mySite']['lat']         # extracts the latitude from the session
    siteLongitude = session['mySite']['lon']        # extracts the longitude from the session
    siteElevation = session['mySite']['elevation']  # extracts the elevation from the session
    
    # Form variables initialization
    secondsSize = 60
    minutesSize = 60
    hoursSize = 24
    sizeBinX = 5
    targetTempSize = 21
    preCoolingSize = array(([0, 15, 30, 45, 60]))
    restPeriodSize = array(([15, 30, 45, 60]))
    modifySuccess = 0
    adminName = []
    date = []
    scheduleData = []

    scheduleData = scheduleFileDataExtraction(siteName)                 # Funtion created to extract the information from the CSV file needed for the schedule 
    adminName = str(scheduleData['adminName'])                             # extracts the administrator's name from the schedule
    date = str(scheduleData['date'])                                         # extracts the date from the schedule
    time = str(scheduleData['time'])                                         # extracts the time from the schedule
    exposureSeconds = int(scheduleData['exposureSeconds'])
    exposureMinutes = int(scheduleData['exposureMinutes'])
    restPeriodTimeSeconds = scheduleData['restPeriodTimeSeconds']
    binX = scheduleData['binX']
    binY= scheduleData['binY']
    targetTemperature = scheduleData['targetTemperature']
    preCoolingMinutes = scheduleData['preCoolingMinutes']
    binning=binX
        
    # If submission is entered: 
    if "submit" in request.form:
        # For the navigation bar in all the pages
        site = sitesInNetwork()                             # calls the function to find the sites in the network using the addsite.csv file
        for sites in site:
            if 'mysite'+sites['id'] in request.form:        # iterates through all the sites  
                session['current_id'] = sites['id']         # the site id from the current session
                session['current_site'] = sites['name']     # the site name from the current session
                return redirect(url_for('users.site_desc')) # redirects the user to the site chosen from the navigation bar to the "Site System Monitoring" page

        siteName = session['mySite']['name'] 
        
        # Obtaining information
        adminName = request.form['adminName']
        date = request.form['date']
        time = request.form['time']

        # Giving nonmeaning variables a value
        startTime = 000000
        stopTime = 235959
        
        # Obtaining the "exposure time" from the form
        exposureTimeMinutes = int(request.form['exposureMinutes'])
        exposureTimeSeconds = int(request.form['exposureSeconds'])
        
        exposureTime = exposureTimeSeconds + exposureTimeMinutes*100

        # Obtaining the "Rest Period Time" from the form
        restPeriodTimeSeconds = int(request.form['restPeriod'])
        # Renaming variable for better track keeping
        restPeriodTime = restPeriodTimeSeconds
        
        # Obtaining the "Binning X" from the form in the X direction
        binX = int(request.form['binX'])
        binY = binX
        
        # Obtaining the "Target Temperature" from the form
        targetTemperature = int(request.form['targetTemp'])
        
        # Obtaining the "Pre-Cooling" from the form
        preCoolingMinutes = int(request.form['preCooling'])

        # Setting Pre-Cooling Duration 
        if preCoolingMinutes == 0:
            preCoolingDuration = preCoolingMinutes
            preCooling = 0                                      # Turns OFF Pre-Cooling

        if preCoolingMinutes > 0: 
            preCoolingDuration = preCoolingMinutes
            preCooling = 1                                      # Turns ON Pre-Cooling

        column1 = ['Site', 'Admin', 'Date', 'Time', 'Start Time', 'Stop Time', 'Exposure Time', 'Rest Period', 'BinX', 'BinY', 'Target Temperature', 'Pre-Cooling', 'Pre-Cooling Duration']
        column2 = [siteName, adminName, date, time, startTime, stopTime, exposureTime, restPeriodTime, binX, binY, targetTemperature, preCooling, preCoolingDuration ]
        
        try:
            # Opening the schedule file that includes all of the informaiton
            scheduleFile = open(SITES_DIR + siteName + SCHEDULE_DIR + 'Schedule.csv' , 'wb')
            g = csv.writer(scheduleFile);

        
            for i in range(len(column1)):
                g.writerow([column1[i],  column2[i]])

            scheduleFile.close()
            modifySuccess = 1
        except:
            flash("Submission Fail")
        
        scheduleData = scheduleFileDataExtraction(siteName)                 # Funtion created to extract the information from the CSV file needed for the schedule 
        adminName = str(scheduleData['adminName'])                          # extracts the administrator's name from the schedule
        date = str(scheduleData['date'])                                    # extracts the date from the schedule
        time = str(scheduleData['time'])                                    # extracts the time from the schedule
        exposureSeconds = int(scheduleData['exposureSeconds'])
        exposureMinutes = int(scheduleData['exposureMinutes'])
        restPeriodTimeSeconds = scheduleData['restPeriodTimeSeconds']
        binX = scheduleData['binX']
        binY= scheduleData['binY']
        targetTemperature = scheduleData['targetTemperature']
        preCoolingMinutes = scheduleData['preCoolingMinutes']
        binning=binX


        return render_template('users/Private/modifyschedule.html', site=site, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation, secondsSize=secondsSize, hoursSize=hoursSize, 
            minutesSize=minutesSize, sizeBinX=sizeBinX, targetTempSize=targetTempSize, preCoolingSize=preCoolingSize, modifySuccess=modifySuccess, restPeriodSize=restPeriodSize, exposureMinutes=exposureMinutes, exposureSeconds=exposureSeconds,
            targetTemperature=targetTemperature, preCoolingMinutes=preCoolingMinutes, binning=binning, 
            restPeriodTimeSeconds=restPeriodTimeSeconds, date=date, time=time, adminName=adminName)

    return render_template('users/Private/modifyschedule.html',  site=site, siteName=siteName, siteLatitude=siteLatitude, siteLongitude=siteLongitude, siteElevation=siteElevation, secondsSize=secondsSize, hoursSize=hoursSize, 
        minutesSize=minutesSize, sizeBinX=sizeBinX, targetTempSize=targetTempSize, preCoolingSize=preCoolingSize, modifySuccess=modifySuccess, restPeriodSize=restPeriodSize, exposureMinutes=exposureMinutes, exposureSeconds=exposureSeconds,
        targetTemperature=targetTemperature, preCoolingMinutes=preCoolingMinutes, binning=binning, 
        restPeriodTimeSeconds=restPeriodTimeSeconds, date=date, time=time, adminName=adminName)
