# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 16:41:15 2018

Script to read 3-column csv files of lat, lon, mmi and convert to shakemap-
compliant *_dat.xml file.

Usage:
    python mmicsv2xml.py <path/to/param/file>
    
    e.g.
        python mmicsv2xml.py 201605201814/raw/peterman.param
        
    Format of the param file, e.g.:
        
        locstring      = Petermann Ranges, NT
        eqla           = -25.579
        eqlo           = 129.832
        eqdep          = 0.0
        eqmw           = 6.1
        yyyymmddHHMMSS = 20160520181402
        csvfile        = Pettermann_ranges_2016_intensities.csv   

Output xml file is written to event's "current" folder, e.g.

    <INSTALL_DIR>/data/201605201814/current/20160520181402_dat.xml

@author: u56903
"""

from roman import toRoman
from mapping_tools import distance
from numpy import array, arange, mean, where
from mmi_tools import mmi2pgm_worden12, cmpsps2g
from sys import argv
from os import path, sep, mkdir
import time

# set param file:
#paramFile = argv[1]
#inFolder = path.split(paramFile)[0]

# load param file
#lines = open(paramFile).readlines()

# set params
locstring      = 'NT-SA-WA Border Region'
eqla           = -26.0
eqlo           = 129.0
eqdp          = 5.0
eqmw           = 5.
yyyymmddHHMMSS = '25991212010101'
csvfile        = 'test_data.csv'
 
# set output folder
outFolder = ''

'''
csvfile format with one header line:
    LAT, LON, MMI
'''

# get station data
lines = open(csvfile).readlines()
lat = []
lon = []
mmi = []
for line in lines[1:]:
    dat = line.strip().split(',')
    if lines[0].startswith('ORIGINTIME,'):
        lat.append(float(dat[6])) # keep as string
        lon.append(float(dat[5])) # keep as string
        mmi.append(float(dat[7]))
    else:
        try:
            lat.append(float(dat[0])) # keep as string
            lon.append(float(dat[1])) # keep as string
            mmi.append(float(dat[2])) 
        except:
            print('Missing data!')
    
# grid data and take average
mmi = array(mmi)
lat = array(lat)
lon = array(lon)

# set grid size (in degrees)
dd = 0.05
d2 = dd / 2.

# set lat/lon ranges +/- 6 degrees
lonrng = arange(eqlo-6., eqlo+6., dd)
latrng = arange(eqla-6., eqla+6., dd)

# set min number of obs/per cell required for average mmi - depends on total obs
minObs = 0
if len(mmi) > 100:
    minObs = 1
elif minObs > 500:
    minObs = 2

# loop thru grid cells values and take average
avmmi = []
avlo  = []
avla  = []
for lo in lonrng:
     for la in latrng:
         idx = where((lon >= lo-d2) & (lon < lo+d2) & \
                     (lat >= la-d2) & (lat < la+d2))[0]
         
         # add to list greater than minObs
         if len(idx) > minObs:
             avmmi.append(mean(mmi[idx]))
             
             # use data lat/lon
             if len(idx) == 1:
                 avlo.append(lon[idx])
                 avla.append(lat[idx])
             else:
                 avlo.append(lo)
                 avla.append(la)

# make earthquake header
smtxt = '<shakemap-data code_version="4.0 GSM" map_version="1">\n'
earthquake = ''.join(('<earthquake id="', str(yyyymmddHHMMSS), '" lat="', str(eqla), '" lon="', \
                       str(eqlo), '" mag="', str(eqmw), '" year="', str(yyyymmddHHMMSS)[0:4], '" month="', str(yyyymmddHHMMSS)[4:6], \
                       '" day="', str(yyyymmddHHMMSS)[6:8], '" hour="', str(yyyymmddHHMMSS)[8:10], '" minute="', str(yyyymmddHHMMSS)[10:12], \
                       '" second="', str(yyyymmddHHMMSS)[12:], '" timezone="GMT" depth="', str(eqdp), \
                       '" locstring="', locstring, '" created="', str(int(time.time())), '"/>\n'))

endtxt = '</stationlist>\n</shakemap-data>'

smtxt += earthquake + '<stationlist created="' + str(int(time.time())) + '">\n'

# start alternate raw_mmi file
raw_mmi = ','.join((yyyymmddHHMMSS[:-2], str(eqlo), str(eqla), str(eqdp), str(eqmw), locstring.replace(',', ''))) + '\n'

# now loop thru station data and get GM values
i = 0
# get stn loc
for la, lo, mi in zip(avla, avlo, avmmi):
    if mi > 0:
        code = 'OBS_'+str(i+1)
        source = 'AU'
        netid  = 'Intensity'
        name = 'Unknown (Intensity ' + toRoman(round(mi)) +')'
        
        station = '"'.join(('<station code=', code, ' name=', name, ' insttype="Observed', \
                            ' lat=', str('%0.4f' % float(la)), ' lon=', str('%0.4f' % float(lo)), ' source=', source, \
                            ' netid=', netid, ' commtype="Intensity" intensity=', str('%0.1f' % mi),'>\n'))
    
        smtxt += station
        
        raw_mmi += ','.join((str('%0.4f' % float(lo)), str('%0.4f' % float(la)), str('%0.1f' % mi))) + '\n'
    
        # get distance
        rngkm = distance(float(la), float(lo), float(eqla), float(eqlo))[0]       
        
        # make mmi comp
        pga = cmpsps2g(mmi2pgm_worden12([mi], 'pga', float(eqmw), rngkm)[0])
        pgv = (mmi2pgm_worden12([mi], 'pgv', eqmw, rngkm)[0])
        
        acc = str('%0.4f' % pga)
        vel = str('%0.4f' % pgv)
        smtxt += '"'.join(('    <comp name="DERIVED">\n        <acc value=', acc, '/>\n        <vel value=', \
                          vel, '/>\n    </comp>\n</station>\n'))
           
        i += 1

# end text
smtxt += endtxt

# check to see if exists
currentFolder = path.join(outFolder, 'current')
if path.isdir(currentFolder) == False:
    mkdir(currentFolder)

# write to xml file
f = open(path.join(outFolder, 'current', str(yyyymmddHHMMSS)+'_obs_dat.xml'), 'w')
f.write(smtxt)
f.close()

###############################################################################
# export simplified file

# write to xml file
f = open(yyyymmddHHMMSS+'_raw_mmi.csv', 'w')
f.write(raw_mmi)
f.close()
