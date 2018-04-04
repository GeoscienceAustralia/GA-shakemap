# -*- coding: utf-8 -*-
"""
Created on Thu Mar 01 16:41:15 2018

Script to read 3-column csv files of lat, lon, mmi and convert to shakemap-
compliant *_dat.xml file

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

@author: u56903
"""

from roman import toRoman
from mapping_tools import distance
from numpy import array, arange, mean, where
from mmi_tools import mmi2pgm_worden12, cmpsps2g
from sys import argv
from os import path, sep
import time

# set param file:
paramFile = argv[1]
inFolder = path.split(paramFile)[0]

# load param file
lines = open(paramFile).readlines()

# set params
locstring  = lines[0].split('=')[-1].strip()
eqla    = float(lines[1].split('=')[-1].strip())
eqlo    = float(lines[2].split('=')[-1].strip())
eqdp    = float(lines[3].split('=')[-1].strip())
eqmw    = float(lines[4].split('=')[-1].strip())
yyyymmddHHMMSS = lines[5].split('=')[-1].strip()
csvfile = lines[6].split('=')[-1].strip()

# set output folder
outFolder = sep.join(paramFile.split(sep)[:3])

'''
csvfile format with one header line:
    LAT, LON, MMI
'''

# get station data
lines = open(path.join(inFolder, csvfile)).readlines()[1:]
lat = []
lon = []
mmi = []
for line in lines:
    dat = line.strip().split(',')
    lat.append(float(dat[0])) # keep as string
    lon.append(float(dat[1])) # keep as string
    mmi.append(float(dat[2]))
    
# grid data and take average
mmi = array(mmi)
lat = array(lat)
lon = array(lon)

dd = 0.05
d2 = dd / 2.

lonrng = arange(eqlo-6., eqlo+6., dd)
latrng = arange(eqla-6., eqla+6., dd)

# set min number of obs required for average mmi
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

# now loop thru station data and get GM values
i = 0
# get stn loc
for la, lo, mi in zip(avla, avlo, avmmi):
    code = 'OBS_'+str(i+1)
    source = 'AU'
    netid  = 'Intensity'
    name = 'Unknown (Intensity ' + toRoman(round(mi)) +')'
    
    station = '"'.join(('<station code=', code, ' name=', name, ' insttype="Observed', \
                        ' lat=', str('%0.4f' % float(la)), ' lon=', str('%0.4f' % float(lo)), ' source=', source, \
                        ' netid=', netid, ' commtype="Intensity" intensity=', str('%0.1f' % mi),'>\n'))

    smtxt += station

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

# write to text
f = open(path.join(outFolder, 'current', str(yyyymmddHHMMSS)+'_dat.xml'), 'wb')
f.write(smtxt)
f.close()