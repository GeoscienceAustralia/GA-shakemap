basfrom matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, where
#from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from mapping_tools import get_map_polygons, mask_outside_polygons
from sys import argv
from os import getcwd

import matplotlib as mpl
mpl.style.use('classic')

plt.rcParams['pdf.fonttype'] = 42

'''
 run map_lm_intensities.py Lake_Muir_16_09_geolocated_MMI.csv '2018-09-16'
'''
##########################################################################################
# settings
##########################################################################################

csvmmi = argv[1]
#evid = argv[2] # fmt= YYYY-MM-DD

'''
if evid.startswith('2018-09'):
    eqlat = -34.390
    eqlon = 116.799
    mlstr = '5.7'
elif evid.startswith('2018-11'):
    eqlat = -34.423
    eqlon = 116.787
    mlstr = '5.4'
'''
'''
# Carnarvon
eqla = -23.304
eqlo = 112.623
degrng = 7.
latoff = -3.
lonoff = 5.
mstr = '5.9'
place = 'Offshore Carnarvon, WA'
evid = '2018-12-16'

'''
# Banda Sea
eqla = -7.577
eqlo = 128.762
degrng = 5.
latoff = -4.
lonoff = 3.5
mstr = '6.5'
place = 'Banda Sea'
evid = '2018-12-01'


# set grid size (in degrees)
if degrng <= 3:
    dd = 0.08
    d2 = dd / 2.
    d3 = d2 * 0.8
    latrng = arange(eqla-2*degrng, eqla+2*degrng, dd*0.8) # kluge to get approx square
    pltbuffer = 0.1
    txtoff = 0.025
elif degrng > 3 and degrng < 6:
    dd = 0.15
    d2 = dd / 2.
    d3 = d2 * 0.85
    latrng = arange(eqla-2*degrng, eqla+2*degrng, dd*0.85) # kluge to get approx square
    pltbuffer = 0.5
    txtoff = 0.07
elif degrng >= 6:
    dd = 0.2
    d2 = dd / 2.
    d3 = d2 * 0.85
    latrng = arange(eqla-2*degrng, eqla+2*degrng, dd*0.85) # kluge to get approx square
    pltbuffer = 0.7
    txtoff = 0.1

# set lat/lon ranges +/- 6 degrees
lonrng = arange(eqlo-2*degrng, eqlo+2*degrng, dd)

minObs = 1 # just for MMI plotting

##########################################################################################
# parse mmi csv
##########################################################################################

lines = open(csvmmi).readlines()
mmi = []
lat = []
lon = []
for line in lines:
    dat = line.strip().split(',')
    mmi.append(float(dat[2]))
    lon.append(float(dat[1]))
    lat.append(float(dat[0]))
    

mmi = array(mmi)
lon = array(lon)
lat = array(lat)

##########################################################################################
# set up map
##########################################################################################

# set map on fly
if degrng > 2 and degrng < 5:
    parSpace = 1.
    merSpace = 2.
elif degrng >= 5.:
    parSpace = 2.
    merSpace = 2.

# make bounds on the fly
urcrnrlat = eqla + degrng + latoff
llcrnrlat = eqla - degrng + latoff
urcrnrlon = eqlo + degrng + lonoff
llcrnrlon = eqlo - degrng + lonoff

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=100)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()

m.fillcontinents(color='0.9',lake_color='lightskyblue')
m.drawparallels(arange(-90.,90.,parSpace), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
m.drawmeridians(arange(0.,360.,merSpace), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
#m.drawmapscale(145.4, -39.2, 145.4, -39.2, 100, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# read topo
##########################################################################################
#mdiv = 500.

'''
print 'Reading topo file...'
netcdffile = '//Users//tallen//Documents//DATA//SRTM03//victoria_srtm03.grd'
nc = NetCDFFile(netcdffile)

zscale =2. #gray
zscale = 75. #colour
#zscale = 10. #colour
data = nc.variables['z'][:] / zscale
lons = nc.variables['x'][:]
lats = nc.variables['y'][:]

# transform to metres
mdiv = 500.
nx = int((m.xmax-m.xmin)/mdiv)+1
ny = int((m.ymax-m.ymin)/mdiv)+1
topodat = m.transform_scalar(data,lons,lats,nx,ny)

##########################################################################################
# plot intensity grid from shakemap
##########################################################################################

print 'Reading MMI file...'
nc = NetCDFFile('mi.grd')

#zscale =1. #colour
data = nc.variables['z'][:]
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/mdiv)+1
ny = int((m.ymax-m.ymin)/mdiv)+1

mmidat = m.transform_scalar(data,lons,lats,nx,ny)

print 'Getting colormap...'
# get colormap
cptfile = 'mi_pop.cpt'
cmap, zvals = cpt2colormap(cptfile, 256)
#cmap=plt.get_cmap('Spectral',256)

# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 0)
norm = mpl.colors.Normalize(vmin=0.5, vmax=10.5)#myb
rgb = ls.shade_rgb(cmap(norm(mmidat)), topodat, blend_mode='hsv', vert_exag=1) 
im = m.imshow(rgb)
'''
##########################################################################################
# get land & lake polygons for masking
##########################################################################################
polys = get_map_polygons(m)

#mask_outside_polygon(polys[1][::-1], ax=None)
mask_outside_polygons(polys, 'lightskyblue', plt)

# get lake ploygons
polygons = []
for polygon in m.lakepolygons:
    poly = polygon.get_coords()
    plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
    polygons.append(poly)
    
##########################################################################################
# grid & plot MMI locs
##########################################################################################
# get colormap
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/mi_pop.cpt'
else:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//mi_pop.cpt'
ncols = 10
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))
    
# loop thru grid cells values and take average
avmmi = []
avlo  = []
avla  = []
for lo in lonrng:
    for la in latrng:
        idx = where((lon >= lo-d2) & (lon < lo+d2) & \
                    (lat >= la-d3) & (lat < la+d3))[0]
        
        # add to list greater than minObs
        if len(idx) >= minObs:
            avmmi.append(mean(mmi[idx]))
            avlo.append(lo)
            avla.append(la)
            
            # now plot
            pltx = [lo-d2, lo+d2, lo+d2, lo-d2, lo-d2]
            plty = [la+d3, la+d3, la-d3, la-d3, la+d3]
            
            x, y = m(pltx, plty)
            colidx = int(round(mean(mmi[idx])))
            c= tuple(cs[colidx][:-1])
            plt.fill(x, y, fc=c, ec='0.55')

"""
for i in range(0, len(mmi)): 
    #get colour idx
    
    colidx = int(round(mmi[i]))
    x, y = m(lon[i], lat[i])
    zo = sortidx[i] + 20
    #print i, mmi[i], colidx
    # set color tuple
    c= tuple(cs[colidx][:-1])
    plt.plot(x, y, 'o', markerfacecolor=c, markeredgecolor='k', markeredgewidth=0.25, \
             markersize=11, zorder=zo, alpha=1.)
"""
##########################################################################################
# add cities
##########################################################################################
'''
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
clat = [-31.95, -35.02, -33.65, -33.35, -33.36, -33.97, -31.65, -34.24, -34.96]
clon = [115.86, 117.88, 115.35, 115.66, 116.15, 115.05, 116.66, 116.15, 117.36] #144.]
cloc = ['Perth','Albany', 'Busselton', 'Bunbury', 'Collie', 'Margaret River', 'Northam', 'Manjimup', 'Denmark']
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='k', markeredgecolor='k', markeredgewidth=0.5, markersize=6, zorder=11000)

# label cities
off = 0.25
for i, loc in enumerate(cloc):
    if i == 1:
        x, y = m(clon[i]+0.03, clat[i]-0.03)
        plt.text(x, y, loc, size=14, ha='left', va='top', weight='normal', path_effects=path_effects, zorder=10000+i, fontsize=10)
    elif i == 8:
        x, y = m(clon[i]-0.05, clat[i]-0.03)
        plt.text(x, y, loc, size=14, ha='right', va='top', weight='normal', path_effects=path_effects, zorder=10000+i, fontsize=10)
    elif i == 0 or i == 0 or i == 3 or i == 7:
        x, y = m(clon[i]-0.05, clat[i]+0.03)
        plt.text(x, y, loc, size=14, ha='right', weight='normal', path_effects=path_effects, zorder=10000+i)
    else:
        x, y = m(clon[i]+0.03, clat[i]+0.03)
        plt.text(x, y, loc, size=14, ha='left', weight='normal', path_effects=path_effects, zorder=10000+i, fontsize=10)
'''
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]

# parse AU cities
cityFile = '../GA-shakemap/mapping/cities1000_au_ascii.txt'

lines = open(cityFile).readlines()

# label cities
clatList = []
clonList = []

# for ordering in plot priority
pltlo = []
pltla = []
pltpop = []
pltloc = []

for line in lines:
    dat = line.strip().split('\t')
    clat = float(dat[4])
    clon = float(dat[5])
    loc = dat[1]
    
    off = 0.15
    if clat > llcrnrlat and clat < urcrnrlat-0.15 \
       and clon > llcrnrlon and clon < urcrnrlon:
        
        # add city to list
        pltlo.append(clon)
        pltla.append(clat)
        pltpop.append(float(dat[14]))
        pltloc.append(dat[1])

# order locs in trems of poplation
sortidx = argsort(pltpop)[::-1] # and reverse

# now loop through ordered cities - use max 15
i = 0
for si in sortidx:
    # plt max 15 locs
    if i < 15:
        pltCity = True
        print pltloc[si]
        
        for clol, clal in zip(clonList, clatList):
            if abs(pltla[si] - clal) < pltbuffer and abs(pltlo[si] - clol) < (2.5*pltbuffer):
                #print pltloc[si]
                pltCity = False
                
        # build list of locs
        clatList.append(pltla[si])
        clonList.append(pltlo[si])
        minLatDiff = 0.
        minLonDiff = 0.
        
        if pltCity == True:
            x, y = m(pltlo[si], pltla[si])
            plt.plot(x, y, 'o', markerfacecolor='k', markeredgecolor='k', markeredgewidth=0.5, markersize=6, zorder=11000)
    
            x, y = m(pltlo[si]+txtoff, pltla[si]+txtoff)
            plt.text(x, y, pltloc[si], size=14, ha='left', weight='normal', path_effects=path_effects)
            
            i += 1
            
# make random marker for label
#plt.text(x, y, s, bbox=dict(facecolor='red', alpha=0.5))     

##########################################################################################
# annotate earthquake
##########################################################################################

x, y = m(eqlo, eqla)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=25, zorder=10000)

plt.title(' '.join((evid, 'M', mstr, place)), fontsize=20)
##########################################################################################
# make colorbar
##########################################################################################

# set colourbar
plt.gcf().subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.34,0.05,0.33,0.03]) # setup colorbar axes.

norm = mpl.colors.Normalize(vmin=0.5, vmax=10.5)#myb
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = range(1,11)
rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=19)
    
plt.savefig('_'.join((evid.replace('-',''), place.split(',')[0].replace(' ', '_'), \
            'gridded_mmi_data_gridded.png')), format='png', bbox_inches='tight', dpi=300)
plt.show()
