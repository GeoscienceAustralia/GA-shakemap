'''
script to map shakemap grid.xml data

Usage:
    python map_shakemap.py <path/to/grid.xml>
    
    e.g.:
        
    run map_shakemap.py ../../testing/grid.xml
'''

from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, mgrid, ogrid
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from mapping_tools import get_map_polygons, mask_outside_polygons
from shakemap_tools import parse_dataxml
from sys import argv

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse shakemap grid
##########################################################################################

xmlPath = argv[1]

event, gridspec, fields, xmldata = parse_dataxml(xmlPath)

# set data fields
xmllons = xmldata[:,0]
xmllats = xmldata[:,1]

# find MMI col
keys = fields.keys()
for key in keys:
    if fields[key] == 'mmi':
        mmicol = int(key[-1]) - 1

mmi = xmldata[:,mmicol]

##########################################################################################
# set up map
##########################################################################################
# get map extents
urcrnrlat = gridspec['lat_max']
llcrnrlat = gridspec['lat_min']
urcrnrlon = gridspec['lon_max']
llcrnrlon = gridspec['lon_min']

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
#m.drawmapboundary(fill_color='blue', zorder=50)
#m.fillcontinents(color='coral',lake_color='aqua')
m.drawparallels(arange(-90.,90.,0.5), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
m.drawmeridians(arange(0.,360.,0.5), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
'''
m.drawmapscale(llcrnrlat+0.2, llcrnrlat+0.1, llcrnrlat+0.2, llcrnrlat+0.1, 100, \
               fontsize = 16, barstyle='fancy', zorder=100)
'''
##########################################################################################
# read topo
##########################################################################################
mdiv = 500.

print 'Reading topo file...'
netcdffile = '../../../DATA/GEBCO/au_gebco.nc'
nc = NetCDFFile(netcdffile)

zscale =2. #gray
zscale = 40. #colour

'''
# if using SRTM 30m
data = nc.variables['z'][:] / zscale
lons = nc.variables['x'][:]
lats = nc.variables['y'][:]
'''
# if using GEBCO 30 arcsec
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/mdiv)+1
ny = int((m.ymax-m.ymin)/mdiv)+1
topodat = m.transform_scalar(data,lons,lats,nx,ny)

##########################################################################################
# plot intensity grid from shakemap
##########################################################################################

print 'Resampling data...'
N = 500j
extent = (llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat)

xs,ys = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
resampled = griddata(xmllons, xmllats, mmi, xs, ys, interp='linear')

# get 1D lats and lons for map transform
lons = ogrid[extent[0]:extent[1]:N]
lats = ogrid[extent[2]:extent[3]:N]

mmidat = m.transform_scalar(resampled.T,lons,lats,nx,ny)

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
# plot MMI locs
##########################################################################################
"""
from mmi_tools import parse_usgs_dyfi_geocoded
hnd = []

gammifile = 'GA_MMI_0.05d.csv'
lines = open(gammifile).readlines()
gammi = []
galat = []
galon = []
for line in lines:
    dat = line.strip().split(',')
    gammi.append(float(dat[2]))
    galon.append(float(dat[0]))
    galat.append(float(dat[1]))
    
x, y = m(galon, galat)
'''
for i, mmi in enumerate(gammi):
    h1 = plt.plot(x[i], y[i], marker=r"$ {} $".format(str(int(round(mmi)))), ms=10, color='r', mec='none')
    #plt.plot(x[i], y[i], marker=str(int(round(mmi))), ms=10, color='r')
'''    
h1 = plt.plot(x, y, 'kx', lw=2., ms=7)

melmmifile = 'MEL_MMI_0.05d.csv'
lines = open(melmmifile).readlines()
melmmi = []
mellat = []
mellon = []
for line in lines:
    dat = line.strip().split(',')
    melmmi.append(float(dat[2]))
    mellon.append(float(dat[0]))
    mellat.append(float(dat[1]))
    
x, y = m(mellon, mellat)
h2 = plt.plot(x, y, 'ko', markerfacecolor='None', lw=2., ms=7)

'''
dyfifile = 'usgs_geocoded_cdi.txt'
dyfidict = parse_usgs_dyfi_geocoded(dyfifile)
#dyfifile = 'usgs_zip_cdi.txt'
#dyfidict = parse_usgs_dyfi_zip(dyfifile)
dyfimmi  = []
dyfilat = []
dyfilon = []
for rec in dyfidict:
    if rec['nresp'] >= 1:
       dyfimmi.append(rec['cdi'])
       dyfilat.append(rec['lat'])
       dyfilon.append(rec['lon'])
'''
dyfimmifile = 'USGS_MMI_0.05d.csv'
lines = open(dyfimmifile).readlines()
dyfimmi = []
dyfilat = []
dyfilon = []
for line in lines:
    dat = line.strip().split(',')
    dyfimmi.append(float(dat[2]))
    dyfilon.append(float(dat[0]))
    dyfilat.append(float(dat[1]))

x, y = m(dyfilon, dyfilat)
h3 = plt.plot(x, y, 'k+', lw=2., ms=9)

plt.legend([h1[0], h2[0], h3[0]], ['AU MMI', 'MEL MMI', 'USGS DYFI'], \
            fontsize=15, loc=4, numpoints=1)
"""
##########################################################################################
# add cities
##########################################################################################
"""
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
clat = [-37.814, -38.197, -38.125, -38.235, -38.161, -38.102, -38.606] #-37.]
clon = [144.964, 146.541, 146.263, 146.395, 145.932, 147.065, 145.588] #144.]
cloc = ['Melbourne','Traralgon', 'Moe', 'Morwell', 'Warragul', 'Sale', 'Wonthaggi']
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2, markersize=11)

# label cities
off = 0.25
for i, loc in enumerate(cloc):
    if i == 3:
        x, y = m(clon[i]+0.025, clat[i]-0.025)
        plt.text(x, y, loc, size=18, ha='left', va='top', weight='normal', path_effects=path_effects)
    elif i == 4:
        x, y = m(clon[i]-0.048, clat[i]+0.025)
        plt.text(x, y, loc, size=18, ha='right', weight='normal', path_effects=path_effects)
    else:
        x, y = m(clon[i]+0.025, clat[i]+0.025)
        plt.text(x, y, loc, size=18, ha='left', weight='normal', path_effects=path_effects)
"""
##########################################################################################
# annotate earthquake
##########################################################################################

eqlat = event['lat']
eqlon = event['lon']

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=25)
                  
##########################################################################################
# add title
##########################################################################################

titlestr = ' '.join(('MW',str(event['magnitude']), event['event_timestamp'])) + '\n' \
                     + event['event_description']
plt.title(titlestr, fontsize=16)                     
##########################################################################################
# make colorbar
##########################################################################################

# set colourbar
plt.gcf().subplots_adjust(bottom=0.13)
cax = fig.add_axes([0.34,0.05,0.33,0.03]) # setup colorbar axes.

cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = range(1,11)
rom_num = ['I', 'II', 'III', 'VI', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=15)

##########################################################################################
# save file
##########################################################################################
    
pngFile = event['event_description'].split(',')[0].replace(' ','_')+'.png'
plt.savefig(pngFile, format='png', bbox_inches='tight', dpi=300)
plt.show()
