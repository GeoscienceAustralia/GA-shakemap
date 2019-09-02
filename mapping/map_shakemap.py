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
from mapping_tools import get_map_polygons, mask_outside_polygons, annotate_cities
from shakemap_tools import parse_dataxml
from sys import argv
from os import getcwd, mkdir, path
import shapefile

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')
cwd = getcwd()

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
lonspan = urcrnrlon - llcrnrlon

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
if lonspan <= 3.0:
    tickspace = 0.5
    scale_len = 50
    latoff = 0.12
    lonoff = 0.3
elif lonspan <= 4.0:
    tickspace = 0.5
    scale_len = 100
    latoff = 0.0
    lonoff = 0.0
elif lonspan > 4.0:
    tickspace = 1.0
    scale_len = 200
    latoff = 0.0
    lonoff = 0.0
m.drawparallels(arange(-90.,90.,tickspace), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
m.drawmeridians(arange(0.,360.,tickspace), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)

m.drawmapscale(llcrnrlon+lonoff, llcrnrlat+latoff, llcrnrlon+lonoff, llcrnrlat+latoff, scale_len, \
               fontsize = 14, barstyle='fancy', zorder=100)

##########################################################################################
# read topo
##########################################################################################

mdiv = 500.

print 'Reading topo file...'
netcdffile = '../../../DATA/GEBCO/au_gebco.nc'
nc = NetCDFFile(netcdffile)

zscale =2. #gray
if lonspan <= 4.0:
    zscale = 30. #colour
else:
    zscale = 30. #colour

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

# nas interprets grids differently
if cwd.startswith('/nas'):
    mmidat = m.transform_scalar(resampled.T,lons,lats,nx,ny)
else:
    mmidat = m.transform_scalar(resampled,lons,lats,nx,ny)

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
'''
gammifile = '../data/201605201814/raw/Pettermann_ranges_2016_intensities.csv'
lines = open(gammifile).readlines()[1:]
gammi = []
galat = []
galon = []
for line in lines:
    dat = line.strip().split(',')
    gammi.append(float(dat[2]))
    galon.append(float(dat[1]))
    galat.append(float(dat[0]))
    
x, y = m(galon, galat)
m.plot(x, y, 'o', mec='k', mfc='None', ms=10, lw=0.5, label='AU MMI')
plt.legend(fontsize=15, loc=4, numpoints=1)
'''
##########################################################################################
# add contours
##########################################################################################
'''
# resomple to smooth contours
N = 100j
xs2,ys2 = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
smoothed = griddata(xmllons, xmllats, mmi, xs2, ys2, interp='linear')

# set levels
levels = arange(2.5, 9.6, 1.)

# plt contours
x, y = m(xs2, ys2)
csm = plt.contour(x, y, smoothed, levels, colors='w', linewidth=3)
plt.clabel(csm, inline=1, fontsize=10, fmt='%0.1f')
'''

##########################################################################################
# add cities
##########################################################################################
numCities=15
annotate_cities(numCities, plt, m)
##########################################################################################
# annotate earthquake
##########################################################################################

eqlat = event['lat']
eqlon = event['lon']

x, y = m(eqlon, eqlat)
#plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=25)
                  
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
rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=15)

##########################################################################################
# save file
##########################################################################################
    
pngFile = event['event_description'].split(',')[0].replace(' ','_')+'.png'
plt.savefig(pngFile, format='png', bbox_inches='tight', dpi=300)

##########################################################################################
# make shapefile of contour lines
##########################################################################################

# check to see if shapefile contours exists
if path.isdir('contours') == False:
    mkdir('contours')
    
# make list of levels - old levels array([0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.12, 0.18, 0.24])
allLevels = [arange(0.5, 11., 1.)]
 
levelNames = ['mmi_contours'] #, 'lev_0_01', 'lev_0_02', 'lev_0_05'


# resomple to smooth contours
N = 80j
xs2,ys2 = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
smoothed = griddata(xmllons, xmllats, mmi, xs2, ys2, interp='linear')

resampled = smoothed
xs = xs2
ys = ys2          

# loop thru levels
for levels, levelName in zip(allLevels, levelNames):
    
    # setup shapefile
    '''
    outshp = path.join('contours', '_'.join((levelNames[0].replace(' ','_'), key, \
                       levelName, 'contours.shp')))
    '''
    outshp = 'mmi_contours.shp'

    # set shapefile to write to
    w = shapefile.Writer(shapefile.POLYLINE)
    w.field('MMI','F', 5, 2)
        
    # have to re-contour using un-transformed lat/lons
    cs = plt.contour(xs, ys, resampled, levels, colors='k')
    
    plt.close(fig)
    
    # loop through contour levels
    for l, lev in enumerate(cs.levels):
        contours = cs.collections[l].get_paths()
        
        # now loop through multiple paths within level
        if len(contours) > 0:
            for cnt in contours:
                
                # add polyline to shapefile
                w.line(parts=[cnt.vertices], shapeType=shapefile.POLYLINE)
                
                # add level attribute
                w.record(lev)

    # now save area shapefile
    w.save(outshp)
    
    # write projection file
    prjfile = outshp.strip().split('.shp')[0]+'.prj'
    f = open(prjfile, 'wb')
    f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
    f.close()

plt.show()
