from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour

plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
#108/152/-44/-8
urcrnrlat = -8.
llcrnrlat = -44.
urcrnrlon = 152.
llcrnrlon = 107.5
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='i',area_thresh=2000.)

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print 'Reading netCDF file...'
#nc = NetCDFFile('//Users//tallen//Documents//DATA//GMT//GEBCO//Australia_30c.nc')
#nc = NetCDFFile('//Users//tallen//Documents//DATA//GMT//GEBCO//Australia_30c.nc')
#nc = NetCDFFile('usgs_vs30.grd') # usgs
#nc = NetCDFFile('asscm_wii_vs30.grd.hold') # non-modified
nc = NetCDFFile('asscm_wii_vs30.900.grd') # hardwired to 900?


zscale =20. #gray
zscale =50. #colour
zscale = 1. 
data = nc.variables['z'][:] / zscale
try:
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]
except:
    lons = nc.variables['x'][:]
    lats = nc.variables['y'][:]
    
# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print 'Getting colormap...'
# get colormap
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//mby_topo-bath_mod.cpt'
cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt//YlOrRd_09.cpt'
cmap, zvals = cpt2colormap(cptfile, 10)
cmap = remove_last_cmap_colour(cmap)
#cmap = cm.get_cmap('terrain', 256)

# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 45)
norm = mpl.colors.Normalize(vmin=100/zscale, vmax=1100/zscale)#myb
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)
"""
##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,
ticks = arange(0, 19, 2)
years = float(minyear) + ticks*10
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)
"""
cax = fig.add_axes([0.77,0.3,0.02,0.4]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=1., norm=norm) # 
#cb.set_ticks(years)
#cb.set_ticklabels(labels)
cb.set_label('Vs30 (m/s', rotation=270, labelpad=20, fontsize=15)


plt.savefig('aus_vs30_map.png', format='png', bbox_inches='tight') #, dpi=300)
plt.show()
