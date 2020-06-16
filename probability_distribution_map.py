# coding=utf-8

import os, sys
import numpy as np
import numpy.ma as ma
import glob
import matplotlib
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
import xarray as xr
from datetime import datetime
from netCDF4 import Dataset, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import animateScatter
import time
import utils
import config_marmine as confm
import fast_histogram
import cmocean
import laplacefilter
import mpl_util
from matplotlib.pyplot import cm 
import drawCircleOfDistance
import drawBathymetry

def get_paths(confobj):
    
   # polygons = get_polygons()        
    base= r'results' 
    if confobj.selectyear=='all':
        print("Creating probability map for all simulated years")
        return base+'/MarMine_{}_opendrift*_novertical.nc'.format(confobj.species) 
    else:
        print("Creating probability map for simulated year {}".format(confobj.selectyear))
        return base+'/MarMine_{}_opendrift_*{}_to_*{}_novertical.nc'.format(confobj.species,confobj.selectyear,confobj.selectyear+1) 
  #  return base+'/MarMine_{}_opendrift_{}_to_{}_novertical.nc'.format(confobj.species,confobj.startdate,confobj.enddate) 



def create_map(confobj):
    confobj.ax = plt.subplot(111)
    mymap = Basemap(llcrnrlon=confobj.xmin, llcrnrlat=confobj.ymin,
                    urcrnrlon=confobj.xmax, urcrnrlat=confobj.ymax,resolution='i', 
                    projection='merc')

   # mymap.drawmapboundary(fill_color='#677a7a')
   
    mymap.fillcontinents(color='lightgrey',zorder=2)
    mymap.drawcoastlines()
    confobj.mymap=mymap
    drawBathymetry.addBathymetry(confobj)

def get_pos(confobj):     
    print("Opening file {}".format(confobj.paths)) 
    df = xr.open_mfdataset(confobj.paths, concat_dim='trajectory') 
    d = df.groupby(df.trajectory) #.sel(time=slice(confobj.startdate, confobj.enddate))
    
    return df.lat[24:-1],df.lon[24:-1], df.z[24:-1], df.time[24:-1]

def distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295     #Pi/180
    a = 0.5 - np.cos((lat2 - lat1) * p)/2 + np.cos(lat1 * p) * np.cos(lat2 * p) * (1 - np.cos((lon2 - lon1) * p)) / 2
    return 12742 * np.arcsin(np.sqrt(a))

def createBins(confobj):

    print('func: createBins() => Creating bins for averaging')
   
    dy = distance(confobj.ymin,confobj.xmin,confobj.ymax,confobj.xmin)
    dx = distance(confobj.ymin,confobj.xmin,confobj.ymin,confobj.xmax)

    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))
    print("Distance from minimum to maximim latitude binned area is %s km"%(dy))
    
    degrees_to_radians = np.pi/180.0
    radians_to_degrees = 180.0/np.pi
    earth_radius=6371
    # Distance in km between new bins converted to degrees
    
    r = earth_radius*np.cos(confobj.ymin*degrees_to_radians)
    confobj.deltaX = (confobj.requiredResolution/r)*radians_to_degrees
    confobj.lon_bins = np.arange(np.floor(confobj.xmin),np.ceil(confobj.xmax),confobj.deltaX)
    confobj.dx=len(confobj.lon_bins)

    confobj.deltaY = (confobj.requiredResolution/earth_radius)*radians_to_degrees
    confobj.lat_bins = np.arange(np.floor(confobj.ymin),np.ceil(confobj.ymax),confobj.deltaY)
    confobj.dy=len(confobj.lat_bins)
    print('=> created binned array of domain of grid cell size (%s,%s) with resolution %s'%(confobj.deltaX,confobj.deltaY,confobj.requiredResolution))
    
    

def get_density(lats, lons, nlevels, confobj):

    bins=(confobj.dy, confobj.dx)
    ranges=((confobj.ymin,confobj.ymax), (confobj.xmin,confobj.xmax))
    bins = np.asarray(bins).astype(np.int64)
    ranges = np.asarray(ranges).astype(np.float64)
    edges = (np.linspace(*ranges[0,:], bins[0]+1),np.linspace(*ranges[1,:], bins[1]+1))
    print(np.shape(lats))

    density = fast_histogram.histogram2d(lats, lons, bins=bins, range=ranges)
    print("Done with density")
    H2, *ledges = np.histogram2d(lats, lons, bins=bins, range=ranges)
    print("Done with numpy")

  #  density = histogram2d(lats, lons, range=[[confobj.ymin,confobj.ymax],[confobj.xmin,confobj.xmax]], bins=[confobj.dy,confobj.dx])
    print(np.sum(density-H2))
    density=H2
    total = np.sum(density)
    density=(density/total)*100.
    print("Total number of points {} percentage sum {}".format(total,np.sum(density)))

    density = ma.masked_where(density == 0, density)
    levels = MaxNLocator(nbins=nlevels).tick_values(0.1,2)
    #norm .tick_values(density.min(),density.max())
   
    norm = BoundaryNorm(levels, ncolors=confobj.cmap.N, clip=True)

    # Turn the lon/lat of the bins into 2 dimensional arrays 
    lon_bins_2d, lat_bins_2d = np.meshgrid(confobj.lon_bins, confobj.lat_bins)
  
    return lon_bins_2d,lat_bins_2d,density,norm

def make_map(confobj):
    plt.clf()
    create_map(confobj)
    lats,lons,depths,times= get_pos(confobj)
    createBins(confobj)   
    confobj.results_startdate=times[0].values
    confobj.results_enddate=times[-1].values
    
    if confobj.plot_type == 'heatmap':
        nlevels = 10    
        confobj.cmap = plt.get_cmap(confobj.cmapname)     
     
        lon_bins_2d,lat_bins_2d,density,norm = get_density(lats[:,:], lons[:,:], nlevels, confobj)
        xs, ys = confobj.mymap(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh       
        
        cs = confobj.mymap.pcolormesh(xs, ys, density, cmap=confobj.cmap, norm=norm, edgecolors='face',alpha=0.8,linewidths=0.1)   

        # LOKI distribution/seed area
        X,Y = drawCircleOfDistance.createCircleAroundWithRadius(confobj.st_lats[0],confobj.st_lons[0],confobj.releaseRadius/1000.)
        confobj.mymap.plot(X,Y,latlon=True,marker=None,color='w',linewidth=0.2)

        # LOKI main point
        x,y = confobj.mymap(confobj.st_lons[0], confobj.st_lats[0])
        confobj.mymap.plot(x,y ,marker='D',color='w',markersize=0.1)
        plt.colorbar(cs, shrink=0.7)
        figname = r'MarMine_{}_{}_to_{}_heatmap.png'.format(confobj.species,confobj.results_startdate,confobj.results_enddate) 
     

    elif type == 'scatter':
        x,y = confobj.mymap(lons,lats)
        confobj.mymap.scatter(x,y,alpha = 0.5,c = 'k',s = 10)
        figname = r'MarMine_{}_{}_to_{}_scatter.png'.format(confobj.species,confobj.startdate,confobj.enddate) 
     
    plt.savefig(figname,format = 'png',dpi = 300)



def call_make_map(confobj):
    
    confobj.paths = get_paths(confobj) 
    print(confobj.paths)
    make_map(confobj)

if __name__ == "__main__":
    start_time = time.time()
    
    plotallyears=True
    allyears=[2015,2016,2017,2018,"all"]

    for selectyear in allyears:
        experiments = [1]
        confobj=confm.marmine_conf()
        confobj.selectyear=selectyear

        marine_organism="molusc"
        for experiment in experiments:

            confobj.experiment=experiment
            confobj.species=marine_organism
            confobj.plot_type='heatmap'
            call_make_map(confobj) 

        print("---  It took %s seconds to run the script ---" % (time.time() - start_time))   
