#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.models.oceandrift import OceanDrift
from IBM.pelagicplankton import PelagicPlanktonDrift
from pprint import pprint
from netCDF4 import Dataset, date2num, num2date
from scipy.ndimage.filters import gaussian_filter
import matplotlib
import os

o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information

#######################
# Preparing readers
#######################

base = 'results'
baseout = 'figures'
hexagon = False
startdate = '01062015'
enddate = '30122015'
experiment = 1
filename="results/MarMine_molusc_opendrift_"+str(startdate)+"_to_"+str(enddate)+"_novertical.nc"
print(filename)
reader_basemap = reader_basemap_landmask.Reader(
                       llcrnrlon=2, llcrnrlat=65.,
                       urcrnrlon=20, urcrnrlat=75,
                       resolution='i', projection='merc')
o.add_reader(reader_basemap)

plotfilenameAnime = 'Figures/test.mp4'
plotfilenameColor = 'Figures/test_color.png'
plotfilename = 'Figures/test.png'

if os.path.exists(filename):
    print("=> Opening input file: {}".format(filename)
    )
    o.io_import_file(filename)
    # o.add_reader([reader_basemap]) #Do not include basemap when stranding is deactivated

   # o.plot_vertical_distribution()
    o.plot(linecolor='z',lvmin=-150, lvmax=0,filename=plotfilenameColor)
    o.plot(filename=plotfilename)
    o.animation(filename=plotfilenameAnime)
else:
    print("=> File does not exist {}".format(filename))