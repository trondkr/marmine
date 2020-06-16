#
#
# Config object for MarMine
import time, calendar
from datetime import datetime, timedelta
import numpy as np

__author__ = 'Trond Kristiansen'
__email__ = 'me (at) trondkristiansen.com'
__created__ = datetime(2019, 7, 19)
__modified__ = datetime(2019, 7, 19)
__version__ = "1.0"
__status__ = "Development, modified on 19.07.2019"


class marmine_conf():

    def user_defined_inits(self):
        self.experiment=1
        self.plot_type='heatmap'
        self.requiredResolution=30 # km between bins

        self.xmin=-22.0
        self.xmax=23.5
        self.ymin=65.0
        self.ymax=78.0
        self.releaseParticles=8000
        self.releaseRadius=15000
        self.lowDepth=-1000
        self.highDepth=-950
        self.verticalBehavior=False

        if self.experiment == 1: 
            self.startdate=datetime(2014,6,1,0,0,0)
            self.enddate=datetime(2018,6,1,0,0,0)

        self.basedir='/Volumes/DATASETS/NORDIC4KM/netcdf4/' 
        years=(np.linspace(self.startdate.year% 10,self.enddate.year% 10, (self.enddate.year-self.startdate.year)+1,endpoint=True)).astype(int)
        self.pattern='roms_nordic4.an.24h.201[{}-{}]*'.format(years[0],years[-1])

        self.species='molusc'
        self.plot_type='heatmap'
        self.cmapname='RdYlBu_r'
        self.selectyear='all'
        # LOKI station - seed locations
        self.st_lons = [8.1]
        self.st_lats = [73.5]
    

    def define_orgamisms(self):
        return {'molusc': ['temperature', 'salinity', 'ssh', 'uvel', 'vvel']}

    def __init__(self):
        print('\n--------------------------\n')
        print('Started ' + time.ctime(time.time()))

        self.user_defined_inits()
        
        self.paths=None
        self.mymap=None
        self.ax=None
        self.deltaX = None
        self.deltaY = None
        self.dx=None
        self.dy=None
        self.cmap=None
        self.outputFilename=None
        self.results_startdate=None
        self.results_enddate=None

