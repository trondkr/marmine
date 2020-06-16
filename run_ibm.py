from datetime import datetime, timedelta
import numpy as np
from opendrift.readers import reader_basemap_landmask
from IBM.pelagicplankton import PelagicPlanktonDrift

from opendrift.readers import reader_netCDF_CF_NORDIC4KM
import logging
import os
from netCDF4 import Dataset, date2num,num2date
from numpy.random import RandomState
import matplotlib.pyplot as plt
import time
import config_marmine as confm
import random

__author__ = 'Trond Kristiansen'
__email__ = 'me (at) trondkristiansen.com'
__created__ = datetime(2019, 7, 26)
__modified__ = datetime(2019, 7, 26)
__version__ = "1.0"
__status__ = "Development, modified on 25.06.2018, 06.08.2018, 26.07.2019"

def createOutputFilenames(confobj):
    startDate=''
    if confobj.startdate.day<10:
        startDate+='0%s'%(confobj.startdate.day)
    else:
        startDate+='%s'%(confobj.startdate.day)

    if confobj.startdate.month<10:
        startDate+='0%s'%(confobj.startdate.month)
    else:
        startDate+='%s'%(confobj.startdate.month)

    startDate+='%s'%(confobj.startdate.year)

    endDate=''
    if confobj.enddate.day<10:
        endDate+='0%s'%(confobj.enddate.day)
    else:
        endDate+='%s'%(confobj.enddate.day)

    if confobj.enddate.month<10:
        endDate+='0%s'%(confobj.enddate.month)
    else:
        endDate+='%s'%(confobj.enddate.month)

    endDate+='%s'%(confobj.enddate.year)
 
    # File naming
    if confobj.verticalBehavior:
        outputFilename='results/MarMine_%s_opendrift_%s_to_%s_vertical.nc'%(confobj.species,startDate,endDate)
    else:
        outputFilename='results/MarMine_%s_opendrift_%s_to_%s_novertical.nc'%(confobj.species,startDate,endDate)
    if not os.path.exists('results'):
        os.makedirs('results')
    if os.path.exists(outputFilename):
        os.remove(outputFilename)
        
    confobj.outputFilename=outputFilename

   
def createAndRunSimulation(confobj):

    # Setup a new simulation
    o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information
    o.complexIBM=False

    reader_basemap = reader_basemap_landmask.Reader(llcrnrlon=confobj.xmin, llcrnrlat=confobj.ymin,
                    urcrnrlon=confobj.xmax, urcrnrlat=confobj.ymax,resolution='i', 
                    projection='merc')
                       
    o.add_reader([reader_basemap]) 
    print(confobj.basedir+confobj.pattern)
    reader_physics = reader_netCDF_CF_NORDIC4KM.Reader(confobj.basedir+confobj.pattern)
    o.add_reader([reader_physics ]) 

    #######################
    #Adjusting configuration
    #######################
    if confobj.verticalBehavior:
        o.set_config('processes:turbulentmixing', False)
    else:
        o.set_config('processes:turbulentmixing',  False)
    
    o.set_config('processes:verticaladvection', False)
    o.set_config('turbulentmixing:diffusivitymodel','windspeed_Sundby1983')
    
    o.set_config('turbulentmixing:TSprofiles', False)
    o.set_config('drift:scheme', 'runge-kutta4')
    o.set_config('drift:max_age_seconds', 7776000)

    o.set_config('general:coastline_action', 'previous') #Prevent stranding, jump back to previous position
    o.set_config('general:basemap_resolution', 'c')


    #######################
    # IBM configuration   
    #######################
    o.set_config('biology:constantIngestion', 0.75)
    o.set_config('biology:activemetabOn', 1)
    o.set_config('biology:cod', True)
    o.set_config('biology:haddock', False)
    o.set_config('biology:attenuationCoefficient',0.18)
    if confobj.verticalBehavior:
        o.set_config('biology:fractionOfTimestepSwimming',0.15) # Pause-swim behavior
    else:
        o.set_config('biology:fractionOfTimestepSwimming',0.00) # Pause-swim behavior
    o.set_config('biology:lowerStomachLim',0.3) #Min. stomach fullness needed to actively swim down
    

    z_levels=random.randrange(confobj.lowDepth,confobj.highDepth,confobj.releaseParticles) 
    print(z_levels)
    print('Seeding {} elements within a radius of {} m:'.format(confobj.releaseParticles, confobj.releaseRadius))
   
    print("Releasing {} larvae between {} and {}".format(confobj.species,confobj.startdate,confobj.enddate))
    o.seed_elements(lon=confobj.st_lons, 
                    lat=confobj.st_lats, 
                    number=confobj.releaseParticles, 
                    radius=[confobj.releaseRadius], cone=False, 
                    time=[confobj.startdate,confobj.enddate], 
                    z=z_levels)
    
    print('Elements scheduled for {} : {}'.format(confobj.species,o.elements_scheduled))

    o.run(end_time=confobj.enddate, 
          time_step=timedelta(hours=1),
          time_step_output=timedelta(hours=2),
          outfile=confobj.outputFilename)
    print(o)


if __name__ == "__main__":
    start_time = time.time()

    confobj=confm.marmine_conf()
    
    experiments = [1]
    
    for experiment in experiments:
        confobj.experiment=experiment
        createOutputFilenames(confobj)
        print("Result files will be stored as:\nnetCDF=> {}".format(confobj.outputFilename))
        createAndRunSimulation(confobj)

    print("---  It took %s seconds to run the script ---" % (time.time() - start_time)) 