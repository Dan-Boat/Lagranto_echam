#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 19:13:05 2021

@author: dboateng
"""

import numpy as np
import matplotlib
import matplotlib.pylab as plt  
from dypy.plotting import Mapfigure
from mpl_toolkits.basemap import Basemap
import os 

# initialize the variable for reading trajectories from an ascii file
from dypy.lagranto import Tra
trajs1 = Tra()
trajs2 = Tra()
trajs3 = Tra()

# 
# load the trajectories
trajs1.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/June/wcb.1')
trajs2.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/July/wcb.1')
trajs3.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/August/wcb.1')
path_to_store = "/home/dboateng/source_codes/lagranto/new/a003/"
# list of variables for plotting
variables =["p"]
#variables = ['p', 'Q', 'T', 'TH']
#titles = ['Elevation, hPa', 'Specific humidity, g/kg', 'Temperature, K', 'Potential temperature, K']
titles = ['Elevation, hPa'] #'Temperature, K', "Relative Humidity"]
for i in range(len(variables)):
    # set up the map and projection
    fig = plt.figure(figsize = (10,8))
    m = Mapfigure(width=9500000,height=6500000,
            resolution='l',projection='laea',\
            lat_ts=50,lat_0=52,lon_0=0.)
    m.drawmap(nbrem=10, nbrep=10)    
    # plot the trajectories
    cont = m.plot_traj(trajs1, variables[i], cmap = 'Blues')
    cont = m.plot_traj(trajs2, variables[i], cmap = 'Greys')
    cont = m.plot_traj(trajs3, variables[i], cmap = 'Reds')
    # display the colorbar
    cb = fig.colorbar(cont)    
    # set font size in captions
    for font_objects in cb.ax.yaxis.get_ticklabels():
        font_objects.set_size(20)
    # set figure title
    plt.title(titles[i], fontsize=24)
    # save the figure
    figname = variables[i] + '.svg'
    fig.savefig(os.path.join(path_to_store, figname), bbox_inches='tight')
    
    