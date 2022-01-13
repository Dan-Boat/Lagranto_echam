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
trajs4 = Tra()
trajs5 = Tra()
trajs0 = Tra()
# trajs7 = Tra()
# trajs8 = Tra()
# trajs9 = Tra()
# trajs10 = Tra()
# trajs11 = Tra()
# trajs12 = Tra()

# 
# load the trajectories
#trajs1.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/June/wcb.1')
trajs0.load_ascii('/home/dboateng/source_codes/lagranto/new/a002/August/Trace/wcb_0.1')
trajs1.load_ascii('/home/dboateng/source_codes/lagranto/new/a002/August/Trace/wcb_1.1')
trajs2.load_ascii('/home/dboateng/source_codes/lagranto/new/a002/August/Trace/wcb_2.1')
trajs3.load_ascii('/home/dboateng/source_codes/lagranto/new/a002/August/Trace/wcb_3.1')
trajs4.load_ascii('/home/dboateng/source_codes/lagranto/new/a002/August/Trace/wcb_4.1')
trajs5.load_ascii('/home/dboateng/source_codes/lagranto/new/a002/August/Trace/wcb_5.1')
# trajs7.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/August/Trace/wcb_7.1')
# trajs8.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/August/Trace/wcb_8.1')
# trajs9.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/August/Trace/wcb_9.1')

# trajs10.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/August/wcb.1')
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
    cont = m.plot_traj(trajs1, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs2, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs3, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs4, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs5, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs0, variables[i], cmap = 'Reds')

    
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
    
    