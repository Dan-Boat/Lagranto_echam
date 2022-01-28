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
data_path = "/home/dboateng/source_codes/lagranto/new/a001/June/Trace/Stuttgart"

files = ["wcb_1.1", "wcb_2.1"]
filename = os.path.join(data_path, files[0])
trajs = Tra(filename)
new_trajs = [Tra(os.path.join(data_path, f)) for f in files]
all_trajs = trajs.concatenate(new_trajs)

trajs1 = Tra()
trajs2 = Tra()
trajs3 = Tra()
trajs4 = Tra()
trajs5 = Tra()
trajs6 = Tra()
trajs7 = Tra()
trajs8 = Tra()
trajs9 = Tra()
trajs10 = Tra()
trajs0 = Tra()
# trajs12 = Tra()

traj1 = Tra()
traj2 = Tra()
traj3 = Tra()
traj4 = Tra()
traj5 = Tra()
traj6 = Tra()
traj7 = Tra()
traj8 = Tra()
traj9 = Tra()
traj10 = Tra()


# 
# load the trajectories
#trajs1.load_ascii('/home/dboateng/source_codes/lagranto/new/a003/June/wcb.1')
trajs0.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/Stuttgart/wcb_1.1')
trajs1.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/Stuttgart/wcb_2.1')


# trajs2.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_2.1')
# trajs3.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_3.1')
# trajs4.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_4.1')
# trajs5.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_5.1')
# trajs6.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_6.1')
# trajs7.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_7.1')
# trajs8.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_8.1')
# trajs9.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_9.1')
# trajs10.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_10.1')


# traj1.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_1.1')
# traj2.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_2.1')
# traj3.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_3.1')
# traj4.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_4.1')
# traj5.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_5.1')
# traj6.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_6.1')
# traj7.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_7.1')
# traj8.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_8.1')
# traj9.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_9.1')
# traj10.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/July/Trace/wcb_10.1')

#trajs5.load_ascii('/home/dboateng/source_codes/lagranto/new/a001/June/Trace/wcb_5.1')
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
    cont = m.plot_traj(trajs6, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs7, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs8, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs9, variables[i], cmap = 'Reds')
    cont = m.plot_traj(trajs10, variables[i], cmap = 'Reds')
    
    ont = m.plot_traj(traj1, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj2, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj3, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj4, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj5, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj6, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj7, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj8, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj9, variables[i], cmap = 'Blues')
    cont = m.plot_traj(traj10, variables[i], cmap = 'Blues')

    
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
    
    