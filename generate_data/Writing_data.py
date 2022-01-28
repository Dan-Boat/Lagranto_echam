#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 11:15:31 2021

@author: dboateng
This script generates the individual data required for trajectory analysis with lagranto package. Note that the package 
require that all the individaul time point variables must be in a singe netcdf file. Therefore the processed echam data 
can be pass the this script which would generate them depending on the user declarations. Therefore, understanding of the 
function is required. Pls read the documentation of the variables or input parameters
"""

#import modules 
import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os 





def echam4lagranto(path_to_data, filename, month=None, season=None, day=None, year=2000,
                   hrs=[0, 6, 12, 18], path_to_save=None, month_num=None, season_num=None,
                   start_date="2000-06-01"):
    """
    

    Parameters
    ----------
    path_to_data : TYPE:str
        DESCRIPTION. The path to where the processed data is stored
    filename : TYPE: str
        DESCRIPTION. The filename of the processed echam data
    month : TYPE: int, optional
        DESCRIPTION. The default is None. The number for the month to generate
    season : TYPE: int, optional
        DESCRIPTION. The default is None. The number of the season to generate
    day : TYPE: int, optional
        DESCRIPTION. The default is None. The number for day to generate (can be defined 
                                                                          in a loop)
    year : TYPE: int, optional
        DESCRIPTION. The default is 2000. Use the default because the langranto has issues with time 
        range outside the Gregorian. Since the individual time is in single files, using synthetic year 
        would not affect the calculations
    hrs : TYPE: list, optional
        DESCRIPTION. The default is [0, 6, 12, 18]. Use the default because of the 6h output from echam
    path_to_save : TYPE: str, optional
        DESCRIPTION. The default is None. path to store the data to generate (must be in the lagranto folder)
    month_num : TYPE: int, optional
        DESCRIPTION. The default is None. Month number for saving the data
    season_num : TYPE: int, optional
        DESCRIPTION. The default is None. Season number for filename
    start_date : TYPE: date, optional
        DESCRIPTION. The default is "2000-06-01". Start date for the synthetic time 

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    Netcdf file in your path to save directory

    """
    
    data = xr.open_dataset(os.path.join(path_to_data, filename), decode_cf=True, use_cftime=True)
    
    # time for replacing year
    period = len(data.time)
    time = xr.cftime_range(start=start_date, periods=period, freq="6H", calendar="gregorian")
    
    data["time"] = time
    
    #rename mlev to lev
    data = data.rename({"mlev":"lev"})
    
    #rename variables
    if all(hasattr(data, attr) for attr in ["lsp", "aps", "aprl", "aprc", "u", "v", "omega", "svo", "relhum",
                                            "geopoth", "st", "q"]):
        
        data = data.rename({"lsp":"LSP", "st":"T", "aps":"PS", "omega":"OMEGA", "u":"U", "v":"V",
                                        "geopoth":"z", "svo":"SVO", "relhum":"RH", "aprl":"LP", "aprc":"CP", "q":"Q"})
    elif all(hasattr(data, attr) for attr in ["lsp", "aps", "aprl", "aprc", "u", "v", "omega", "svo", "relhum",
                                            "geopoth", "st"]):
        
        data = data.rename({"lsp":"LSP", "st":"T", "aps":"PS", "omega":"OMEGA", "u":"U", "v":"V",
                                        "geopoth":"z", "svo":"SVO", "relhum":"RH", "aprl":"LP", "aprc":"CP"})
    else: 
        print("More variables or less are not detected in the meta data")
        
        
    # converting surface presure to hPa
    data["PS"] = data.PS / 100 # Pa --> hPa
    data["PS"].attrs["units"] = "hPa"
    data["PS"].attrs["long_name"] = "Surface pressure"
    data["PS"].attrs["code"] = "143"
    data["PS"].attrs["table"] = "128"
    
    
    if month is not None:
        data_m = data.groupby(group="time.month", restore_coord_dims=True)
        data_m = data_m[month].groupby(group="time.day")
    elif month is None:
        if season is None:
            data_m = data.groupby(group="time.day", restore_coord_dims=True)
        else:
            data_m = data.groupby(group="time.season", restore_coord_dims=True)
        
        if day is not None:
            data_d = data_m[day].groupby(group="time.hour")
        else:
            raise ValueError("day must be defined b/n 1-30")
        
        for i,hr in enumerate(hrs):
            if path_to_save is not None:
                if month_num < 10:
                    mname= "0" +str(month_num)
                else:
                    mname = str(month_num)
                if day < 10:
                    dname = "0" +str(day)
                else:
                    dname = str(day)
                if hr < 10:
                    hname = "0" +str(hr)
                else:
                    hname = str(hr)
                data_d[hr].to_netcdf(os.path.join(path_to_save, "P"+ str(year)+ mname + dname +"_"+ hname))
                data_d[hr].to_netcdf(os.path.join(path_to_save, "S"+ str(year)+ mname + dname +"_"+ hname))
            
            else:
                print("saving in the path of runing script")
                data_d[hr].to_netcdf("P"+str(year)+str(month_num)+str(day)+"_"+str(hr))
                data_d[hr].to_netcdf("S"+str(year)+str(month_num)+str(day)+"_"+str(hr))


###############################################################################
#Generating data for 10 days using a monthly input data (This can be adjusted depending on the 
# prcessed echam data)   
# defining paths 

path_to_echam = "/home/dboateng/Model_output_pst/a002_hpc-bw_e5w2.3_t159_PI_Alps_east_100_t159l31.6h/output_processed/6h_MEANS"
path_to_store = "/home/dboateng/source_codes/lagranto/new/a002/August"
filename = "a002_1003_1017_08_lterm.01.nc"
         
for i in range(1,30,1):
    print ("extracting for day", i)
    echam4lagranto(path_to_data=path_to_echam, filename=filename, month=None, day=i, path_to_save=path_to_store, month_num=8, 
                  start_date="2000-08-01", year=2000)
    
path_to_store = "/home/dboateng/source_codes/lagranto/new/a002/July"
filename = "a002_1003_1017_07_lterm.01.nc"
         
for i in range(1,30,1):
    print ("extracting for day", i)
    echam4lagranto(path_to_data=path_to_echam, filename=filename, month=None, day=i, path_to_save=path_to_store, month_num=7, 
                  start_date="2000-07-01", year=2000)
    
path_to_store = "/home/dboateng/source_codes/lagranto/new/a002/June"
filename = "a002_1003_1017_06_lterm.01.nc"
         
for i in range(1,30,1):
    print ("extracting for day", i)
    echam4lagranto(path_to_data=path_to_echam, filename=filename, month=None, day=i, path_to_save=path_to_store, month_num=6, 
                  start_date="2000-06-01", year=2000)
    