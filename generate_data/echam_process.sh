#!/bin/bash
#
# This script basically uses nco package to prepare echam model output for lagranto (daniel.boateng@uni-tuebingen.de)
#--------------------------------------------------------------
# It requires that nco is installed on your system. 

#path to data 

mpath='/home/dboateng/Model_output_pst/a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h/output_processed/6h_MEANS'
filename='a001_1003_1017_07_lterm.01.nc'

#inverting latitudes

ncpdq -O -a -lat ${mpath}/${filename} ${mpath}/${filename}

# rotating of longitudes 

ncks -O -msa -d lon,180.,360. -d lon,0.,179. ${mpath}/${filename} ${mpath}/${filename} 

# changing lats to -180 to 180

ncap2 -O -s 'where(lon>=180.0) lon = lon-360.0' ${mpath}/${filename} ${mpath}/${filename}
