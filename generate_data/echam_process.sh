#!/bin/bash
#
# This script basically uses nco package to prepare echam model output for lagranto (daniel.boateng@uni-tuebingen.de)
#--------------------------------------------------------------
# It requires that nco is installed on your system. 

#inverting latitudes

ncpdq -O -a -lat echam_input output

# rotating of longitudes 

ncks -O -msa -d lon,180.,360. -d lon,0.,179. echam_input output 

# changing lats to -180 to 180

ncap2 -O -s 'where(lon>=180.0) lon = lon-360.0' echam_input output
