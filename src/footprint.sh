#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp footprint short
  echo  
  exit 0
endif

set inpfile = $1
set outfile = $2
set mode    = $3
if ( "${mode}" == "proxy" ) then
   set clon   = $4
   set clat   = $5
   set radius = $6
endif

set dim=`${LAGRANTO}/goodies/trainfo.sh ${inpfile} dim` 

\rm -f traj2num.param
echo \"${inpfile}\"  >! footprint.param
echo \"${outfile}\"  >> footprint.param
echo ${dim}          >> footprint.param
echo \"${mode}\"     >> footprint.param

if ( "${mode}" == "proxy" ) then
  echo ${clon}      >> footprint.param
  echo ${clat}      >> footprint.param
  echo ${radius}    >> footprint.param
endif

${LAGRANTO}/tools/footprint

\rm -f footprint.param

exit 0

