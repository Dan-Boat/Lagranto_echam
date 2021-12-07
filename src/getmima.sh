#!/bin/csh

# Write usage information
if ( ${#argv} == 0) then
  echo 
  ${LAGRANTO}/bin/lagrantohelp getmima short
  echo  
  exit 0
endif

${LAGRANTO}/getmima $*

exit 0

