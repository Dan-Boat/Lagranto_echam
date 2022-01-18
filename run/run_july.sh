# test script for runing backtrajectory 
# defining paths to move data intp src 
mpath="/home/dboateng/source_codes/lagranto/new/a002"
month="July"
lagranto="/home/dboateng/source_codes/lagranto/new"
#
startdate="20000721_18"
enddate="20000711_18"
LOC=Stuttgart
#
# export path 
export LAGRANTO="/home/dboateng/source_codes/lagranto/new"

startdates=("20000721_18"
"20000720_18"
"20000719_18"
"20000718_18"
"20000717_18"
"20000716_18"
"20000715_18"
"20000714_18"
"20000713_18"
"20000712_18"
"20000711_18")

enddates=("20000711_18"
"20000710_18"
"20000709_18"
"20000708_18"
"20000707_18"
"20000706_18"
"20000705_18"
"20000704_18"
"20000703_18"
"20000702_18"
"20000701_18")
#
# move files to lagranto
mv ${mpath}/${month}/P* $lagranto
mv ${mpath}/${month}/S* $lagranto

# create starting points 
# every grid [12E, 17E, 46N, 48.5N]
# every  hPa

#csh create_startf.sh $startdate startf.2 'box.grid(12,17,46,48.5) @list(850)@hPa' -changet
csh create_startf.sh ${startdate} startf.2 'point(9.09,48.77) @list(850)@hPa' -changet

# calculate trajectories. -j flag prevents stops at the ground level
csh caltra.sh $startdate $enddate startf.2 traj_0.4 -j  

#loop through the startdates

for date in {1..10}; do
    echo ${startdates[dates]}
    csh create_startf.sh ${startdates[date]} startf.2 'point(9.09,48.77) @list(850)@hPa' -changet
    csh caltra.sh ${startdates[date]} ${enddates[date]} startf.2 traj_$date.4 -j 
    #csh mergetra.sh traj_i.4 traj.4 traj.4
    csh select.sh traj_$date.4 wcb_$date.1 "IN:LAT:-30,80:ALL"
    csh trace.sh wcb_$date.1 wcb_$date.1

done 
# select the trajectory within certain coordinates
csh select.sh traj_0.4 wcb_0.1 "IN:LAT:-30,80:ALL"

# tracing additional variables aside pressure (variables must be defined in tracevars in the runing directory)
csh trace.sh wcb_0.1 wcb_0.1

# move files back to the original directory
echo "--------cleaning up ......"
mv P* ${mpath}/${month}
mv S* ${mpath}/${month}

mv wcb_0.1 traj_0.4 startf.2 ${mpath}/${month}

echo "----------------creating directory--------------"
mkdir ${mpath}/${month}/Trace/${LOC}
mkdir ${mpath}/${month}/Traj/${LOC}

echo "----------moving files ---------------"
mv wcb_* ${mpath}/${month}/Trace/${LOC}
mv traj_* ${mpath}/${month}/Traj/${LOC}

echo "-------DOne-------"
