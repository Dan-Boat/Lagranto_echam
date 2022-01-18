# test script for runing backtrajectory 
# defining paths to move data intp src 
mpath="/home/dboateng/source_codes/lagranto/new/a002"
month="August"
lagranto="/home/dboateng/source_codes/lagranto/new"
#
startdate="20000821_18"
enddate="20000811_18"
#
# export path 
export LAGRANTO="/home/dboateng/source_codes/lagranto/new"

startdates=("20000821_18"
"20000820_18"
"20000819_18"
"20000818_18"
"20000817_18"
"20000816_18"
"20000815_18"
"20000814_18"
"20000813_18"
"20000812_18"
"20000811_18")

enddates=("20000811_18"
"20000810_18"
"20000809_18"
"20000808_18"
"20000807_18"
"20000806_18"
"20000805_18"
"20000804_18"
"20000803_18"
"20000802_18"
"20000801_18")
#
# move files to lagranto
mv ${mpath}/${month}/P* $lagranto
mv ${mpath}/${month}/S* $lagranto

# create starting points 
# every grid [12E, 17E, 46N, 48.5N]
# every  hPa

#csh create_startf.sh $startdate startf.2 'box.grid(12,17,46,48.5) @list(850)@hPa' -changet
csh create_startf.sh ${startdate} startf.2 'point(11.57,48.15) @list(850)@hPa' -changet

# calculate trajectories. -j flag prevents stops at the ground level
csh caltra.sh $startdate $enddate startf.2 traj.4 -j  

#loop through the startdates

for date in {0..5}; do
    echo ${startdates[dates]}
    csh create_startf.sh ${startdates[date]} startf.2 'point(11.57,48.15) @list(850)@hPa' -changet
    csh caltra.sh ${startdates[date]} ${enddates[date]} startf.2 traj_$date.4 -j 
    #csh mergetra.sh traj_i.4 traj.4 traj.4
    csh select.sh traj_$date.4 wcb_$date.1 "IN:LAT:-30,80:ALL"
    csh trace.sh wcb_$date.1 wcb_$date.1

done 
# select the trajectory within certain coordinates
csh select.sh traj.4 wcb.1 "IN:LAT:-30,80:ALL"

# tracing additional variables aside pressure (variables must be defined in tracevars in the runing directory)
csh trace.sh wcb.1 wcb.1

# move files back to the original directory
echo "--------cleaning up ......"
mv P* ${mpath}/${month}
mv S* ${mpath}/${month}

mv wcb.1 traj.4 startf.2 ${mpath}/${month}

echo "----------------creating directory--------------"
mkdir ${mpath}/${month}/Trace
mkdir ${mpath}/${month}/Traj

echo "----------moving files ---------------"
mv wcb_* ${mpath}/${month}/Trace
mv traj_* ${mpath}/${month}/Traj

echo "-------DOne-------"
