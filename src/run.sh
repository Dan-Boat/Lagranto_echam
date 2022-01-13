# test script for runing backtrajectory 
# defining paths to move data intp src 
mpath="/home/dboateng/source_codes/lagranto/new/a003"
month="August"
lagranto="/home/dboateng/source_codes/lagranto/new"
#
startdate="20000815_18"
enddate="20000809_18"
#
# export path 
export LAGRANTO="/home/dboateng/source_codes/lagranto/new"

startdates=("20000815_12"
"20000815_06"
"20000815_00"
"20000814_18"
"20000814_12"
"20000814_06"
"20000814_00"
"20000813_18"
"20000813_12"
"20000813_06"
"20000813_00"
"20000812_18"
"20000812_12"
"20000812_06"
"20000812_00"
"20000811_18"
"20000811_12"
"20000811_06"
"20000810_18"
"20000810_12"
"20000810_06"
"20000810_00")

enddates=("20000810_12"
"20000810_06"
"20000810_00"
"20000809_18"
"20000809_12"
"20000809_06"
"20000809_00"
"20000808_18"
"20000808_12"
"20000808_06"
"20000808_00"
"20000807_18"
"20000807_12"
"20000807_06"
"20000807_00"
"20000806_18"
"20000806_12"
"20000806_06"
"20000806_18"
"20000806_12"
"20000806_06"
"20000806_00")
#
# move files to lagranto
mv ${mpath}/${month}/P* $lagranto
mv ${mpath}/${month}/S* $lagranto

# create starting points 
# every grid [12E, 17E, 46N, 48.5N]
# every  hPa

#csh create_startf.sh $startdate startf.2 'box.grid(12,17,46,48.5) @list(850)@hPa' -changet
csh create_startf.sh ${startdate} startf.2 'point(15,47) @list(850)@hPa' -changet

# calculate trajectories. -j flag prevents stops at the ground level
csh caltra.sh $startdate $enddate startf.2 traj.1 -j  

#loop through the startdates

for date in {0..23}; do
    echo $date
    csh create_startf.sh ${startdates[date]} startf.2 'point(15,47) @list(850)@hPa' -changet
    csh caltra.sh ${startdates[date]} ${enddates[date]} startf.2 traj_$date.1 -j 
    csh select.sh traj_$date.1 wcb_$date.1 "IN:LAT:-30,80:ALL"
    csh trace.sh wcb_$date.1 wcb_$date.1

done 
# select the trajectory within certain coordinates
csh select.sh traj.1 wcb.1 "IN:LAT:-30,80:ALL"

# tracing additional variables aside pressure (variables must be defined in tracevars in the runing directory)
csh trace.sh wcb.1 wcb.1

# move files back to the original directory
echo "--------cleaning up ......"
mv P* ${mpath}/${month}
mv S* ${mpath}/${month}

mv wcb.1 traj.1 startf.2 ${mpath}/${month}
mv wcb_* ${mpath}/${month}/Trace
mv traj_* ${mpath}/${month}/Traj

echo "-------DOne-------"
