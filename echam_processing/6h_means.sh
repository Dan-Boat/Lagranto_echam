# 
# This script is meant to calculate the 6h means of model output variables extracted for trajectory analysis from ECHAM
# Setting paths and specifications
mpath=/esd/esd02/data/climate_models/echam/echam_output/ESD
folder=a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h
add=("u"
"v"
"aps"
"svo"
"lsp"
"omega"
"temp"
"geopoth"
"relhum"
"aprc"
"aprl")

YRI=1003
YRF=1017
varname=("u"
"v"
"aps"
"svo"
"lsp"
"omega"
"st"
"geopoth"
"relhum"
"aprl"
"aprc")

EXP=a001
#
#
# loading cdo libraries 
module load cdo_1.9.5

ipath=${mpath}/${folder}/output_raw/extract_vars/
outpath=${mpath}/${folder}/output_processed/
echo ${ipath}
echo "----creating path for 6h_means"
mkdir ${outpath}6h_means
echo "-----folder created......"
opath=${outpath}6h_means
#
for var in {1..11}; do
	echo "-----Averaging for variable:........."
	echo ${add[var]}
	echo $opath
	for i in {1..12}; do 
		if ((${i} < 10)); then 
			MON=0${i}
		else
			MON=${i}
		fi 
		echo $MON
		files_in_string=''
		YR=${YRI}
		while ((${YR} <= ${YRF} )); do 
			echo "-------select varname to be processed....."
			cdo select,name=${varname[var]} ${ipath}${EXP}_${YR}${MON}_${add[var]}.01.nc ${opath}/${EXP}_$(($YR))${MON}_${add[var]}.01.nc
			files_in_string=${files_in_string}" "${opath}/${EXP}_$(($YR))${MON}_${add[var]}.01.nc

			YR=$(($YR+1))
		done
	#	
		cdo -P 8 ensmean ${files_in_string} ${opath}/${EXP}_${YRI}_${YRF}_${MON}_${add[var]}.01.nc
		rm -f ${files_in_string}
	done
	#
	cdo mergetime ${opath}/${EXP}_${YRI}_${YRF}_*_${add[var]}.01.nc ${opath}/${EXP}_${YRI}_${YRF}_6h_lterm_${add[var]}.01.nc
done

