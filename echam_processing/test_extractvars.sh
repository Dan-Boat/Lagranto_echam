# EXTRACTING ADDITIONAL VARS (GEOSP, UV, OMEGA, SLP) FROM OUTPUT RAW FILES (@daniel.boateng@uni-tuebingen.de)
module load afterburner_4.8.0
#
#
#
# USER SPECIFICATION 
YRI=1003
YRF=1017
YR=$YRI

mPATH=/esd/esd02/data/climate_models/echam/echam_output/ESD
EXP=a001
FOLDER_NAME=a001_hpc-bw_e5w2.3_t159_PI_Alps_east_300_t159l31.6h

mkdir ${mPATH}/${FOLDER_NAME}/output_raw/extract_vars
oPATH=${mPATH}/${FOLDER_NAME}/output_raw/extract_vars
processed=${mPATH}/${FOLDER_NAME}/output_processed/6h_MEANS

#STARTING lOOP
while ((${YR} <= ${YRF})); do
        for i in {1..12}; do
                if ((${i} < 10)); then
                       MON=0${i}
                else
                       MON=${i}
                fi
		#extracting of relevant variables 
		# cdo -s -f nc -selname, st ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_temp.01.nc
		# cdo -s -f nc -selname, relhum ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_relhum.01.nc
		# cdo -s -f nc -selname, geosp ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_geosp.01.nc
		# cdo -s -f nc -selname, lsp ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_lsp.01.nc
		# cdo -s -f nc -selname, svo ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_svo.01.nc
		# cdo -s -f nc -selname, aps ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_aps.01.nc
		# cdo -s -f nc -selname, aprl ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_aprl.01.nc
		# cdo -s -f nc -selname, aprc ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_aprc.01.nc
		# cdo -s -f nc -selname, evap ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_evap.01.nc
		# #
		# cdo -s -f nc -chvar,aprl,aprt -chcode,142,260 -add ${oPATH}/${EXP}_${YR}${MON}_aprl.01.nc ${oPATH}/${EXP}_${YR}${MON}_aprc.01.nc ${oPATH}/${EXP}_${YR}${MON}_aprt.01.nc
		# cdo -s -f nc -chvar,aprt,pe -chcode,260,265 -add ${oPATH}/${EXP}_${YR}${MON}_aprt.01.nc ${oPATH}/${EXP}_${YR}${MON}_evap.01.nc ${oPATH}/${EXP}_${YR}${MON}_pe.01.nc
		#
		# geospotential height at pressure levels 
        	after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_geopoth.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=156
&END
EOF
			#log of surface pressure
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_lsp.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=152
&END
EOF
			# temperature
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_temp.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=130
&END
EOF
			# surface pressure 
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_aps.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=134
&END
EOF
			# relative humidiy 
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_relhum.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=157
&END
EOF
			# vorticity 
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_svo.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=138
&END
EOF
			# zonal velocity : u
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_u.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=131
&END
EOF
			# meridional velocity : v
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_v.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=132
&END
EOF
			# vertical velocity: omega
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_omega.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=135
&END
EOF

			#convective precipitation
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_aprc.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=142
&END
EOF
			# large scale precipitation
			after ${mPATH}/${FOLDER_NAME}/output_raw/${EXP}_${YR}${MON}.01.nc ${oPATH}/${EXP}_${YR}${MON}_aprl.01.nc <<EOF
&SELECT
TYPE=20, FORMAT=2,
CODE=143
&END
EOF
        done # done for month loop
        YR=$(($YR+1))
	echo ${YR}
done  #done for year loop 



