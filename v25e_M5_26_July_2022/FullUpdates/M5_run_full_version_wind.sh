#!/bin/bash

# Laten we samen nog even de settings in de yaml file doorlopen voordat we de full update aanzetten? Jouke
yamlfile=/home/besselaa/Jouke/WindGrid/scripts/FullUpdates/eobs_fg_no_homog.yaml
progdir=/home/besselaa/Jouke/WindGrid/scripts/FullUpdates
#datadir=/data2/Else/Jouke/WindInput/Stations/
# datadir should be the same as stn_dir in eobs_fg_no_homog.yaml!!!!!

version=25.0e

yearstart=1980
#yearend=1989
yearend=2021
month="1;12"

# Kind of function to be able to run a group of steps in parallel
runit() {
 local years=$1
 local part=$2

 echo $years" start: "
 date

 # 3. preprocess stations
 Rscript ${progdir}/M5_preprocess.R yaml=${yamlfile} year="${years}" month=${month} &> Progress/output_preprocess_${part}.txt

 # 4. monthly background grids
 Rscript ${progdir}/M5_forward_sel.R yaml=${yamlfile} year=${years} month=${month} &> Progress/output_forward_sel_${part}.txt

 # 5: gridding anomalies (output in R format)
 Rscript ${progdir}/M5_gprAnomaly.R yaml=${yamlfile} year=${years} month=${month}  &> Progress/output_gprAnomaly_${part}.txt

 # 6: create netcdf-files
 Rscript ${progdir}/M5_postprocess.R yaml=${yamlfile} year=${years} month=${month} &> Progress/output_postprocess_${part}.txt
 
 echo ${years}" done:"
 date
   
}

 

for start in $(seq ${yearstart} 10 ${yearend}); do
    end=$(($start + 9))
#    end=$(($start + 1 ))
    part=${start}-${end}
    period="${start};${end}"
    
# echo $period

    # Does all gridding steps in parallel for groups of years
    runit "$period" "$part" &
done

wait

exit
