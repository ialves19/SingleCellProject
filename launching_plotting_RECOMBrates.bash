#!/bin/bash

#setting the start of the job
res1=$(date +%s.%N)

module load rstats/3.6
#this script should be run as : qsub -P awadallalab -V -cwd -b y -pe smp 4 -l h_vmem=8G ./launching_plotting_RECOMBrates.bash 21
#it takes 1 arg: chrNb 

chmod +x $HOME/recombinationRates_and_plotting.R
Rscript $HOME/recombinationRates_and_plotting.R $1

#timing the job
res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)

echo "Total runtime: $dd days $dh hrs $dm min $ds secs"

