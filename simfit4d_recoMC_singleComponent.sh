#!/bin/bash

par=1
tag=1
year1=2016
year2=2017
year3=0
datalike=1

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit4d ]; then mkdir logs_simFit4d; fi
if [ ! -d simFit4dResults ]; then mkdir simFit4dResults; fi
if [ ! -d plotSimFit4d_d ]; then mkdir plotSimFit4d_d; fi

# Compile dictionary and macro
make AngDict
make RooDoubleCBDict
make simfit4d_recoMC_singleComponent

while [ ${tag} -ge 0 ]; do
    while read -a line; do
	bin=${line[0]}
	
	./simfit4d_recoMC_singleComponent ${bin} ${par} ${tag} 1 1 ${datalike} ${year1} ${year2} ${year3}\
	    &> logs_simFit4d/simfit4d_recoMC_singleComponent_${bin}_${par}_${tag}_${year1}_${year2}_${year3}_dataStat${datalike}.out &

    done < ../confSF/KDE_SF.list
    tag=$(( ${tag} - 1 ))
done
