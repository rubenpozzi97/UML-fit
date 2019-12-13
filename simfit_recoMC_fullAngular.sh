#!/bin/bash

par=1
year1=2017
year2=0
year3=0
datalike=1

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
make AngDict
make simfit_recoMC_fullAngular

while read -a line; do
	bin=${line[0]}
	
	./simfit_recoMC_fullAngular ${bin} ${par} 1 1 ${datalike} ${year1} ${year2} ${year3}\
	    &> logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${year1}_${year2}_${year3}_dataStat${datalike}.out &

done < ../confSF/KDE_SF.list
