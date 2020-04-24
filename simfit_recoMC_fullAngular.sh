#!/bin/bash

par=1
datalike=${1}

pow=${2}

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
make simfit_recoMC_fullAngular

while read -a line; do
	bin=${line[0]}
	
	for year in {2016..2018}; do

	    ./simfit_recoMC_fullAngular ${bin} ${par} ${pow} 0 1 ${datalike} ${year} \
		&> logs_simFit/simfit_recoMC_fullAngular_ad-10-3-e8_${bin}_${par}_${pow}_${datalike}_${year}.out &

	done

	./simfit_recoMC_fullAngular ${bin} ${par} ${pow} 0 1 ${datalike} 2016 2017 2018 \
	    &> logs_simFit/simfit_recoMC_fullAngular_ad-10-3-e8_${bin}_${par}_${pow}_${datalike}_2016_2017_2018.out &

done < ../confSF/KDE_SF.list
