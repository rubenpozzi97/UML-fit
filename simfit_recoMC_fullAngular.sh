#!/bin/bash

par=1
pow=1

multi=${1}
nsam=${2}

plot=0
save=1

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

	    ./simfit_recoMC_fullAngular ${bin} ${par} ${pow} ${multi} ${nsam} ${plot} ${save} ${year} \
		&>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${pow}_${multi}_${nsam}_${year}.out &

	done

	./simfit_recoMC_fullAngular ${bin} ${par} ${pow} ${multi} ${nsam} ${plot} ${save} 2016 2017 2018 \
	    &>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${pow}_${multi}_${nsam}_2016_2017_2018.out &

done < ../confSF/KDE_SF.list
