#!/bin/bash

par=1

multi=0
nsam=${1}

plot=0
save=1

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults4d ]; then mkdir -p simFitResults4d; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullAngularMass; then
 
    while read -a line; do
	bin=${line[0]}
	
	# for year in {2016..2018}; do
	
	#     ./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 0 ${plot} ${save} ${year} \
	# 	&>logs_simFit/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_${year}.out &
	
	# done

	./simfit_recoMC_fullAngularMass ${bin} ${par} ${multi} ${nsam} 0 ${plot} ${save} 2016 2017 2018 \
	    &>logs_simFit/simfit_recoMC_fullAngularMass_${bin}_${par}_${multi}_${nsam}_2016_2017_2018.out &

    done < ../confSF/KDE_SF.list

fi
