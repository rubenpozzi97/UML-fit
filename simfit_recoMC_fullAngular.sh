#!/bin/bash

par=1

multi=0
nsam=${1}

plot=0
save=1

f1=100
f4=10

b1=50
b4=50

m1=0
m4=0

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullAngular; then
 
    while read -a line; do
	bin=${line[0]}
	
	# for year in {2016..2018}; do
	
	#     ./simfit_recoMC_fullAngular ${bin} ${par} ${pow} ${multi} ${nsam} ${plot} ${save} ${year} \
	# 	&>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${pow}_${multi}_${nsam}_${year}.out &
	
	# done

	./simfit_recoMC_fullAngular ${bin} ${par} ${f1} ${f4} ${b1} ${b4} ${m1} ${m4} ${multi} ${nsam} 0 ${plot} ${save} 2016 2017 2018 \
	    &>logs_simFit/simfit_recoMC_fullAngular_randLik_${bin}_${par}_${f1}_${f4}_${b1}_${b4}_${m1}_${m4}_${multi}_${nsam}_2016_2017_2018.out &

    done < ../confSF/KDE_SF.list

fi
