#!/bin/bash

set -e

par=1

multi=0
nsam=0

plot=0
save=0

f1=${1}
f4=${2}

b1=${3}
b4=${4}

m1=${5}
m4=${6}

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d simFitResults ]; then mkdir simFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
make simfit_recoMC_fullAngular

while read -a line; do
	bin=${line[0]}
	
	# for year in {2016..2018}; do

	#     ./simfit_recoMC_fullAngular ${bin} ${par} ${pow} ${multi} ${nsam} ${plot} ${save} ${year} \
	# 	&>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${pow}_${multi}_${nsam}_${year}.out &

	# done

	./simfit_recoMC_fullAngular ${bin} ${par} ${f1} ${f4} ${b1} ${b4} ${m1} ${m4} ${multi} ${nsam} ${plot} ${save} 2016 2017 2018 \
	    &>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${f1}_${f4}_${b1}_${b4}_${m1}_${m4}_${multi}_${nsam}_2016_2017_2018.out &

done < ../confSF/KDE_SF.list
