#!/bin/bash

par=1

multi=1
nsam=100

plot=0
save=0

f1=${1}
f4=${2}
f5=${3}

b1=${4}
b4=${5}
b5=${6}

m1=${7}
m4=${8}
m5=${9}

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

	./simfit_recoMC_fullAngular ${bin} ${par} ${f1} ${f4} ${f5} ${b1} ${b4} ${b5} ${m1} ${m4} ${m5} ${multi} ${nsam} ${plot} ${save} 2016 2017 2018 \
	    &>logs_simFit/simfit_recoMC_fullAngular_${bin}_${par}_${f1}_${f4}_${f5}_${b1}_${b4}_${b5}_${m1}_${m4}_${m5}_${multi}_${nsam}_2016_2017_2018.out &

done < ../confSF/KDE_SF.list
