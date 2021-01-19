#!/bin/bash

seed=${1}
nsam=${2}

save=1

f1=100
f4=10

b1=50
b4=50

m1=0
m4=0

# Create directories for fit logs, results and plots
if [ ! -d logs_toyFit ]; then mkdir logs_toyFit; fi
if [ ! -d toyFitResults ]; then mkdir toyFitResults; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
if make simfit_toy_fullAngular; then
 
    while read -a line; do
	bin=${line[0]}
	gen1=${line[1]}
	gen2=${line[2]}
	gen3=${line[3]}
	gen4=${line[4]}
	gen5=${line[5]}
	gen6=${line[6]}
	gen7=${line[7]}
	gen8=${line[8]}
	
	./simfit_toy_fullAngular ${bin} ${gen1} ${gen2} ${gen3} ${gen4} ${gen5} ${gen6} ${gen7} ${gen8} ${f1} ${f4} ${b1} ${b4} ${m1} ${m4} ${seed} ${nsam} 0 ${save} 2016 2017 2018 \
	    &>logs_toyFit/simfit_toy_fullAngular_${bin}_${gen1}_${gen2}_${gen3}_${gen4}_${gen5}_${gen6}_${gen7}_${gen8}_${f1}_${f4}_${b1}_${b4}_${m1}_${m4}_${seed}_${nsam}_2016_2017_2018.out &

    done < ../confSF/toy_coord.list

fi
