#!/bin/bash

par=1
datalike=1

pow=${1}
sh1=${2}
sh5=${3}
coef=${4}
int=${5}

# If you omit the integral value, it is automatically read
if [ -z ${int} ]; then
    log=logs_intBound/integrate_boundary_${pow}_${sh1}_${sh5}_${coef}.out
    if [ ! -f ${log} ]; then
	echo "File not found:" ${log}
	return 1
    fi
    arr=($(tail -n1 ${log}));
    int=${arr[6]}
    if [ -z ${int} ]; then
	echo "File not correct:" ${log}
        return 1
    fi
fi

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

	    ./simfit_recoMC_fullAngular ${bin} ${pow} ${sh1} ${sh5} ${coef} ${int} ${year} \
		&> logs_simFit/simfit_recoMC_fullAngular_${bin}_${year}_${pow}_${sh1}_${sh5}_${coef}.out &

	done

	./simfit_recoMC_fullAngular ${bin} ${pow} ${sh1} ${sh5} ${coef} ${int} 2016 2017 2018 \
	    &> logs_simFit/simfit_recoMC_fullAngular_${bin}_2016_2017_2018_${pow}_${sh1}_${sh5}_${coef}.out &

done < ../confSF/KDE_SF.list
