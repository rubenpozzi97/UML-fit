#!/bin/bash

# Create directories for fit logs, results and plots
if [ ! -d logs_intBound ]; then mkdir logs_intBound; fi

# Compile dictionary and macro
make integrate_boundary

pow=2

for sh1 in {10000..10000..900}; do
    for sh5 in {20..50..10}; do
	coef=50
	while [ ${coef} -le 100 ]; do

	    echo "Submitting" ${pow} ${sh1} ${sh5} ${coef}

	    ./integrate_boundary ${pow} ${sh1} ${sh5} ${coef} \
	    	&> logs_intBound/integrate_boundary_${pow}_${sh1}_${sh5}_${coef}.out &

	    coef=$((${coef}*2))

	done
    done
done
