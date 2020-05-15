#!/bin/bash

par=1
y6=${1}
y7=${2}
y8=${3}

gen=0

bana=0

# Create directories for fit logs, results and plots
if [ ! -d logs_simFit ]; then mkdir logs_simFit; fi
if [ ! -d plotSimFit_d ]; then mkdir plotSimFit_d; fi

# Compile dictionary and macro
# make AngDict
make plotMultiResults

while read -a line; do
	bin=${line[0]}
	
	./plotMultiResults ${bin} ${par} ${y6} ${y7} ${y8} ${gen} ${bana} \
	    &>logs_simFit/plotMultiResults_${bin}_${par}_${y6}_${y7}_${y8}_${gen}_${bana}.out &

done < ../confSF/KDE_SF.list
