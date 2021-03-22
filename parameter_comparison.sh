#!/bin/bash

par=2
nsam=${0}
constrain=0 # constrain = 0 (unconstrained) , constrain = 1 (cosntrained)
comp=0 # comp = 0 (CT) , comp = 1 (WT) , comp > 1 (both)

# Create directories for fit logs, results and plots
if [ ! -d logs_params ]; then mkdir logs_params; fi
if [ ! -d params_results ]; then mkdir -p params_results; fi
if [ ! -d plot_params ]; then mkdir plot_params; fi

# Compile dictionary and macro
# make AngDict
if make parameter_comparison; then

     bin=5
     year=2018

     ./parameter_comparison ${bin} ${par} ${nsam} ${constrain} ${comp} ${year}\
         &>logs_params/parameter_comparison_${bin}_${par}_${nsam}_${constrain}_${comp}_${year}.out &

#    while read -a line; do
#        bin=${line[0]}

#         for year in {2016..2018}; do

#             ./parameter_comparison ${bin} ${par} ${nsam} ${constrain} ${comp} ${year}\
#              &>logs_params/parameter_comparison_${bin}_${par}_${nsam}_${constrain}_${comp}_${year}.out &

#         done

#    done < ../confSF/KDE_SF.list

fi
    

