#!/bin/bash

par=2
multi=0
nsam=${0}
constrain=0 # constrain = 0 (unconstrained) , constrain = 1 (cosntrained)
comp=0 # comp = 0 (CT) , comp = 1 (WT) , comp > 1 (both)

plot=1
save=1

# Create directories for fit logs, results and plots
if [ ! -d logs_simFitMass ]; then mkdir logs_simFitMass; fi
if [ ! -d simFitMassResults ]; then mkdir -p simFitMassResults; fi
if [ ! -d plotSimMassFit ]; then mkdir plotSimMassFit; fi
if [ ! -d plotSimMassFit_odd+even ]; then mkdir plotSimMassFit_odd+even; fi

# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullMass; then
    
#     bin=5
#     year=2018

#     ./simfit_recoMC_fullMass ${bin} ${par} ${multi} ${nsam} ${plot} ${save} ${constrain} ${comp} ${year}\
#         &>logs_simFitMass/simfit_recoMC_fullMass_${bin}_${par}_${multi}_${nsam}_${constrain}_${comp}_${year}.out &
    
    while read -a line; do
        bin=${line[0]}

         for year in {2016..2018}; do

             ./simfit_recoMC_fullMass ${bin} ${par} ${multi} ${nsam} ${plot} ${save} ${constrain} ${comp} ${year}\
              &>logs_simFitMass/simfit_recoMC_fullMass_${bin}_${par}_${multi}_${nsam}_${constrain}_${comp}_${year}.out &

         done

        #./simfit_recoMC_fullMass ${bin} ${par} ${multi} ${nsam} ${plot} ${save} ${constrain} ${comp} 2016 2017 2018\
         #   &>logs_simFitMass/simfit_recoMC_fullMass_${bin}_${par}_${multi}_${nsam}_${constrain}_${comp}_2016_2017_2018.out &

    done < ../confSF/KDE_SF.list

fi
