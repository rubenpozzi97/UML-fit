#!/bin/bash

par=2 # par = 0 (even efficiency) , par = 1 (odd efficiency) , par > 1 (both parity datasets)
multi=0
nsam=0
constrain=1 # constrain = 0 (unconstrained) , constrain = 1 (cosntrained), constrain = 2 (fixed to MC)
comp=2 # comp = 0 (CT) , comp = 1 (WT) , comp > 1 (both)
dat=1 # dat = 0 (MC) , dat = 1 (data)
pdf_model=0 # pdf_model = 0 (nominal fit)

# if constrain=1, then comp=2 (constrained fit)
# else if constrain=0, then comp=0 or comp=1 (unconstrained fit)
# else if constrain=2, then dat=1 (fit with MC fixed and scale factor)
# if dat = 1, then comp = 2

plot=1
save=1

# Create directories for fit logs, results and plots
if [ ! -d logs_simFitMass ]; then mkdir logs_simFitMass; fi
if [ ! -d simFitMassResults ]; then mkdir -p simFitMassResults; fi
if [ ! -d simFitPullResults_dataStat ]; then mkdir -p simFitPullResults_dataStat; fi
if [ ! -d simFitPullResults_DATA ]; then mkdir -p simFitPullResults_DATA; fi
if [ ! -d simFitPullResults_SF ]; then mkdir -p simFitPullResults_SF; fi
if [ ! -d plotSimMassFit ]; then mkdir plotSimMassFit; fi
if [ ! -d plotSimMassFit_constrain ]; then mkdir plotSimMassFit_constrain; fi
if [ ! -d plotSimMassFit_odd+even ]; then mkdir plotSimMassFit_odd+even; fi
if [ ! -d plotSimMassFit_odd+even_constrain ]; then mkdir plotSimMassFit_odd+even_constrain; fi
if [ ! -d plotSimMassFit_odd+even_factor ]; then mkdir plotSimMassFit_odd+even_factor; fi
if [ ! -d plotSimMassFit_dataStat ]; then mkdir plotSimMassFit_dataStat; fi
if [ ! -d plotSimMassFit_dataStat_multiSample ]; then mkdir plotSimMassFit_dataStat_multiSample; fi
if [ ! -d plotSimMassFit_pulls_dataStat ]; then mkdir plotSimMassFit_pulls_dataStat; fi
if [ ! -d plotSimMassFit_pulls_DATA ]; then mkdir plotSimMassFit_pulls_DATA; fi
if [ ! -d plotSimMassFit_pulls_SF ]; then mkdir plotSimMassFit_pulls_SF; fi
if [ ! -d plotSimMassFit_DATA ]; then mkdir plotSimMassFit_DATA; fi


# Compile dictionary and macro
# make AngDict
if make simfit_recoMC_fullMass; then
    
     bin=4
     for year in {2016,2017,2018}; do

     ./simfit_recoMC_fullMass ${bin} ${par} ${multi} ${nsam} ${plot} ${save} ${constrain} ${comp} ${dat} ${pdf_model} ${year}\
         &>logs_simFitMass/simfit_recoMC_fullMass_${bin}_${par}_${multi}_${nsam}_${constrain}_${comp}_${dat}_${pdf_model}_${year}.out &
    
     done

#    while read -a line; do
#        bin=${line[0]}

#         for year in {2016..2018}; do

#            ./simfit_recoMC_fullMass ${bin} ${par} ${multi} ${nsam} ${plot} ${save} ${constrain} ${comp} ${dat} ${pdf_model} ${year}\
#              &>logs_simFitMass/simfit_recoMC_fullMass_${bin}_${par}_${multi}_${nsam}_${constrain}_${comp}_${dat}_${pdf_model}_${year}.out &

#         done

        #./simfit_recoMC_fullMass ${bin} ${par} ${multi} ${nsam} ${plot} ${save} ${constrain} ${comp} ${dat} ${pdf_model} 2016 2017 2018\
         #   &>logs_simFitMass/simfit_recoMC_fullMass_${bin}_${par}_${multi}_${nsam}_${constrain}_${comp}_${dat}_${pdf_model}_2016_2017_2018.out &

#    done < ../confSF/KDE_SF.list

fi
