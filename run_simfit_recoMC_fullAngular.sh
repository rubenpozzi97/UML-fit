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

ibin=${10}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-fit-intBound-integration
export CMSSWDIR=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/CMSSW_10_4_0/src
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE

export WORKDIR=$PWD
cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

echo setting HOME to $HOME 
echo setting CMSSWDIR to $CMSSWDIR

cd $WORKDIR

nbin=0
while read -a line; do
    abin[$nbin]=${line[0]}
    nbin=$((nbin+1))
done < $HOME/../confSF/KDE_SF.list
bin=${abin[$ibin]}

echo 'now submitting for bin ' ${bin}

if [ ! -r $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root not found
    exit 1
fi
if [ ! -r $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root not found
    exit 1
fi
if [ ! -r $HOME/simfit_recoMC_fullAngular_parSub ]; then
    echo $HOME/simfit_recoMC_fullAngular_parSub not found
    exit 1
fi

cp $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root .
cp $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root .
cp $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root .
cp $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root .
cp $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root .
cp $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root .
cp $HOME/simfit_recoMC_fullAngular_parSub .

echo ./simfit_recoMC_fullAngular_parSub ${bin} ${par} ${f1} ${f4} ${f5} ${b1} ${b4} ${b5} ${m1} ${m4} ${m5} ${multi} ${nsam} ${plot} ${save} 2016 2017 2018
./simfit_recoMC_fullAngular_parSub ${bin} ${par} ${f1} ${f4} ${f5} ${b1} ${b4} ${b5} ${m1} ${m4} ${m5} ${multi} ${nsam} ${plot} ${save} 2016 2017 2018

cp recoBoundDist* $HOME/plotSimFit_d/

rm simfit_recoMC_fullAngular_parSub
rm recoMCDataset_b*
rm KDEeff_b*
rm recoBoundDist*
rm simFitResult*
