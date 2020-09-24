#!/bin/bash

seed=${1}
nsam=20

save=1

f1=100
f4=10

b1=50
b4=50

m1=0
m4=0

ibin=${2}

export HOME=/afs/cern.ch/work/a/aboletti/private/Kstmumu-Run2/UML-custMinos
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
    agen1[$nbin]=${line[1]}
    agen2[$nbin]=${line[2]}
    agen3[$nbin]=${line[3]}
    agen4[$nbin]=${line[4]}
    agen5[$nbin]=${line[5]}
    agen6[$nbin]=${line[6]}
    agen7[$nbin]=${line[7]}
    agen8[$nbin]=${line[8]}
    nbin=$((nbin+1))
done < $HOME/../confSF/toy_coord.list
bin=${abin[$ibin]}
gen1=${agen1[$ibin]}
gen2=${agen2[$ibin]}
gen3=${agen3[$ibin]}
gen4=${agen4[$ibin]}
gen5=${agen5[$ibin]}
gen6=${agen6[$ibin]}
gen7=${agen7[$ibin]}
gen8=${agen8[$ibin]}

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
if [ ! -r $HOME/simfit_toy_fullAngular ]; then
    echo $HOME/simfit_toy_fullAngular not found
    exit 1
fi

cp $SAMPLEDIR/2016/lmnr/recoMCDataset_b${bin}_2016.root .
cp $SAMPLEDIR/2017/lmnr/recoMCDataset_b${bin}_2017.root .
cp $SAMPLEDIR/2018/lmnr/recoMCDataset_b${bin}_2018.root .
cp $SAMPLEDIR/2016/lmnr/KDEeff_b${bin}_od_2016.root .
cp $SAMPLEDIR/2017/lmnr/KDEeff_b${bin}_od_2017.root .
cp $SAMPLEDIR/2018/lmnr/KDEeff_b${bin}_od_2018.root .
cp $HOME/simfit_toy_fullAngular .

mkdir toyFitResults

echo 
./simfit_toy_fullAngular ${bin} ${gen1} ${gen2} ${gen3} ${gen4} ${gen5} ${gen6} ${gen7} ${gen8} ${f1} ${f4} ${b1} ${b4} ${m1} ${m4} ${seed} ${nsam} 1 ${save} 2016 2017 2018

cp toyFitResults/* $HOME/toyFitResults/

rm -rf toyFitResults

rm simfit_toy_fullAngular
rm recoMCDataset_b*
rm KDEeff_b*
