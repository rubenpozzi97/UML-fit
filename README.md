# UML-fit

# How to run

**Table of contents**  
[Setup working area](#setup)  
[Create datasets](#createDatasets)  
[Fit to signal MC](#fitMC)  

<a name="setup"/>

## Setup working area
Make sure to be in an area with ROOT v6.12 (or later) available. The code could work with previous versions, but it is not guaranteed.
If you are on CERN lxplus machines, you can achieve it by running:
```sh
export SCRAM_ARCH=slc7_amd64_gcc820
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src/ && cmsenv && cd ../..
```
Clone this branch in the working directory:
```sh
git clone -b master git@github.com:CMSKStarMuMu/UML-fit.git
cd UML-fit
```
<a name="createDatasets"/>

## Create datasets
If needed, change the [location of the ntuples](createDataset.cc#L55-L62), which need to be produced with the code in the [B0KstMuMuNtuple repository](https://github.com/CMSKStarMuMu/B0KstMuMuNtuple).  
Then, you can produce the files which will contain the needed datasets, that is, correctly and wrongly tagged reco'ed events from the MC. The datasets will include the following variables: ctK, ctL, phi, mass, rand.
rand is a random variable uniformly generated between 0 and 1, to be used to select a given statistics from the MC sample.  
Example on how to run for 2017 ntuples:
```sh
root -q -b 'createDataset.cc(7)'
```


<a name="performFit"/>

## Perform fits to MC dataset with PDF*eff

The fit is performed by the simfit_recoMC_singleComponent code.
It requires in input:
* a root file containing the workspace with the roodatasets to be fitted (produced in the step above)
* the files containing the KDE modeling of the efficiency and, if possible, the corresponding partial integrals. These are produced by the workflow described in the eff-KDE repository. 

Compile and run with:
```sh
source simfit_recoMC_singleComponent.sh
```
where you have to set the datasets to be considered (set year = 0 to not include the dataset). 
The variable "datalike" sets the statistics to be considered:  
* datalike = 0 means consider the full MC stat (half of it actually)  
* datalike = 1 means consider a data-like statistics.  

The code will produce a root file `simFitResults/fitResult_recoMC_singleComponentXXXX.root` containing the RooFitResult objects, where XXXX describes the considered datasets.
Corresponding fit projection plots are created in `plotSimFit_d/fitResult_recoMC_singleComponent_*.pdf`.

### Plot and compare fit results
```sh
root -b -q 'plotSimFitResults.cc(1)' # for fit with odd efficiency on even dataset and full MC stat
```
while use
```
root -b -q 'plotSimFitResults.cc(1,-1,true, false, false,true)' 
```
to plot the results for the data-like stat fits.
