#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

#include <RooRealVar.h>
#include <RooFitResult.h>

using namespace std;
using namespace RooFit;

static const int nPars = 8;
string parName [nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

void plotPulls (int q2Bin, int parity=1)
{
  string shortString = Form("b%ip%i",q2Bin,parity);
  /*
  string fResultName = Form("simFitResults/simFitResult_recoMC_fullAngular201620172018_MCStat_b%i.root", q2Bin);
  auto fResult = TFile::Open(fResultName.c_str());
  if (!fResult || fResult->IsZombie()) {
    cout<<"File not found: "<<fResultName<<endl;
    return;
  }
  string fitResName = "simFitResult_"+shortString+"subs0";
  auto fitRes = (RooFitResult*)fResult->Get(fitResName.c_str());
  if (!fitRes || fitRes->IsZombie()) {
    cout<<"Fit result "<<fitResName<<" not found in file "<<fResultName<<endl;
    return;
  }
  */
  TChain MINOS_output("MINOS_output","");
  // string filename = Form("simFitResults/simFitResults_2_6/simFitResult_recoMC_fullAngular201620172018_dataStat_b%i.root",q2Bin);
  string filename = Form("simFitResults/simFitResult_recoMC_fullAngular201620172018_dataStat_b%i_*.root",q2Bin);
  MINOS_output.Add(filename.c_str());

  cout<<MINOS_output.GetEntries()<<endl;

  vector<double> vRef (nPars);
  vector<double> vBest(nPars);
  vector<double> vHigh(nPars);
  vector<double> vLow (nPars);
  vector<TH1D> vPull(nPars);
  for (int iPar=0; iPar<nPars; ++iPar) {

    MINOS_output.SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vBest[iPar]);
    MINOS_output.SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh[iPar]);
    MINOS_output.SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow [iPar]);

    vPull[iPar] = TH1D(Form("hPull%s",parName[iPar].c_str()),
		       Form("%s pull distribution - q2 bin %i",parName[iPar].c_str(),q2Bin),
		       24,-4.0,4.0);

    // vRef[iPar] = ((RooRealVar*)fitRes->floatParsFinal().find(parName[iPar].c_str()))->getValV();
    vRef[iPar] = 0;

  }

  for (int iEn=0; iEn<MINOS_output.GetEntries(); ++iEn) {
    MINOS_output.GetEntry(iEn);
    for (int iPar=0; iPar<nPars; ++iPar)
      if (vBest[iPar]>vRef[iPar])
	vPull[iPar].Fill((vBest[iPar]-vRef[iPar])/(vBest[iPar]-vLow[iPar]));
      else
	vPull[iPar].Fill((vBest[iPar]-vRef[iPar])/(vHigh[iPar]-vBest[iPar]));
  }

  TCanvas cPulls("cPulls",Form("Pull distributions for q2 bin %i",q2Bin),2000,1000);
  cPulls.Divide(4,2);

  for (int iPar=0; iPar<nPars; ++iPar) {

    cPulls.cd(iPar+1);
    vPull[iPar].Draw();

  }

  cPulls.SaveAs(Form("plotSimFit_d/pullDistr_%s.pdf",shortString.c_str()));

}
