#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

using namespace std;

static const int nPars = 8;
string parName [nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

void plotPulls (int q2Bin, int parity=1)
{

  vector<double> vRef (nPars);
  fstream fs ("../confSF/toy_coord.list", fstream::in);
  vector<string> splitLine;
  string fullLine = "";
  getline(fs,fullLine);
  while (!fullLine.empty()) {
    stringstream ss(fullLine);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> vstrings(begin, end);
    if (vstrings.size()>0 && vstrings[0]==Form("%i",q2Bin)) {
      splitLine = vstrings;
      break;
    }
    getline(fs,fullLine);
  }
  if (splitLine.empty()) {
    cout<<"q2 bin "<<q2Bin<<" not found in cofig file"<<endl;
    return;
  }
  if (splitLine.size()!=1+nPars) {
    cout<<"Number of parameters not correct: "<<nPars<<" expected, "<<splitLine.size()-1<<" read"<<endl;
    return;
  }
  for (int iPar=0; iPar<nPars; ++iPar) {
    vRef[iPar] = atof(splitLine[iPar+1].c_str());
    // cout<<parName[iPar]<<":\t"<<vRef[iPar]<<endl;
  }

  TChain MINOS_output("MINOS_output","");
  // string filename = Form("simFitResults/simFitResults_2_6/simFitResult_recoMC_fullAngular201620172018_dataStat_b%i.root",q2Bin);
  // string filename = Form("simFitResults/simFitResult_recoMC_fullAngular201620172018_dataStat_b%i_*.root",q2Bin);
  string confString = Form("b%i_",q2Bin);
  for (int iPar=0; iPar<nPars; ++iPar)
    confString = confString + Form((iPar>0?"-%.3f":"%.3f"),vRef[iPar]);
  string filename = "toyFitResults/simFitResult_toy_fullAngular_201620172018_"+confString+"_s*.root";
  MINOS_output.Add(filename.c_str());

  cout<<MINOS_output.GetEntries()<<endl;

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

  int covLow  = vPull[0].FindBin(-0.999);
  int covHigh = vPull[0].FindBin( 0.999);
  cout<<"Coverage values:"<<endl;

  for (int iPar=0; iPar<nPars; ++iPar) {

    cPulls.cd(iPar+1);
    vPull[iPar].Draw();

    double inside = vPull[iPar].Integral(covLow,covHigh);
    double all = vPull[iPar].Integral();
    double cover = inside / all;
    double coverErr = sqrt( inside * (all-inside) / all / all / all );
    printf("%s:\t%.1f +/- %.1f %%\n",parName[iPar].c_str(),100*cover,100*coverErr);

  }

  cPulls.SaveAs(Form("plotSimFit_d/pullDistr_%s.pdf",confString.c_str()));

}
