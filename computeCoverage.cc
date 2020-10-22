#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

using namespace std;

string label = "toy2";

static const int nPars = 8;
string parName [nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

void computeCoverage (int toyIndex, int parity=1)
{

  vector< vector<double> > vCover (nPars);
  vector< vector<double> > vCoErr (nPars);

  vector< vector<double> > vX (nPars);
  vector< vector<double> > vXe (nPars);

  vector< char* > binLab (0);

  fstream fs ("../confSF/toy_coord.list", fstream::in);
  string fullLine = "";

  for (int iToy=0; (toyIndex<0 || iToy<=toyIndex); ++iToy) {
    getline(fs,fullLine);
    if (fullLine.empty()) break;
    if (iToy<toyIndex) continue;
    cout<<"Running "<<fullLine<<endl;

    stringstream ss(fullLine);
    istream_iterator<string> begin(ss);
    istream_iterator<string> end;
    vector<string> splitLine(begin, end);

    if (splitLine.size()!=1+nPars) {
      cout<<"Number of parameters not correct: "<<nPars<<" expected, "<<splitLine.size()-1<<" read"<<endl;
      return;
    }

    int q2Bin = atoi(splitLine[0].c_str());
    vector<double> vRef (nPars);
    for (int iPar=0; iPar<nPars; ++iPar)
      vRef[iPar] = atof(splitLine[iPar+1].c_str());

    binLab.push_back(Form("bin %i",q2Bin));

    TChain MINOS_output("MINOS_output","");
    string confString = Form("b%i_",q2Bin);
    for (int iPar=0; iPar<nPars; ++iPar)
      confString = confString + Form((iPar>0?"-%.3f":"%.3f"),vRef[iPar]);
    string filename = Form("toyFitResults_b%i/simFitResult_",q2Bin)+label+"_fullAngular_201620172018_"+confString+"_s*.root";
    MINOS_output.Add(filename.c_str());

    int ntoys = MINOS_output.GetEntries();
    cout<<"Toys: "<<ntoys<<endl;

    vector<double> vHigh(nPars);
    vector<double> vLow (nPars);
    vector<int> nIn (nPars);
    for (int iPar=0; iPar<nPars; ++iPar) {
      
      MINOS_output.SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh[iPar]);
      MINOS_output.SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow [iPar]);
      
      nIn[iPar] = 0;

    }

    for (int iEn=0; iEn<MINOS_output.GetEntries(); ++iEn) {
      MINOS_output.GetEntry(iEn);
      for (int iPar=0; iPar<nPars; ++iPar)
	if ( vRef[iPar] < vHigh[iPar] && vRef[iPar] > vLow[iPar] )
	  ++nIn[iPar];
    }

    for (int iPar=0; iPar<nPars; ++iPar) {

      vCover[iPar].push_back( 100.0 * nIn[iPar] / ntoys );
      vCoErr[iPar].push_back( 100.0 * sqrt( 1.0* nIn[iPar] * (ntoys-nIn[iPar]) / ntoys / ntoys / ntoys ) );
      printf("%.1f ",vCover[iPar].back());

      vX[iPar].push_back( vX[iPar].size() + ( (2*iPar+1) / (2.0*nPars) ) );
      vXe[iPar].push_back( 0 );

    }

    cout<<endl;

  }

  if (vX[0].empty()) {
    cout<<"No corresponding toys found"<<endl;
    return;
  }

  TCanvas cCover("cCover","Coverage values",2000,1000);

  vector<TGraphErrors*> grCover (nPars);
  for (int iPar=0; iPar<nPars; ++iPar) {
    grCover[iPar] = new TGraphErrors(vX[iPar].size(),&vX[iPar][0],&vCover[iPar][0],&vXe[iPar][0],&vCoErr[iPar][0]);
    grCover[iPar]->SetName(Form("grCov%s",parName[iPar].c_str()));
    grCover[iPar]->SetMarkerStyle(20);
  }

  gStyle->SetOptStat(0);
  cCover.cd()->SetTicky(2);
  auto hLab = new TH1S ("hLab","Toy coverage study;;coverage",binLab.size(),0,binLab.size());
  for (int iBin=0; iBin<binLab.size(); ++iBin)
    hLab->GetXaxis()->SetBinLabel(iBin+1,binLab[iBin]);
  hLab->SetMinimum(63);
  hLab->SetMaximum(78);
  hLab->GetYaxis()->SetTickLength(0.006);
  hLab->Draw();

  TLine line (0,68.27,binLab.size(),68.27);
  line.SetLineStyle(2);
  line.SetLineColor(15);
  line.Draw();

  for (int iPar=0; iPar<nPars; ++iPar)
    grCover[iPar]->Draw("P");

  cCover.SaveAs(Form("plotSimFit_d/coverages_%s.pdf",label.c_str()));

}
