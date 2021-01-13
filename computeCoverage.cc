#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

using namespace std;

string label = "toy-b2tuned";
bool isExtended = true;

static const int nPars = 8;
string parName [nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};
string parTitle[nPars] = {"F_{L}","P_{1}","P_{2}","P_{3}","P'_{4}","P'_{5}","P'_{6}","P'_{8}"};

void computeCoverage (int Index, int parity=1)
{

  int binIndex = -1;
  int toyIndex = Index;
  if (isExtended) {
    binIndex = Index;
    toyIndex = -1;
  }    

  vector< vector<double> > vCover (nPars);
  vector< vector<double> > vCoErr (nPars);

  vector< vector<double> > vX (nPars);
  vector< vector<double> > vXe (nPars);

  vector< string > binLab (0);

  fstream fs ("../confSF/toy_coord.list", fstream::in);
  string fullLine = "";

  uint iLabel = 0;
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
    if (binIndex>=0 && binIndex!=q2Bin) continue;
    vector<double> vRef (nPars);
    for (int iPar=0; iPar<nPars; ++iPar)
      vRef[iPar] = atof(splitLine[iPar+1].c_str());

    string binLabel = Form("Bin %i",q2Bin);
    if (isExtended) binLabel = parTitle[iLabel/2] + (iLabel%2==0?" high":" low");
    binLab.push_back(binLabel);
    ++iLabel;

    TChain MINOS_output("MINOS_output","");
    string confString = Form("b%i_",q2Bin);
    for (int iPar=0; iPar<nPars; ++iPar)
      confString = confString + Form((iPar>0?"-%.3f":"%.3f"),vRef[iPar]);
    // string filename = Form("toyFitResults_b%i/simFitResult_",q2Bin)+label+"_fullAngular_201620172018_"+confString+"_s*.root";
    string filename = Form("toyFitResults4_b%i/simFitResult_toy2_fullAngular_201620172018_",q2Bin)+confString+"_s*.root";
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
  string plotTitle = "Toy coverage study;;coverage";
  if (isExtended) plotTitle = Form("Toy coverage study - q2 bin %i;;coverage",binIndex);
  auto hLab = new TH1S ("hLab",plotTitle.c_str(),binLab.size(),0,binLab.size());
  for (int iBin=0; iBin<binLab.size(); ++iBin)
    hLab->GetXaxis()->SetBinLabel(iBin+1,binLab[iBin].c_str());
  hLab->SetMinimum(58);
  hLab->SetMaximum(88);
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
