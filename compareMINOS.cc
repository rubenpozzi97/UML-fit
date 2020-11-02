#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TLine.h>

using namespace std;

// string label1 = "6";
// string label2 = "7";

static const int nPars = 8;
string parName [nPars] = {"Fl","P1","P2","P3","P4p","P5p","P6p","P8p"};

int nSeeds = 50;

bool verb = false;

void compareMINOS (string label1, string label2, int toyIndex=0, int parity=1)
{

  fstream fs ("../confSF/toy_coord.list", fstream::in);
  string fullLine = "";

  char h1[7];
  char h2[7];
  char l1[7];
  char l2[7];

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

    string confString = Form("b%i_",q2Bin);
    for (int iPar=0; iPar<nPars; ++iPar)
      confString = confString + Form((iPar>0?"-%.3f":"%.3f"),vRef[iPar]);

    int nErrNToy = 0;
    int nEnt = 0;
    int nErrDiff = 0;
    int nChecked = 0;
    vector<int> nDiff (0);

    for (int iSeed=0; iSeed<nSeeds; ++iSeed) {
      string filename1 = Form("toyFitResults_b%i_v%s/simFitResult_toy%s_fullAngular_201620172018_%s_s%i.root",q2Bin,label1.c_str(),label1.c_str(),confString.c_str(),iSeed+1);
      string filename2 = Form("toyFitResults_b%i_v%s/simFitResult_toy%s_fullAngular_201620172018_%s_s%i.root",q2Bin,label2.c_str(),label2.c_str(),confString.c_str(),iSeed+1);

      TChain MINOS_output1("MINOS_output","");
      TChain MINOS_output2("MINOS_output","");
      MINOS_output1.Add(filename1.c_str());
      MINOS_output2.Add(filename2.c_str());

      int ntoys1 = MINOS_output1.GetEntries();
      int ntoys2 = MINOS_output2.GetEntries();
      if ( ntoys1==0 || ntoys2==0 ) continue;

      if ( ntoys1 != ntoys2 ) {
	++nErrNToy;
	if (verb) cout<<iToy<<"\t"<<iSeed+1<<"\tdifferent nr. of converging toys"<<endl;
	continue;
      }

      if (verb) cout<<iToy<<"\t"<<iSeed+1<<"\tcomparing "<<ntoys1<<" results"<<endl;

      vector<double> vBest1(nPars);
      vector<double> vHigh1(nPars);
      vector<double> vLow1 (nPars);
      vector<double> vBest2(nPars);
      vector<double> vHigh2(nPars);
      vector<double> vLow2 (nPars);

      for (int iPar=0; iPar<nPars; ++iPar) {
      
	MINOS_output1.SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vBest1[iPar]);
	MINOS_output1.SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh1[iPar]);
	MINOS_output1.SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow1 [iPar]);
      
	MINOS_output2.SetBranchAddress(Form("%s_best",parName[iPar].c_str()),&vBest2[iPar]);
	MINOS_output2.SetBranchAddress(Form("%s_high",parName[iPar].c_str()),&vHigh2[iPar]);
	MINOS_output2.SetBranchAddress(Form("%s_low" ,parName[iPar].c_str()),&vLow2 [iPar]);
      
      }

      for (int iEn=0; iEn<ntoys1; ++iEn) {
	MINOS_output1.GetEntry(iEn);
	MINOS_output2.GetEntry(iEn);
	bool differentResult = false;
	++nEnt;
	for (int iPar=0; iPar<nPars; ++iPar)
	  if ( fabs(vBest1[iPar]-vBest2[iPar])>0.0005 ) {
	    differentResult = true;
	    ++nErrDiff;
	    if (verb) cout<<iToy<<"\t"<<iSeed+1<<"\t"<<iEn<<"\tdifferent results"<<endl;
	    break;
	  }
	if (differentResult) continue;
	for (int iPar=0; iPar<nPars; ++iPar) {
	  sprintf(h1,"%.3f",vHigh1[iPar]);
	  sprintf(h2,"%.3f",vHigh2[iPar]);
	  sprintf(l1,"%.3f",vLow1[iPar]);
	  sprintf(l2,"%.3f",vLow2[iPar]);
	  // if ( fabs(vHigh1[iPar]-vHigh2[iPar])>0.00049 || 
	  //      fabs(vLow1 [iPar]-vLow2 [iPar])>0.00049 )
	  // if ( ((int)(1000*vHigh1[iPar])) != ((int)(1000*vHigh2[iPar])) || 
	  //      ((int)(1000*vLow1 [iPar])) != ((int)(1000*vLow2 [iPar])) ) 
	  ++nChecked;
	  // if ( strcmp(h1,h2) != 0 || strcmp(l1,l2) != 0 ) {
	    if (verb) printf("%i\t%i\t%i\t%i\tdiff: %.3f -> %.3f\t%.3f -> %.3f\n",iToy,iSeed+1,iEn,iPar,
	    	   vHigh1[iPar],vHigh2[iPar],
	    	   vLow1 [iPar],vLow2 [iPar]);
	    double maxDiff = 1e3*TMath::Max(fabs(atof(h1)-atof(h2)),fabs(atof(h1)-atof(h2)));
	    uint hDiff = (uint)round(maxDiff);
	    // if (hDiff==0) continue;
	    while (nDiff.size() < hDiff+1) nDiff.push_back(0);
	    ++nDiff[hDiff];
	  // }

	}

      }

    }

    cout<<"Bad trees "<<nErrNToy<<"/"<<nSeeds<<" bad results "<<nErrDiff<<"/"<<nEnt<<endl;
    cout<<nChecked<<" ->\t";
    for (int i=0; i<nDiff.size(); ++i) cout<<nDiff[i]<<" ("<<i<<")\t";
    cout<<endl;
    for (int i=1; i<nDiff.size(); ++i) printf("%.1f\t",100.0*nDiff[i]/nChecked);
    cout<<endl;
    
  }

}






