#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <RooFitResult.h>


void branching_fraction(int year){

  if(year < 2016 || year > 2018){return;}

  const int n_q2Bin = 8;
  double q2Bin[n_q2Bin+1] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};

  //bin width
  double bin_width[n_q2Bin];
  for(int i = 0; i < n_q2Bin; i++){
    bin_width[i] = q2Bin[i+1]-q2Bin[i];
  }

  //efficiencies
  TFile* f_eff = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",year));
  TH1F* eff = (TH1F*)f_eff->Get("eff_x_acc");

  double eff_x_acc[n_q2Bin];
  for(int i = 0; i < n_q2Bin; i++){
    eff_x_acc[i] = eff->GetBinContent(i);
  }

  //yields (from fit to data with constraints)
  double yields[n_q2Bin];
  for(int i = 0; i < n_q2Bin; i++){
    TFile* f_yield = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass%i_DATA_b%ip2c1_subs0CT+WT.root",year,i));
    RooFitResult* fitresult = (RooFitResult*)f_yield->Get(Form("simFitResult_b%ip2c1subs0",i));
    RooRealVar* sig_yield = (RooRealVar*)fitresult->floatParsFinal().find(Form("sig_yield^{%i}",year));
   
    yields[i] = sig_yield->getVal();

    delete f_yield;
  }

  //normalisation channel branching fraction (B0->K*0 J/Psi)
  //const double B_norm = 0.00127;
  const double B_norm = (0.132/100)*(5.96/100);

  //branching fraction
  double branch[n_q2Bin];
  for(int i = 0; i < n_q2Bin; i++){
    if( (i != 4) && (i != 6) ){branch[i] = (yields[i]/yields[4])*(eff_x_acc[4]/eff_x_acc[i])*(B_norm/bin_width[4]);}
  
  //cout << "bin = " << i << " branch = " << branch[i] << endl;
  }

  //branching fraction (resonant channels)
  double branch_resonant = (yields[6]/yields[4])*(eff_x_acc[4]/eff_x_acc[6])*(B_norm/bin_width[4]);
  cout << "year " << year << " - " << "dB_Psi/dq2 = " << branch_resonant << endl;

  return;
}
