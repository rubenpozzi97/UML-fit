#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <RooDataSet.h>
#include <RooArgList.h>
#include <RooGlobalFunc.h>
#include <map>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooAbsPdf.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooFormulaVar.h>
#include <TH1.h>
#include <TH1D.h>
#include <RooAbsReal.h>
#include <iostream>
#include <TLegend.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <string>
#include <fstream>
#include <TF1.h>
#include <TNtupleD.h>
#include <TLatex.h>
#include <TStyle.h>
#include <vector>
#include <sstream>
#include <iomanip>
#include <RooFit.h>
#include <RooCmdArg.h>
#include <RooPlot.h>
#include <RooArgList.h>

using namespace std;
using namespace RooFit;

void addVariables(RooWorkspace& w);
void createDatasets(bool data, int year, int q2Bin, RooWorkspace& w);
std::vector<TH1D*> sideband_subtraction(RooWorkspace& w, RooWorkspace& w_fit, int n, int n_var, int year, int q2Bin);
TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n, int year, int q2Bin);


void MC_validation(bool create, bool data, int year, int q2Bin){

  RooWorkspace* w = new RooWorkspace("w");
  //addVariables(*w);
  if(create){
    addVariables(*w);
    createDatasets(data, year, q2Bin, *w);
    return;
  }

  // import dataset to workspace
  TString input_file_data = Form("~/public/UML-fit/Datasets/Dataset_DATA_201%i_b%i.root",year,q2Bin);
  TFile* f_input_data = new TFile(input_file_data);
  RooWorkspace* w_input_data = (RooWorkspace*)f_input_data->Get("w");
  RooDataSet* data_DATA = (RooDataSet*)w_input_data->data(Form("data_DATA_201%i_b%i",year,q2Bin));
  w->import(*data_DATA, Rename("data"));
  w->Print();

  int n_bins = 40;

  TString variables[] = {"mumuMass", "mmkMass", "mmpiMass", "pionPt", "kaonPt", "kstarmass", "kstTrkpEta", "kstTrkmEta", "mupEta", "mumEta", "mumuEta", "kstEta", "bEta", "kstTrkpPhi", "kstTrkmPhi", "mupPhi", "mumPhi", "mumuPhi", "kstPhi", "bPhi", "kstTrkpPt", "kstTrkmPt", "mupPt", "mumPt", "mumuPt", "kstPt", "bPt"};

  int n_var = sizeof(variables)/sizeof(variables[0]);

  // Access RooFitResult and fit PDFs produced with simfit_recoMC_fullMass.cc code
  TString input_fit_result = Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c1m0_subs0CT+WT.root",year,q2Bin);
  TFile* f_input_fit_result = new TFile(input_fit_result);
  RooFitResult* fitresult = (RooFitResult*)f_input_fit_result->Get(Form("simFitResult_b%ip2c1m0subs0",q2Bin)); 
  RooWorkspace* w_fit = (RooWorkspace*)f_input_fit_result->Get("ws_pars");
  fitresult->Print("v");

  std::vector<TH1D*> histos_sideband_sub;
  std::vector<TH1D*> histos_splot;
  std::vector<TH1D*> histos_mc;

  // Sideband subtraction
  histos_sideband_sub = sideband_subtraction(*w, *w_fit, n_bins, n_var, year, q2Bin);

  return;
}


std::vector<TH1D*> sideband_subtraction(RooWorkspace& w, RooWorkspace& w_fit, int n, int n_var, int year, int q2Bin){

  RooDataSet* data = (RooDataSet*)w.data("data");
  data->Print();

  //RooAbsPdf* BkModel = w_fit.pdf(Form("CB_bkg_PDF_201%i",year));
  RooRealVar* mass = (w.var("mass"));
  RooRealVar* lamdba = (w_fit.var(Form("lambda^{201%i}",year)));
  RooExponential* BkModel = new RooExponential("BkModel", "BkModel", *mass, *lamdba);

  vector<RooRealVar> variables;

  variables.push_back(*(w.var("mass")));
  variables.push_back(*(w.var("mumuMass")));
  variables.push_back(*(w.var("mmkMass")));
  variables.push_back(*(w.var("mmpiMass")));
  variables.push_back(*(w.var("pionPt")));
  variables.push_back(*(w.var("kaonPt")));
  variables.push_back(*(w.var("kstarmass")));
  variables.push_back(*(w.var("kstTrkpEta")));
  variables.push_back(*(w.var("kstTrkmEta")));
  variables.push_back(*(w.var("mupEta")));
  variables.push_back(*(w.var("mumEta")));
  variables.push_back(*(w.var("mumuEta")));
  variables.push_back(*(w.var("kstEta")));
  variables.push_back(*(w.var("bEta")));
  variables.push_back(*(w.var("kstTrkpPhi")));
  variables.push_back(*(w.var("kstTrkmPhi")));
  variables.push_back(*(w.var("mupPhi")));
  variables.push_back(*(w.var("mumPhi")));
  variables.push_back(*(w.var("mumuPhi")));
  variables.push_back(*(w.var("kstPhi")));
  variables.push_back(*(w.var("bPhi")));
  variables.push_back(*(w.var("kstTrkpPt")));
  variables.push_back(*(w.var("kstTrkmPt")));
  variables.push_back(*(w.var("mupPt")));
  variables.push_back(*(w.var("mumPt")));
  variables.push_back(*(w.var("mumuPt")));
  variables.push_back(*(w.var("kstPt")));
  variables.push_back(*(w.var("bPt")));

  RooDataSet* reduceddata_side;
  RooDataSet* reduceddata_central;

  double left;
  double right;

  if( ((year == 6) &&  (q2Bin <= 3)) ||  ((year == 7) && (q2Bin <=2)) || ((year == 8) && (q2Bin <= 3)) ){left = 5.11;}
  else{left = 5.12;}

  if( ((year == 6) && (q2Bin <= 1)) || ((year == 7) && (q2Bin <= 3)) || ((year == 8) && (q2Bin <= 3)) ){right = 5.45;}
  else{right = 5.44;}

  mass->setRange("right",right,mass->getMax());
  mass->setRange("left",mass->getMin(),left);
  mass->setRange("peak",left,right);
  mass->setRange("peakright",left,mass->getMax());
  mass->setRange("total", mass->getMin(), mass->getMax()); 

  reduceddata_side =  (RooDataSet*)data->reduce(Form("mass>%lf || mass<%lf", right, left));
  reduceddata_central = (RooDataSet*)data->reduce(Form("mass>%lf",left));
  reduceddata_central = (RooDataSet*)reduceddata_central->reduce(Form("mass<%lf",right));

  /*
  TH1D* sidebands = (TH1D*)reduceddata_side->createHistogram("mass1", *mass, Binning(100, 5.03, mass->getMax()));
  TH1D* peak = (TH1D*)reduceddata_central->createHistogram("mass2", *mass, Binning(100, 5.03, mass->getMax()));
  TH1D* total = (TH1D*)data->createHistogram("mass3", *mass, Binning(100, 5.03, mass->getMax()));

  sidebands->SetLineColor(kRed);
  peak->SetLineColor(kBlue);
  total->SetLineColor(kBlack);

  TCanvas c;
  c.cd();
  total->Draw();
  sidebands->Draw("same");
  peak->Draw("same");
  c.SaveAs("~/public/UML-fit/results/datasets.gif");
  */
 
  //Integrating the background distribution 
  RooAbsReal* int_fit_side_right = BkModel->createIntegral(*mass, *mass, "right");
  RooAbsReal* int_fit_side_left = BkModel->createIntegral(*mass, *mass, "left");
  RooAbsReal* int_fit_peak = BkModel->createIntegral(*mass, *mass, "peak");

  cout << "Integral left band = " << int_fit_side_left->getVal() << endl;
  cout << "Integral right band = " << int_fit_side_right->getVal() << endl;
  cout << "Integral peak = " << int_fit_peak->getVal() << endl;
  cout << "normalisation = " << (BkModel->createIntegral(*mass, *mass, "total"))->getVal() << endl;

  double factor = (int_fit_peak->getVal())/(int_fit_side_right->getVal() + int_fit_side_left->getVal());
  std::cout << std::endl << "Factor: " << factor << std::endl;

  std::vector<TH1D*> histos;

  histos.push_back(create_histogram(variables[1], "mumuMass",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[2], "mmkMass",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[3], "mmpiMass",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[4], "pionPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[5], "kaonPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[6], "kstarmass",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[7], "kstTrkpEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[8], "kstTrkmEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[9], "mupEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[10], "mumEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[11], "mumuEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[12], "kstEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[13], "bEta",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[14], "kstTrkpPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[15], "kstTrkmPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[16], "mupPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[17], "mumPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[18], "mumuPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[19], "kstPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[20], "bPhi",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[21], "kstTrkpPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[22], "kstTrkmPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[23], "mupPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[24], "mumPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[25], "mumuPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[26], "kstPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));
  histos.push_back(create_histogram(variables[27], "bPt",factor, reduceddata_side, reduceddata_central, data, n, year, q2Bin));

  return histos;

}

TH1D* create_histogram(RooRealVar var,TString name, double factor, RooDataSet* reduced, RooDataSet* central, RooDataSet* total, int n, int year, int q2Bin){

  std::cout<< "n in create_histogram = "<< n <<std::endl;
  TH1D* dist_side = (TH1D*)reduced->createHistogram("dist_side",var, Binning(n, var.getMin(), var.getMax()));
  dist_side->SetMarkerColor(kRed);
  dist_side->SetLineColor(kRed);
  dist_side->SetNameTitle("dist_side", "");

  TH1D* hist_dist_peak = (TH1D*)central->createHistogram(var.GetName(), var, Binning(n, var.getMin(), var.getMax()));
  TH1D* dist_peak = new TH1D(*hist_dist_peak);
  dist_peak->SetMarkerColor(kBlue);
  dist_peak->SetLineColor(kBlue);
  dist_peak->SetNameTitle("dist_peak", "");

  hist_dist_peak->SetMarkerColor(kBlack);
  hist_dist_peak->SetLineColor(kBlack);
  hist_dist_peak->SetNameTitle("dist_total", "");

  dist_peak->Add(dist_side, -factor);
  dist_side->Scale(factor);

  dist_peak->SetStats(0);
  dist_side->SetStats(0);
  hist_dist_peak->SetStats(0);
  TCanvas c;

  hist_dist_peak->Draw();
  dist_side->Draw("same");
  dist_peak->Draw("same");

  dist_peak->SetXTitle(var.GetName());
  dist_side->SetXTitle(var.GetName());
  hist_dist_peak->SetXTitle(var.GetName());

  hist_dist_peak->GetYaxis()->SetRangeUser(0, 1.3*hist_dist_peak->GetMaximum());

  TLegend *leg = new TLegend (0.7, 0.9, 0.9, 1.0);
  leg->AddEntry("dist_peak", "Signal", "l");
  leg->AddEntry("dist_side", "Background", "l");
  leg->AddEntry("hist_dist_peak", "Total", "l");
  leg->Draw("same");

  std::cout<<"name: "<<var.GetName()<<std::endl;
  std::cout<<"histo name: "<<hist_dist_peak->GetName()<<std::endl;

  c.SaveAs( "~/public/UML-fit/results/sideband_sub/"+name+Form("sideband_sub_201%i_b%i.gif",year,q2Bin) );
  return dist_peak;

}


void createDatasets(bool data, int year, int q2Bin, RooWorkspace& w){

  static const int nBins = 9;
  float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

  double PDGB0Mass = 5.27958;
  double PDGJpsiMass = 3.096916;
  double PDGPsiPrimeMass = 3.686109;
  double PDGKstMass = 0.896;

  bool isJpsi = false;
  bool isPsi  = false;
  bool isLMNR = false;
  
  cout << "q2Bin = " << q2Bin << endl;

  if(q2Bin==4){isJpsi = true;}
  else if(q2Bin==6){isPsi = true;}
  else{isLMNR = true;}

  // Load ntuples
  TFile*f;

  if(data){f = new TFile(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iData_All_finalSelection.root",year,year));}
  else{
    if(isLMNR){f = new TFile(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iMC_LMNR.root", year, year));}
    else if(isJpsi){f = new TFile(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iMC_JPSI.root", year, year));}
    else if(isPsi){f = new TFile(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iMC_PSI.root", year, year));}
  }

  TTree* t = (TTree*)f->Get("ntuple");

  // variables for selections
  double recoDimuMass;
  t->SetBranchAddress( "mumuMass", &recoDimuMass );

  double recoB0Mass;
  t->SetBranchAddress( "tagged_mass", &recoB0Mass );

  bool passB0Psi_lmnr, passB0Psi_jpsi, passB0Psi_psip;
  t->SetBranchAddress( "passB0Psi_lmnr", &passB0Psi_lmnr );
  t->SetBranchAddress( "passB0Psi_jpsi", &passB0Psi_jpsi );
  t->SetBranchAddress( "passB0Psi_psip", &passB0Psi_psip );

  double tagB0;
  t->SetBranchAddress( "tagB0"    , &tagB0     );

  Long64_t eventN;
  t->SetBranchAddress( "eventN", &eventN     );

  double wt_mass, wt_kstarmass, kaonPt, pionPt, mmpiMass, mmkMass;
  t->SetBranchAddress( "wt_mass",      &wt_mass      );
  t->SetBranchAddress( "wt_kstarmass", &wt_kstarmass );
  t->SetBranchAddress( "kaonPt",       &kaonPt       );
  t->SetBranchAddress( "pionPt",       &pionPt       );
  t->SetBranchAddress( "mmpiMass",     &mmpiMass     );
  t->SetBranchAddress( "mmkMass",      &mmkMass      );

  double x0Cut=-0.4;
  double y0Cut= 0.3;
  double x1Cut= 0.6;
  double y1Cut=-0.1;

  double x_0Cut=3;
  double y_0Cut=3.8;
  double x_1Cut=3.6;
  double y_1Cut=4.8;

  double CutX1=3.2;
  double CutX2=3.6;
  double CutY1=4.7;
  double CutY2=4.9;

  double PUweight = 1.;
  //t->SetBranchAddress( "weight", &PUweight );

  // Workspace variables
  double bPt, kstarmass, kstTrkpEta, kstTrkmEta, mupEta, mumEta, mumuEta, kstEta, bEta, kstTrkpPhi, kstTrkmPhi, mupPhi, mumPhi, mumuPhi, kstPhi, bPhi, kstTrkpPt, kstTrkmPt, mupPt, mumPt, mumuPt, kstPt;

  t->SetBranchAddress("bPt", &bPt);
  t->SetBranchAddress("kstarmass", &kstarmass);
  t->SetBranchAddress("kstTrkpEta", &kstTrkpEta);
  t->SetBranchAddress("kstTrkmEta", &kstTrkmEta);
  t->SetBranchAddress("mupEta", &mupEta);
  t->SetBranchAddress("mumEta", &mumEta);
  t->SetBranchAddress("mumuEta", &mumuEta);
  t->SetBranchAddress("kstEta", &kstEta);
  t->SetBranchAddress("bEta", &bEta);
  t->SetBranchAddress("kstTrkpPhi", &kstTrkpPhi);
  t->SetBranchAddress("kstTrkmPhi", &kstTrkmPhi);
  t->SetBranchAddress("mupPhi", &mupPhi);
  t->SetBranchAddress("mumPhi", &mumPhi);
  t->SetBranchAddress("mumuPhi", &mumuPhi);
  t->SetBranchAddress("kstPhi", &kstPhi);
  t->SetBranchAddress("bPhi", &bPhi);
  t->SetBranchAddress("kstTrkpPt", &kstTrkpPt);
  t->SetBranchAddress("kstTrkmPt", &kstTrkmPt);
  t->SetBranchAddress("mupPt", &mupPt);
  t->SetBranchAddress("mumPt", &mumPt);
  t->SetBranchAddress("mumuPt", &mumuPt);
  t->SetBranchAddress("kstPt", &kstPt);

  // loop over tree events
  RooRealVar* Bmass = (RooRealVar*)w.var("mass");
  RooRealVar* BkaonPt = (RooRealVar*)w.var("kaonPt");
  RooRealVar* BpionPt = (RooRealVar*)w.var("pionPt");
  RooRealVar* BmmpiMass = (RooRealVar*)w.var("mmpiMass");
  RooRealVar* BmmkMass = (RooRealVar*)w.var("mmkMass");
  RooRealVar* BmumuMass = (RooRealVar*)w.var("mumuMass");
  RooRealVar* BbPt = (RooRealVar*)w.var("bPt");
  RooRealVar* Bkstarmass = (RooRealVar*)w.var("kstarmass");
  RooRealVar* BkstTrkpEta = (RooRealVar*)w.var("kstTrkpEta");
  RooRealVar* BkstTrkmEta = (RooRealVar*)w.var("kstTrkmEta");
  RooRealVar* BmupEta = (RooRealVar*)w.var("mupEta");
  RooRealVar* BmumEta = (RooRealVar*)w.var("mumEta");
  RooRealVar* BmumuEta = (RooRealVar*)w.var("mumuEta");
  RooRealVar* BkstEta = (RooRealVar*)w.var("kstEta");
  RooRealVar* BbEta = (RooRealVar*)w.var("bEta");
  RooRealVar* BkstTrkpPhi = (RooRealVar*)w.var("kstTrkpPhi");
  RooRealVar* BkstTrkmPhi = (RooRealVar*)w.var("kstTrkmPhi");
  RooRealVar* BmupPhi = (RooRealVar*)w.var("mupPhi");
  RooRealVar* BmumPhi = (RooRealVar*)w.var("mumPhi");
  RooRealVar* BmumuPhi = (RooRealVar*)w.var("mumuPhi");
  RooRealVar* BkstPhi = (RooRealVar*)w.var("kstPhi");
  RooRealVar* BbPhi = (RooRealVar*)w.var("bPhi");
  RooRealVar* BkstTrkpPt = (RooRealVar*)w.var("kstTrkpPt");
  RooRealVar* BkstTrkmPt = (RooRealVar*)w.var("kstTrkmPt");
  RooRealVar* BmupPt = (RooRealVar*)w.var("mupPt");
  RooRealVar* BmumPt = (RooRealVar*)w.var("mumPt");
  RooRealVar* BmumuPt = (RooRealVar*)w.var("mumuPt");
  RooRealVar* BkstPt = (RooRealVar*)w.var("kstPt");

  RooArgSet vars(*Bmass,*BkaonPt,*BpionPt,*BmmpiMass,*BmmkMass,*BmumuMass,*BbPt);
  vars.add(*Bkstarmass);
  vars.add(*BkstTrkpEta);
  vars.add(*BkstTrkmEta);
  vars.add(*BmupEta);
  vars.add(*BmumEta);
  vars.add(*BmumuEta);
  vars.add(*BkstEta);
  vars.add(*BbEta);
  vars.add(*BkstTrkpPhi);
  vars.add(*BkstTrkmPhi);
  vars.add(*BmupPhi);
  vars.add(*BmumPhi);
  vars.add(*BmumuPhi);
  vars.add(*BkstPhi);
  vars.add(*BbPhi);
  vars.add(*BkstTrkpPt);
  vars.add(*BkstTrkmPt);
  vars.add(*BmupPt);
  vars.add(*BmumPt);
  vars.add(*BmumuPt);
  vars.add(*BkstPt);

  RooDataSet* data_MC = new RooDataSet(Form("data_MC_201%i_b%i",year,q2Bin), "MC events", vars, "weight");
  RooDataSet* data_DATA = new RooDataSet(Form("data_DATA_201%i_b%i",year,q2Bin), "Data events", vars);

  int xBin;

  int counter = 0;
  int numEntries = t->GetEntries();
  for(int evt = 0; evt < numEntries; evt++){
    t->GetEntry(evt);

    // anti-radiation cut
    if (isLMNR && passB0Psi_lmnr == 0) continue;
    else if (isJpsi && passB0Psi_jpsi == 0) continue;
    else if (isPsi  && passB0Psi_psip == 0)  continue;

    // find q2 bin of current candidate
    xBin=-1;
    for (int i=0; i<nBins; ++i)
      if ( ( pow(recoDimuMass,2) < binBorders[i+1] ) &&
           ( pow(recoDimuMass,2) > binBorders[i]   ) ) {
        xBin = i;
        break;
      }
    if (xBin<0) continue;       

    // apply cut for bin 4 
    bool XCut= (( (PDGB0Mass - wt_mass) - y0Cut ) / (y1Cut-y0Cut)) < (((wt_kstarmass-PDGKstMass)-x0Cut) / (x1Cut-x0Cut)) && \
                  kaonPt > pionPt && \
                  (wt_kstarmass-PDGKstMass)>0 && \
                  (mmpiMass > CutX1 && mmpiMass < CutX2) && \
                  (mmkMass >  CutY1 && mmkMass  < CutY2) && \
                  ((mmkMass - y_0Cut) / (y_1Cut - y_0Cut)) > ((mmpiMass-x_0Cut)/(x_1Cut-x_0Cut));

    if (XCut && xBin == 4) continue;      

    if ( evt > 1.0*counter*numEntries/100 ) {
      cout<<counter<<"%"<<endl;
      counter += 10;
    }

    if((recoB0Mass > 5.0) && (recoB0Mass < 5.6)){
      Bmass->setVal(recoB0Mass);
      BkaonPt->setVal(kaonPt);
      BpionPt->setVal(pionPt);
      BmmpiMass->setVal(mmpiMass);
      BmmkMass->setVal(mmkMass);
      BmumuMass->setVal(recoDimuMass);
      BbPt->setVal(bPt);
      Bkstarmass->setVal(kstarmass);
      BkstTrkpEta->setVal(kstTrkpEta);
      BkstTrkmEta->setVal(kstTrkmEta);
      BmupEta->setVal(mupEta);
      BmumEta->setVal(mumEta);
      BmumuEta->setVal(mumuEta);
      BkstEta->setVal(kstEta);
      BbEta->setVal(bEta);
      BkstTrkpPhi->setVal(kstTrkpPhi);
      BkstTrkmPhi->setVal(kstTrkmPhi);
      BmupPhi->setVal(mupPhi);
      BmumPhi->setVal(mumPhi);
      BmumuPhi->setVal(mumuPhi);
      BkstPhi->setVal(kstPhi);
      BbPhi->setVal(bPhi);
      BkstTrkpPt->setVal(kstTrkpPt);
      BkstTrkmPt->setVal(kstTrkmPt);
      BmupPt->setVal(mupPt);
      BmumPt->setVal(mumPt);
      BmumuPt->setVal(mumuPt);
      BkstPt->setVal(kstPt);

      if(data){
        data_DATA->add(*Bmass);
        data_DATA->add(*BkaonPt);   
        data_DATA->add(*BpionPt);
        data_DATA->add(*BmmpiMass);
        data_DATA->add(*BmmkMass);
        data_DATA->add(*BmumuMass);
        data_DATA->add(*BbPt);
        data_DATA->add(*Bkstarmass);
        data_DATA->add(*BkstTrkpEta);
        data_DATA->add(*BkstTrkmEta);
        data_DATA->add(*BmupEta);
        data_DATA->add(*BmumEta);
        data_DATA->add(*BmumuEta);
        data_DATA->add(*BkstEta);
        data_DATA->add(*BbEta);
        data_DATA->add(*BkstTrkpPhi);
        data_DATA->add(*BkstTrkmPhi);
        data_DATA->add(*BmupPhi);
        data_DATA->add(*BmumPhi);
        data_DATA->add(*BmumuPhi);
        data_DATA->add(*BkstPhi);
        data_DATA->add(*BbPhi);
        data_DATA->add(*BkstTrkpPt);
        data_DATA->add(*BkstTrkmPt);
        data_DATA->add(*BmupPt);
        data_DATA->add(*BmumPt);
        data_DATA->add(*BmumuPt);
        data_DATA->add(*BkstPt);

      }
    }

  }


  if(data){
    w.import(*data_DATA);
    
    TFile* f_DATA = new TFile(Form("~/public/UML-fit/Datasets/Dataset_DATA_201%i_b%i.root",year,q2Bin), "RECREATE");
    f_DATA->cd();
    w.Write();
    f_DATA->Close();
  }
  else{
    w.import(*data_MC);
   
    TFile* f_MC = new TFile(Form("~/public/UML-fit/Datasets/Dataset_MC_201%i_b%i.root",year,q2Bin), "RECREATE");
    f_MC->cd();
    w.Write();
    f_MC->Close();
  }
}

void addVariables(RooWorkspace& w){

  RooRealVar mass("mass","mass", 5.0,5.6);
  RooRealVar mumuMass("mumuMass", "mumuMass", 2., 4.);
  RooRealVar mmkMass("mmkMass", "mmkMass", 3., 6.);
  RooRealVar mmpiMass("mmpiMass", "mmpiMass", 1.5, 5.5);
  RooRealVar pionPt("pionPt", "pionPt", 0., 40.);
  RooRealVar kaonPt("kaonPt", "kaonPt", 0., 40.);
  RooRealVar kstarmass("kstarmass", "kstarmass", 0.75, 1.05);
  RooRealVar kstTrkpEta("kstTrkpEta", "kstTrkpEta", -2.4, 2.4);
  RooRealVar kstTrkmEta("kstTrkmEta", "kstTrkmEta", -2.4, 2.4);
  RooRealVar mupEta("mupEta", "mupEta", -2.4, 2.4);
  RooRealVar mumEta("mumEta", "mumEta", -2.4, 2.4);
  RooRealVar mumuEta("mumuEta", "mumuEta", -2.4, 2.4);
  RooRealVar kstEta("kstEta", "kstEta", -2.4, 2.4);
  RooRealVar bEta("bEta", "bEta", -2.4, 2.4);
  RooRealVar kstTrkpPhi("kstTrkpPhi", "kstTrkpPhi", -3.2, 3.2);
  RooRealVar kstTrkmPhi("kstTrkmPhi", "kstTrkmPhi", -3.2, 3.2);
  RooRealVar mupPhi("mupPhi", "mupPhi", -3.2, 3.2);
  RooRealVar mumPhi("mumPhi", "mumPhi", -3.2, 3.2);
  RooRealVar mumuPhi("mumuPhi", "mumuPhi", -3.2, 3.2);
  RooRealVar kstPhi("kstPhi", "kstPhi", -3.2, 3.2);
  RooRealVar bPhi("bPhi", "bPhi", -3.2, 3.2);
  RooRealVar kstTrkpPt("kstTrkpPt", "kstTrkpPt", 0., 40.);
  RooRealVar kstTrkmPt("kstTrkmPt", "kstTrkmPt", 0., 40.);
  RooRealVar mupPt("mupPt", "mupPt", 0., 50.);
  RooRealVar mumPt("mumPt", "mumPt", 0., 60.);
  RooRealVar mumuPt("mumuPt", "mumuPt", 0., 100.);
  RooRealVar kstPt("kstPt", "kstPt", 0., 50.);
  RooRealVar bPt("bPt", "bPt", 0., 150.);

  w.import(mass);
  w.import(mumuMass);
  w.import(mmkMass);
  w.import(mmpiMass);
  w.import(pionPt);
  w.import(kaonPt);
  w.import(kstarmass);
  w.import(kstTrkpEta);
  w.import(kstTrkmEta);
  w.import(mupEta);
  w.import(mumEta);
  w.import(mumuEta);
  w.import(kstEta);
  w.import(bEta);
  w.import(kstTrkpPhi);
  w.import(kstTrkmPhi);
  w.import(mupPhi);
  w.import(mumPhi);
  w.import(mumuPhi);
  w.import(kstPhi);
  w.import(bPhi);
  w.import(kstTrkpPt);
  w.import(kstTrkmPt);
  w.import(mupPt);
  w.import(mumPt);
  w.import(mumuPt);
  w.import(kstPt);
  w.import(bPt);

}

int main(int argc, char** argv){

  bool create = true;
  bool data = true;
  int year = 2016;
  int q2Bin = -1;

  if ( (argc > 1) && (atoi(argv[1]) == 0) ) create = false;
  if ( (argc > 2) && (atoi(argv[2]) == 0) ) data = false;
  if ( argc > 3 ) year   = atoi(argv[3]);
  if ( argc > 4 ) q2Bin   = atoi(argv[4]); 

  cout <<  "create      " << create       << endl;
  cout <<  "data        " << data         << endl;
  cout <<  "year        " << year         << endl;
  cout <<  "q2Bin       " << q2Bin        << endl;

  const int nBins = 8;
  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      MC_validation(create, data, year, q2Bin);
  else
    MC_validation(create, data, year, q2Bin);

  return 0;
}





