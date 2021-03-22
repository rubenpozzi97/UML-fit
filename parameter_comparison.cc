#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>
#include <TPad.h>

#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooHist.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooPlotable.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>
#include <RooAddition.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooConstVar.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include <PdfSigRTMass.h>
#include <PdfSigWTMass.h>
#include <utils.h>
#include <RooArgSet.h>

using namespace RooFit;
using namespace std;

std::map<int,float> scale_to_data;

void parameter_comparison(int q2Bin, int parity, uint nSample, int constrain, int comp, int year){

  string all_years = Form("%i",year);
  string stat = nSample > 0 ? "_dataStat":"_MCStat";
  string component = "";
  string component1 = "";

  if(comp == 0){component = "CT";}
  else if(comp == 1){component = "WT";}

  TString input_file_Maria = "~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass" + all_years + stat + Form("_b%ip%ic%i_", q2Bin, parity,constrain) + component;

  TString input_file_Sara = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_newbdt.root",year);
  
  TFile* f_Maria = TFile::Open(input_file_Maria);
  TFile* f_Sara = TFile::Open(input_file_Sara);

  RooFitResult* fitresult_Maria = (RooFitResult*)f_Maria->Get(Form("simFitResult_b%ip%ic%isubs0",q2Bin,parity,constrain));

  if(comp == 0){component1 = "RT";}
  else if(comp == 1){component1 = "WT";} 

  RooFitResult* fitresult_Sara = (RooFitResult*)f_Sara->Get( ("results_" + component1 + Form("_%i",q2Bin)).c_str() );//change for RT/WT

  // workspace (not used)
  RooWorkspace* w;
  // rooargsets (not used)
  RooArgSet c_vars_rt, c_pdfs_rt;
  RooArgSet c_vars_wt, c_pdfs_wt;

  RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6);
  RooRealVar* rand = new RooRealVar("rand","rand", 0,1);

  RooRealVar* sigma_rt2_m;
  RooRealVar* f1rt_m;
  RooRealVar* sigma_rt2_s;
  RooRealVar* f1rt_s; 

  // Marias's parameters
  RooRealVar* mean_rt_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("mean_{RT}^{%i}",year));
  RooRealVar* sigma_rt_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#sigma_{RT1}^{%i}",year));
  RooRealVar* alpha_rt1_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#alpha_{RT1}^{%i}",year));
  RooRealVar* alpha_rt2_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#alpha_{RT2}^{%i}",year));
  RooRealVar* n_rt1_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("n_{RT1}^{%i}",year));
  RooRealVar* n_rt2_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("n_{RT2}^{%i}",year));
  if(q2Bin >= 5){
    sigma_rt2_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#sigma_{RT2}^{%i}",year));
    f1rt_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("f^{RT%i}",year));
  }

  RooRealVar* mean_wt_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("mean_{WT}^{%i}",year));
  RooRealVar* sigma_wt_m =  (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#sigma_{WT1}^{%i}",year));
  RooRealVar* alpha_wt1_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#alpha_{WT1}^{%i}",year));
  RooRealVar* alpha_wt2_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("#alpha_{WT2}^{%i}",year));
  RooRealVar* n_wt1_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("n_{WT1}^{%i}",year));
  RooRealVar* n_wt2_m = (RooRealVar*)fitresult_Maria->floatParsFinal().find(Form("n_{WT2}^{%i}",year));

  // Sara's parameters
  RooRealVar* mean_rt_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("mean_{RT}^{%i}",q2Bin));
  RooRealVar* sigma_rt_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#sigma_{RT1}^{%i}",q2Bin));
  RooRealVar* alpha_rt1_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#alpha_{RT1}^{%i}",q2Bin));
  RooRealVar* alpha_rt2_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#alpha_{RT2}^{%i}",q2Bin));
  RooRealVar* n_rt1_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("n_{RT1}^{%i}",q2Bin));
  RooRealVar* n_rt2_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("n_{RT2}^{%i}",q2Bin));
  if(q2Bin >= 5){
    sigma_rt2_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#sigma_{RT2}^{%i}",q2Bin));
    f1rt_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("f^{RT%i}",q2Bin));
  }

  RooRealVar* mean_wt_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("mean_{WT}^{%i}",q2Bin));
  RooRealVar* sigma_wt_s =  (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#sigma_{WT1}^{%i}",q2Bin));
  RooRealVar* alpha_wt1_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#alpha_{WT1}^{%i}",q2Bin));
  RooRealVar* alpha_wt2_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("#alpha_{WT2}^{%i}",q2Bin));
  RooRealVar* n_wt1_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("n_{WT1}^{%i}",q2Bin));
  RooRealVar* n_wt2_s = (RooRealVar*)fitresult_Sara->floatParsFinal().find(Form("n_{WT2}^{%i}",q2Bin));
  
  RooAbsPdf* dcb_rt_m;
  RooAbsPdf* dcb_wt_m;
  RooAbsPdf* dcb_rt_s;
  RooAbsPdf* dcb_wt_s;

  if(comp == 0){
    if(q2Bin >= 5){
      dcb_rt_m = createRTMassShape(q2Bin, mass, mean_rt_m, sigma_rt_m, sigma_rt2_m, alpha_rt1_m, alpha_rt2_m, n_rt1_m, n_rt2_m, f1rt_m, w, year, false, c_vars_rt, c_pdfs_rt);  
      dcb_rt_s = createRTMassShape(q2Bin, mass, mean_rt_s, sigma_rt_s, sigma_rt2_s, alpha_rt1_s, alpha_rt2_s, n_rt1_s, n_rt2_s, f1rt_s, w, year, false, c_vars_rt, c_pdfs_rt);
    } 
    else{
      dcb_rt_m = createRTMassShape(q2Bin, mass, mean_rt_m, sigma_rt_m, alpha_rt1_m, alpha_rt2_m, n_rt1_m, n_rt2_m, w, year, false, c_vars_rt, c_pdfs_rt);
      dcb_rt_s = createRTMassShape(q2Bin, mass, mean_rt_s, sigma_rt_s, alpha_rt1_s, alpha_rt2_s, n_rt1_s, n_rt2_s, w, year, false, c_vars_rt, c_pdfs_rt);
    }
  }
  else if(comp == 1){
    dcb_wt_m = createWTMassShape(q2Bin, mass, mean_wt_m, sigma_wt_m, alpha_wt1_m, alpha_wt2_m, n_wt1_m, n_wt2_m, w, year, false, c_vars_wt, c_pdfs_wt );
    dcb_wt_s = createWTMassShape(q2Bin, mass, mean_wt_s, sigma_wt_s, alpha_wt1_s, alpha_wt2_s, n_wt1_s, n_wt2_s, w, year, false, c_vars_wt, c_pdfs_wt );
  }

 RooPlot* massframe = mass->frame(Title(Form("q^{2} bin %i, %i efficiency, year %i",q2Bin,parity,year)));

 if(comp == 0){
   dcb_rt_m->plotOn(massframe, RooFit::Name("Maria RT"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
   dcb_rt_s->plotOn(massframe, RooFit::Name("Sara RT"),Range("all"),LineColor(kBlue),LineStyle(1),LineWidth(2));
   dcb_rt_m->paramOn(massframe,Layout(0.1,0.4,0.9));
   massframe->getAttText()->SetTextSize(0.02) ;
   dcb_rt_s->paramOn(massframe,Layout(0.6,0.9,0.9));
   massframe->getAttText()->SetTextSize(0.02) ;
 }
 else if(comp == 1){
   dcb_wt_m->plotOn(massframe, RooFit::Name("Maria WT"),Range("all"),LineColor(kRed),LineStyle(1),LineWidth(2));
   dcb_wt_s->plotOn(massframe, RooFit::Name("Sara WT"),Range("all"),LineColor(kBlue),LineStyle(1),LineWidth(2));
   dcb_wt_m->paramOn(massframe,Layout(0.1,0.4,0.9));
   massframe->getAttText()->SetTextSize(0.02) ;
   dcb_wt_s->paramOn(massframe,Layout(0.6,0.9,0.9));
   massframe->getAttText()->SetTextSize(0.02) ;

 }

 TCanvas c;
 c.SetTitle("");

 massframe->GetYaxis()->SetTitleOffset(1.3);
 massframe->SetXTitle("mass (GeV)");

 massframe->Draw();

 TLegend *leg = new TLegend (0.6,0.2,0.9,0.4);
 if(comp == 0){
   leg->AddEntry(massframe->findObject("Maria RT"), "Maria RT", "lep");
   leg->AddEntry(massframe->findObject("Sara RT"), "Sara RT", "lep");
 }
 else if(comp == 1){
   leg->AddEntry(massframe->findObject("Maria WT"), "Maria WT", "lep");
   leg->AddEntry(massframe->findObject("Sara WT"), "Sara WT", "lep");
 }
 leg->SetTextSize(0.03);
 leg->Draw("same");

 if(comp == 0){
   c.SaveAs(Form("plot_params/param_comp_RT_b%ip%i_%i.pdf",q2Bin,parity,year));
   c.SaveAs(Form("plot_params/param_comp_RT_b%ip%i_%i.gif",q2Bin,parity,year));}
 else if(comp == 1){
   c.SaveAs(Form("plot_params/param_comp_WT_b%ip%i_%i.pdf",q2Bin,parity,year));
   c.SaveAs(Form("plot_params/param_comp_WT_b%ip%i_%i.gif",q2Bin,parity,year));}

  // Data dataset
  std::vector<RooWorkspace*> wsp;
  std::vector<std::vector<RooDataSet*>> data;
  string filename_data = Form("recoMCDataset_b%i_%i.root", q2Bin, year);
  filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/newphi/", year) + filename_data;

  if (!retrieveWorkspace(filename_data, wsp, Form("ws_b%ip%i", q2Bin, 1-parity), wsp, Form("ws_b%ip%i", q2Bin, 1-parity ), parity))  return;

  uint firstSample = 0;
  uint lastSample = 0;
  RooArgSet reco_vars (*mass, *rand);
 
  string shortString = Form("b%ip%i",q2Bin,parity);

  data.push_back(createDataset(nSample,  firstSample,  lastSample, wsp[0], wsp[0],
                 q2Bin,  parity,  year,
                 reco_vars,  shortString, comp));

  std::map<std::string, RooDataSet*> map;
  map.insert(map.cbegin(), std::pair<const string,RooDataSet*>("data", data[0][0]));

  RooCategory sample ("sample", "sample");
  string isample = "";
  isample.clear(); isample.assign( Form("%i",0) );
  string years = Form("%i",year);
  sample.defineType(("data"+years+"_subs"+isample).c_str());


  RooDataSet allcombData ("allcombData", "combined data",
                          reco_vars,
                          Index(sample),
                          Import(map));

  TCanvas d;
  d.cd();
 
  TPad *p1 = new TPad("p1","p1",0.,0.27,1.,1.);
  p1->SetBorderMode(1);
  p1->SetFrameBorderMode(0);
  p1->SetBorderSize(2);
  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.,0.065,1.,0.24);
  p2->SetBorderMode(1);
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.4);
  p2->Draw();

  TLegend *legs = new TLegend (0.6,0.7,0.9,0.9);

  RooPlot* frames = mass->frame(Title(Form("q^{2} bin %i, %i efficiency, year %i",q2Bin,parity,year)));

  allcombData.plotOn(frames, MarkerColor(kRed+1), LineColor(kRed+1), Binning(80),Name("Data"));

  if(comp == 0){dcb_rt_s->plotOn(frames, Name("Sara RT"),Range("all"),LineColor(kBlue),LineStyle(1),LineWidth(2));}
  else if(comp == 1){dcb_wt_s->plotOn(frames, Name("Sara WT"),Range("all"),LineColor(kBlue),LineStyle(1),LineWidth(2));}

  legs->AddEntry(frames->findObject("Data"),"Data" ,"lep");

  if(comp == 0){legs->AddEntry(frames->findObject("Sara RT"),"Sara RT" ,"lep");}
  else if(comp == 1){legs->AddEntry(frames->findObject("Sara WT"),"Sara WT" ,"lep");}

  p1->cd();
  frames->SetXTitle("mass (GeV)");
  frames->Draw();
  legs->Draw("same");

  /*
  RooHist* pull_hist;
  if(comp == 0){RooHist* pull_hist = frames->pullHist("Data","Sara RT");}
  else if(comp == 1){RooHist* pull_hist = frames->pullHist("Data","Sara WT");}
 
  RooPlot *pull_plot = mass->frame();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  
  pull_plot->SetTitle("");
  pull_plot->GetXaxis()->SetTitle("mass (GeV)");
  pull_plot->GetXaxis()->SetTitleSize(0.15);
  pull_plot->GetXaxis()->SetTitleOffset(0.9);
  pull_plot->GetXaxis()->SetLabelSize(0.15);
  pull_plot->GetXaxis()->SetTickLength(0.1);
  pull_plot->GetYaxis()->SetTitle("Pull");
  pull_plot->GetYaxis()->SetTitleSize(0.15);
  pull_plot->GetYaxis()->SetTitleOffset(0.1);
  pull_plot->GetYaxis()->SetLabelSize(0.15);
  pull_plot->GetYaxis()->SetNdivisions(305);
  
  gPad->Update();
  TLine *line = new TLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);
  line->SetLineStyle(2);
  line->SetLineColor(kBlue);

  p2->cd();
  pull_plot->Draw();
  line->Draw("same");
  */

  if(comp == 0){
    d.SaveAs(Form("params_results/data_sara_RT_b%ip%i_%i.pdf",q2Bin,parity,year));
    d.SaveAs(Form("params_results/data_sara_RT_b%ip%i_%i.gif",q2Bin,parity,year));}
  else if(comp == 1){
    d.SaveAs(Form("params_results/data_sara_WT_b%ip%i_%i.pdf",q2Bin,parity,year));
    d.SaveAs(Form("params_results/data_sara_WT_b%ip%i_%i.gif",q2Bin,parity,year));}

}

int main(int argc, char** argv){

  int q2Bin   = -1;
  int parity  = -1;

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);
  if ( argc > 2 ) parity  = atoi(argv[2]);

  uint nSample = 0;
  if ( argc > 3 ) nSample = atoi(argv[3]);

  bool constrain = true;
  if ( argc > 4 && atoi(argv[4]) == 0 ) constrain = false;

  int comp = 0;
  if ( argc > 5 ) comp = atoi(argv[5]);

  int year = 0;
  if ( argc > 6 && atoi(argv[6]) != 0 ) year = atoi(argv[6]);

  scale_to_data.insert(std::make_pair(2016, 0.006*2 /2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
  scale_to_data.insert(std::make_pair(2017, 0.005*2 /2.05 ));
  scale_to_data.insert(std::make_pair(2018, 0.007*2 /1.9  ));

  parameter_comparison(q2Bin, parity, nSample, constrain, comp, year);

  return 0;
}












