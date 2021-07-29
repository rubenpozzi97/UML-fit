#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TAttFill.h>
#include <TFrame.h>
#include <TStyle.h>

#include <TLatex.h>
#include <list>
#include <map>

#include <RooGenericPdf.h>
#include <RooPolynomial.h>
#include <RooExponential.h>
#include <RooProduct.h>
#include <RooFormulaVar.h>
#include <RooArgList.h>
#include <RooAbsPdf.h>
#include <RooExponential.h>
#include <RooHist.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooPlot.h>
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

#include "PdfSigMass.h"
#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;
std::map<int,float> scale_to_data;

TCanvas* c [4*nBins];

void simfit_recoMC_fullMassBin(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, int constrain, int comp, int dat, int pdf_model, std::vector<int> years)
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%ic%im%i",q2Bin,parity,constrain,pdf_model);
  cout<<"Conf: "<<shortString<<endl;

  string dataString;
  if(dat == 0){dataString = "MC";}
  else if(dat == 1){dataString = "DATA";}

  string all_years = "";
  string year = ""; 
  string isample = "";

  string stat;
  if( (nSample == 0) && (dat == 0) ){stat = "_MCStat";}
  else if( (nSample == 1) && (dat == 0) ){stat = "_dataStat";}
  else if( (nSample == 0) && (dat == 1) ){stat = "_DATA";}
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;

  string component = "";
  if(comp == 0){component = "CT";}
  else if(comp == 1){component = "WT";}
  else{component = "CT+WT";}

  std::vector<RooWorkspace*> wsp, wsp_mcmass, wsp_even, wsp_odd;
  std::vector<std::vector<RooDataSet*>> data;

  std::vector<RooGaussian*> c_deltaPeaks, c_fm;
  RooArgSet c_vars_rt, c_pdfs_rt;
  RooArgSet c_vars_wt, c_pdfs_wt;
  RooArgSet c_vars; 
  RooWorkspace * ws_pars = new RooWorkspace("ws_pars");

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

  RooRealVar* mass = 0;
  if(comp == 0){
    mass = new RooRealVar("mass","mass", 5.,5.6);
  }
  else if(comp == 1){
    mass = new RooRealVar("mass","mass", 4.9,5.7);
  }
  else if((q2Bin == 4) && (pdf_model == 0)){
    mass = new RooRealVar("mass","mass", 5.1,5.45);
  }
  else if(pdf_model == 1){
    mass = new RooRealVar("mass","mass", 5.1,5.6);
  }
  else if(pdf_model == 2){
    mass = new RooRealVar("mass","mass", 5.0,5.5);
  }
  else{
    mass = new RooRealVar("mass","mass", 5.,5.6);
  }
  ws_pars->import(*mass);

  RooRealVar* rand = new RooRealVar("rand","rand", 0,1);
  RooArgSet reco_vars (*mass, *rand);

  RooCategory sample ("sample", "sample");
  for(unsigned int iy = 0; iy < years.size(); iy++){
    year.clear(); year.assign(Form("%i",years[iy]));
    all_years += year;
    for(uint is = firstSample; is <= lastSample; is++){
      isample.clear(); isample.assign( Form("%i",is) );
      sample.defineType(("data"+year+"_subs"+isample).c_str());
    }
  }

  // Construct a simultaneous pdf using category sample as index (simultaneous over the 3 years)
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = ("reco" + dataString + Form("Dataset_b%i_%i.root", q2Bin, years[iy])).c_str(); 
 
    // import data (or MC as data proxy) HERE
    if(parity < 2){if (!retrieveWorkspace(filename_data, wsp, Form("ws_b%ip%i", q2Bin, 1-parity ), wsp, Form("ws_b%ip%i", q2Bin, 1-parity ), parity))  return;}
    else{if (!retrieveWorkspace(filename_data, wsp_even, Form("ws_b%ip%i", q2Bin, 0), wsp_odd, Form("ws_b%ip%i", q2Bin, 1), parity))  return;}

    if(parity < 2){
        data.push_back(createDataset(nSample,  firstSample,  lastSample, wsp[iy], wsp[iy],
                                     q2Bin,  parity,  years[iy],
                                     reco_vars,  shortString, comp, dat));
    }
    else{
        data.push_back(createDataset(nSample,  firstSample,  lastSample, wsp_even[iy], wsp_odd[iy],
                                     q2Bin,  parity,  years[iy],
                                     reco_vars,  shortString, comp, dat));
    }
  
    // import mass PDF from fits to the MC 
    if(constrain == 0){
      if(q2Bin == 4){
        string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_Jpsi_newbdt.root",years[iy]);
        if (!retrieveWorkspace(filename_mc_mass, wsp_mcmass, "w", wsp_mcmass, "w", 0))  return;
      }
      else if(q2Bin == 6){
        string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_Psi_newbdt.root",years[iy]);
        if (!retrieveWorkspace(filename_mc_mass, wsp_mcmass, "w", wsp_mcmass, "w", 0))  return;
      }
      else{
        string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_newbdt.root",years[iy]);
        if (!retrieveWorkspace(filename_mc_mass, wsp_mcmass, "w", wsp_mcmass, "w", 0))  return;
      }  
    }

    TString input_file_RT;
    TString input_file_WT;
    TFile* f_RT;
    TFile* f_WT;
    RooFitResult* fitresult_RT;
    RooFitResult* fitresult_WT;

    if(constrain > 0){
      input_file_RT = Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0CT.root", years[iy], q2Bin);
      input_file_WT = Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0WT.root", years[iy], q2Bin);
      f_RT = TFile::Open(input_file_RT);
      f_WT = TFile::Open(input_file_WT);
      fitresult_RT = (RooFitResult*)f_RT->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));
      fitresult_WT = (RooFitResult*)f_WT->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));
    }

    RooRealVar* mean_rt = 0;
    RooRealVar* sigma_rt = 0;
    RooRealVar* alpha_rt1 = 0;
    RooRealVar* alpha_rt2 = 0;
    RooRealVar* n_rt1 = 0;
    RooRealVar* n_rt2 = 0;
    RooRealVar* sigma_rt2 = 0;
    RooRealVar* f1rt = 0;
    RooRealVar* mean_wt = 0;
    RooFormulaVar* mass_wt = 0;
    RooRealVar* mean_difference = 0;
    RooRealVar* sigma_wt = 0;
    RooRealVar* alpha_wt1 = 0;
    RooRealVar* alpha_wt2 = 0;
    RooRealVar* n_wt1 = 0;
    RooRealVar* n_wt2 = 0;
    RooRealVar* signal_yield;
    RooRealVar* CB_yield;
    RooRealVar* lambda;
    RooRealVar* slope;
    RooRealVar* mean_cb;
    RooRealVar* sigma_cb;
    RooRealVar* cofs;
    RooRealVar* scale;
    RooRealVar* shift;
    RooRealVar* f_erf;
 
    // variables for fitting data with parameters fixed to MC and with scale factor applied to the sigmas
    RooRealVar* factor = new RooRealVar(Form("factor^{%i}",years[iy]),"factor",1,0,2); //scale factor, used when constrain=2
    ws_pars->import(*factor);
    RooProduct* sigma_rt_fix;
    RooProduct* sigma_rt2_fix;
    RooProduct* mean_difference_fix;
    RooProduct* sigma_wt_fix;

    RooAbsPdf* dcb_rt;
    RooProdPdf* c_dcb_rt;
    RooAbsPdf* dcb_wt;
    RooProdPdf* c_dcb_wt;
    RooAddPdf* signal_PDF;
    RooProdPdf* final_PDF;
    RooAddPdf* final_PDF_extended;

    // yields
    double mass_peak = 5.27958; //PDG B0 mass
    int signal_yield_initial;
    int CB_yield_initial;

    if(dat == 0){//MC
      signal_yield_initial = data[iy][0]->numEntries();
      signal_yield = new RooRealVar(Form("sig_yield^{%i}",years[iy]), "signal_yield", signal_yield_initial, 0, 2*signal_yield_initial);
      ws_pars->import(*signal_yield);
    }
    else if(dat == 1){//DATA
      signal_yield_initial = data[iy][0]->sumEntries(TString::Format("abs(mass-%g)<0.10", mass_peak));
      signal_yield = new RooRealVar(Form("sig_yield^{%i}",years[iy]), "signal_yield", signal_yield_initial, 0, data[iy][0]->numEntries());
      CB_yield_initial = data[iy][0]->numEntries() - signal_yield_initial;
      CB_yield = new RooRealVar(Form("CB_yield^{%i}",years[iy]), "CB_yield", CB_yield_initial, 0, data[iy][0]->numEntries());
      ws_pars->import(*signal_yield);
      ws_pars->import(*CB_yield);
    }

    // create RT component 
    if(constrain == 0){
      wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
      mean_rt = new RooRealVar(Form("mean_{RT}^{%i}",years[iy]), "massrt", wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal(), 5, 6, "GeV");
      if(comp == 1){mean_rt->setConstant();}
      sigma_rt = new RooRealVar(Form("#sigma_{RT1}^{%i}",years[iy]), "sigmart1", wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal(), 0, 1,"GeV");
      alpha_rt1 = new RooRealVar(Form("#alpha_{RT1}^{%i}",years[iy]), "alphart1", wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal(), 0, 10);
      n_rt1 = new RooRealVar(Form("n_{RT1}^{%i}",years[iy]), "nrt1", wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal(), 0., 200.);
      if(q2Bin != 7){
        alpha_rt2 = new RooRealVar(Form("#alpha_{RT2}^{%i}",years[iy]), "alphart2", wsp_mcmass[iy]->var(Form("#alpha_{RT2}^{%i}", q2Bin))->getVal(), -10, 10);
        n_rt2 = new RooRealVar(Form("n_{RT2}^{%i}",years[iy]), "nrt2", wsp_mcmass[iy]->var(Form("n_{RT2}^{%i}", q2Bin))->getVal(), 0., 200.);
        ws_pars->import(*alpha_rt2);
        ws_pars->import(*n_rt2);
      }
    }
    else if(constrain > 0){
      mean_rt = new RooRealVar(Form("mean_{RT}^{%i}",years[iy]), "massrt", MC_fit_result(input_file_RT, Form("mean_{RT}^{%i}",years[iy]), q2Bin), 5, 6, "GeV");
      if(comp == 1){mean_rt->setConstant();}
      sigma_rt = new RooRealVar(Form("#sigma_{RT1}^{%i}",years[iy]), "sigmart1", MC_fit_result(input_file_RT, Form("#sigma_{RT1}^{%i}",years[iy]), q2Bin), 0, 1,"GeV");
      alpha_rt1 = new RooRealVar(Form("#alpha_{RT1}^{%i}",years[iy]), "alphart1", MC_fit_result(input_file_RT, Form("#alpha_{RT1}^{%i}",years[iy]), q2Bin), 0, 10);
      n_rt1 = new RooRealVar(Form("n_{RT1}^{%i}",years[iy]), "nrt1", MC_fit_result(input_file_RT, Form("n_{RT1}^{%i}",years[iy]), q2Bin), 0., 200.);
      if(q2Bin != 7){
        alpha_rt2 = new RooRealVar(Form("#alpha_{RT2}^{%i}",years[iy]), "alphart2", MC_fit_result(input_file_RT, Form("#alpha_{RT2}^{%i}",years[iy]), q2Bin), -10, 10);
        n_rt2 = new RooRealVar(Form("n_{RT2}^{%i}",years[iy]), "nrt2", MC_fit_result(input_file_RT, Form("n_{RT2}^{%i}",years[iy]), q2Bin), 0., 200.);      
        ws_pars->import(*alpha_rt2);
        ws_pars->import(*n_rt2);
      }
    }

    ws_pars->import(*mean_rt);
    ws_pars->import(*sigma_rt);
    ws_pars->import(*alpha_rt1);
    ws_pars->import(*n_rt1);

    if( (q2Bin == 4) || (q2Bin == 5) || (q2Bin == 6) ){//CB+CB
      if(constrain == 0){
        sigma_rt2 = new RooRealVar(Form("#sigma_{RT2}^{%i}",years[iy]), "sigmaRT2", wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal(), 0,   0.12, "GeV");
        f1rt = new RooRealVar(Form("f^{RT%i}",years[iy]), "f1rt", wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal(), 0,  1.);
      }
      else if(constrain > 0){
        sigma_rt2 = new RooRealVar(Form("#sigma_{RT2}^{%i}",years[iy]), "sigmaRT2", MC_fit_result(input_file_RT, Form("#sigma_{RT2}^{%i}",years[iy]), q2Bin), 0,   0.12, "GeV");
        f1rt = new RooRealVar(Form("f^{RT%i}",years[iy]), "f1rt", MC_fit_result(input_file_RT, Form("f^{RT%i}",years[iy]), q2Bin), 0,  1.);
      }

      ws_pars->import(*sigma_rt2);
      ws_pars->import(*f1rt);
      alpha_rt2->setRange(-10,0);
     
      if(constrain == 1){//constrained fit 
        dcb_rt = createRTMassShape2(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2, f1rt, input_file_RT, years[iy], true, c_vars_rt, c_pdfs_rt);
      } 
      else if(constrain == 0){//unconstrained fit
        dcb_rt = createRTMassShape2(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2, f1rt, input_file_RT, years[iy], false, c_vars_rt, c_pdfs_rt);
      }
      else if(constrain == 2){//fit parameters fixed to MC, scale factor applied to sigmas and mean difference
        sigma_rt->setConstant();
        sigma_rt2->setConstant();
        alpha_rt1->setConstant();
        alpha_rt2->setConstant();
        n_rt1->setConstant();
        n_rt2->setConstant();
        f1rt->setConstant();

        sigma_rt_fix = new RooProduct(Form("sigma_rt_fix^{%i}",years[iy]),"sigma_rt_fix",RooArgList(*factor,*sigma_rt));         
        sigma_rt2_fix = new RooProduct(Form("sigma_rt2_fix^{%i}",years[iy]),"sigma_rt2_fix",RooArgList(*factor,*sigma_rt2)); 

        dcb_rt = createRTMassShape2(q2Bin, mass, mean_rt, sigma_rt_fix, sigma_rt2_fix, alpha_rt1, alpha_rt2, n_rt1, n_rt2, f1rt, input_file_RT, years[iy], false, c_vars_rt, c_pdfs_rt);
      }
    }

    else if(q2Bin == 7){//CB+gauss
      if(constrain == 0){
        sigma_rt2 = new RooRealVar (Form("#sigma_{RT2}^{%i}",years[iy] ), "sigmaRT2"  ,    wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal(), 0,   0.12, "GeV");
        f1rt      = new RooRealVar (Form("f^{RT%i}",years[iy])          , "f1rt"      ,   wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal(), 0,  1.);
      }
      else if(constrain > 0){
       sigma_rt2 = new RooRealVar(Form("#sigma_{RT2}^{%i}",years[iy]), "sigmaRT2", MC_fit_result(input_file_RT, Form("#sigma_{RT2}^{%i}",years[iy]), q2Bin), 0,   0.12, "GeV");
       f1rt = new RooRealVar(Form("f^{RT%i}",years[iy]), "f1rt", MC_fit_result(input_file_RT, Form("f^{RT%i}",years[iy]), q2Bin), 0,  1.);
      }

      ws_pars->import(*sigma_rt2);
      ws_pars->import(*f1rt);

      if(constrain == 1){//constrained fit
        dcb_rt = createRTMassShape3(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1, f1rt, input_file_RT, years[iy], true, c_vars_rt, c_pdfs_rt);
      }
      else if(constrain == 0){//uncsontrained fit
        dcb_rt = createRTMassShape3(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, n_rt1, f1rt, input_file_RT, years[iy], false, c_vars_rt, c_pdfs_rt);
      }
      else if(constrain == 2){//fit parameters fixed to MC, scale factor applied to sigmas and mean difference
        sigma_rt->setConstant();
        sigma_rt2->setConstant();
        alpha_rt1->setConstant();
        n_rt1->setConstant();
        f1rt->setConstant();
        
        sigma_rt_fix = new RooProduct(Form("sigma_rt_fix^{%i}",years[iy]),"sigma_rt_fix",RooArgList(*factor,*sigma_rt));
        sigma_rt2_fix = new RooProduct(Form("sigma_rt2_fix^{%i}",years[iy]),"sigma_rt2_fix",RooArgList(*factor,*sigma_rt2));

        dcb_rt = createRTMassShape3(q2Bin, mass, mean_rt, sigma_rt_fix, sigma_rt2_fix, alpha_rt1, n_rt1, f1rt, input_file_RT, years[iy], false, c_vars_rt, c_pdfs_rt);
      }
    }

    else{//double CB
      alpha_rt2->setRange(0,10);
       
      if(constrain == 1){//constrained fit 
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2, input_file_RT, years[iy], true, c_vars_rt, c_pdfs_rt);
      }
      else if(constrain == 0){//unconstrained fit
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2, input_file_RT, years[iy], false, c_vars_rt, c_pdfs_rt);
      }
      else if(constrain == 2){//fit parameters fixed to MC, scale factor applied to sigmas and mean difference
        sigma_rt->setConstant();
        alpha_rt1->setConstant();
        alpha_rt2->setConstant();
        n_rt1->setConstant();
        n_rt2->setConstant();
      
        sigma_rt_fix = new RooProduct(Form("sigma_rt_fix^{%i}",years[iy]),"sigma_rt_fix",RooArgList(*factor,*sigma_rt));

        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt_fix, alpha_rt1, alpha_rt2, n_rt1, n_rt2, input_file_RT, years[iy], false, c_vars_rt, c_pdfs_rt);
      }
    }

    // create WT component
    if(constrain == 0){
      wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));
      mean_wt = new RooRealVar(Form("mean_{WT}^{%i}",years[iy]), "masswt", wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal(), 5, 6, "GeV");
      sigma_wt = new RooRealVar(Form("#sigma_{WT1}^{%i}",years[iy]), "sigmawt", wsp_mcmass[iy]->var(Form("#sigma_{WT1}^{%i}", q2Bin))->getVal(), 0, 1, "GeV");
      alpha_wt1 = new RooRealVar(Form("#alpha_{WT1}^{%i}",years[iy]), "alphawt1", wsp_mcmass[iy]->var(Form("#alpha_{WT1}^{%i}", q2Bin))->getVal(), 0, 10);
      alpha_wt2 = new RooRealVar(Form("#alpha_{WT2}^{%i}",years[iy]), "alphawt2", wsp_mcmass[iy]->var(Form("#alpha_{WT2}^{%i}", q2Bin))->getVal(), 0, 10);
      n_wt1 = new RooRealVar(Form("n_{WT1}^{%i}",years[iy]), "nwt1", wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal(), 0., 100.);
      n_wt2 = new RooRealVar(Form("n_{WT2}^{%i}",years[iy]), "nwt2", wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal(), 0., 100.);
    }
    else if(constrain > 0){
      sigma_wt = new RooRealVar(Form("#sigma_{WT1}^{%i}",years[iy]), "sigmawt", MC_fit_result(input_file_WT, Form("#sigma_{WT1}^{%i}",years[iy]), q2Bin), 0, 1, "GeV");
      alpha_wt1 = new RooRealVar(Form("#alpha_{WT1}^{%i}",years[iy]), "alphawt1", MC_fit_result(input_file_WT, Form("#alpha_{WT1}^{%i}",years[iy]), q2Bin), 0, 10);
      alpha_wt2 = new RooRealVar(Form("#alpha_{WT2}^{%i}",years[iy]), "alphawt2", MC_fit_result(input_file_WT, Form("#alpha_{WT2}^{%i}",years[iy]), q2Bin), 0, 10);
      n_wt1 = new RooRealVar(Form("n_{WT1}^{%i}",years[iy]), "nwt1", MC_fit_result(input_file_WT, Form("n_{WT1}^{%i}",years[iy]), q2Bin), 0., 100.);
      n_wt2 = new RooRealVar(Form("n_{WT2}^{%i}",years[iy]), "nwt2", MC_fit_result(input_file_WT, Form("n_{WT2}^{%i}",years[iy]), q2Bin), 0., 100.);
      mean_difference = new RooRealVar (Form("mean_difference^{%i}",years[iy]), "mean_difference", MC_fit_result(input_file_WT, Form("mean_{WT}^{%i}",years[iy]), q2Bin) - MC_fit_result(input_file_RT, Form("mean_{RT}^{%i}",years[iy]), q2Bin), -2, 2, "GeV");
    }

    ws_pars->import(*sigma_wt);
    ws_pars->import(*alpha_wt1);
    ws_pars->import(*alpha_wt2);
    ws_pars->import(*n_wt1);
    ws_pars->import(*n_wt2);
    if(constrain == 1){ws_pars->import(*mean_difference);}

    if(constrain == 1){//constrained fit
      mass_wt = new RooFormulaVar(Form("mass_{WT}^{%i}",years[iy]), "meanwt", "@0+@1", RooArgList(*mean_rt,*mean_difference)); 

      if( ((pdf_model != 0) && (q2Bin == 4)) || ((pdf_model == 3) && (q2Bin != 4) && (q2Bin != 6))){ // fix WT shape of MC to compute systematics
        sigma_wt->setConstant();
        alpha_wt1->setConstant();
        alpha_wt2->setConstant();
        n_wt1->setConstant();
        n_wt2->setConstant();
     
        dcb_wt = createWTMassShape(q2Bin, mass, mass_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2, input_file_WT, years[iy], false, c_vars_wt, c_pdfs_wt);
      }
      else{
        dcb_wt = createWTMassShape(q2Bin, mass, mass_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2, input_file_WT, years[iy], true, c_vars_wt, c_pdfs_wt);
      }
    }
    else if(constrain == 0){//unconstrained fit
      dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2, input_file_WT, years[iy], false, c_vars_wt, c_pdfs_wt);
    }
    else if(constrain == 2){//fit parameters fixed to MC, scale factor applied to sigmas and mean difference
      mean_difference->setConstant();
      mean_difference_fix = new RooProduct(Form("mean_difference_fix^{%i}",years[iy]),"mean_difference_fix",RooArgList(*factor,*mean_difference));
      mass_wt = new RooFormulaVar(Form("mean_{WT}^{%i}",years[iy]), "masswt", "@0+@1", RooArgList(*mean_rt,*mean_difference_fix));

      sigma_wt->setConstant();
      alpha_wt1->setConstant();
      alpha_wt2->setConstant();
      n_wt1->setConstant();
      n_wt2->setConstant(); 

      sigma_wt_fix = new RooProduct(Form("sigma_wt_fix^{%i}",years[iy]),"sigma_wt_fix",RooArgList(*factor,*sigma_wt));

      dcb_wt = createWTMassShape(q2Bin, mass, mass_wt, sigma_wt_fix, alpha_wt1, alpha_wt2, n_wt1, n_wt2, input_file_WT, years[iy], false, c_vars_wt, c_pdfs_wt);
    }  

    if(constrain == 1){ // fit with Gaussian constraints
      /// create constrained PDF for RT mass
      RooArgList constr_rt_list = RooArgList(c_pdfs_rt);
      constr_rt_list.add(*dcb_rt);
      c_dcb_rt = new RooProdPdf(("c_dcb_rt_"+year).c_str(), ("c_dcb_rt_"+year).c_str(), constr_rt_list);
      c_vars.add(c_vars_rt);

      ///mean_difference constraints for data-like MC and data
      RooRealVar* constrained_mean_RT = (RooRealVar*)fitresult_RT->floatParsFinal().find(Form("mean_{RT}^{%i}",years[iy]));
      RooRealVar* constrained_mean_WT = (RooRealVar*)fitresult_WT->floatParsFinal().find(Form("mean_{WT}^{%i}",years[iy]));

      double mean_diff_val = (constrained_mean_WT->getVal()) - (constrained_mean_RT->getVal()); 
      double mean_diff_error = sqrt( pow((constrained_mean_WT->getError()),2) + pow((constrained_mean_RT->getError()),2) );

      cout << "mean_diff_val = " << mean_diff_val << endl;
      cout << "mean_diff_error = " <<  mean_diff_error << endl;

      RooGaussian* mean_constr = new RooGaussian(Form("c_mean_diff_%i", years[iy]) ,
                                              Form("c_mean_diff_%i", years[iy]),
                                              *mean_difference,
                                              RooConst(mean_diff_val),
                                              RooConst(mean_diff_error)
                                              );

      /// create constrained PDF for WT mass   
      RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
      constr_wt_list.add(*mean_constr);
      constr_wt_list.add(*dcb_wt);
      c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list);
      c_vars.add(c_vars_wt);
      c_vars.add(*mean_difference);

      /// create constraint on mFrac (here there is no efficiency, therefore value set to measured value on MC)
      RooRealVar* constrained_yield_RT = (RooRealVar*)fitresult_RT->floatParsFinal().find(Form("sig_yield^{%i}",years[iy]));
      RooRealVar* constrained_yield_WT = (RooRealVar*)fitresult_WT->floatParsFinal().find(Form("sig_yield^{%i}",years[iy])); 
 
      double nrt_mc = constrained_yield_RT->getVal();
      double nwt_mc = constrained_yield_WT->getVal();
      double fraction = nwt_mc / (nrt_mc + nwt_mc);

      RooRealVar* mFrac = new RooRealVar(Form("mFrac^{%i}",years[iy]),"mistag fraction",fraction, 0, 1);
      ws_pars->import(*mFrac);
      c_fm.push_back( new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,
                                      RooConst(fraction),
                                      RooConst(fM_sigmas[years[iy]][q2Bin])
                                      ));

      cout << fraction << "   " << fM_sigmas[years[iy]][q2Bin] << endl;
      c_vars.add(*mFrac);

      signal_PDF = new RooAddPdf(("signal_PDF_"+year).c_str(),
				 ("signal_PDF_"+year).c_str(),
				  RooArgList(*c_dcb_wt,*c_dcb_rt), *mFrac);
      final_PDF = new RooProdPdf(("final_PDF_"+year).c_str(),
                                ("final_PDF_"+year).c_str(),
                                 RooArgList(*signal_PDF, *c_fm[iy]));                         
      if(dat == 0){
        final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
      					   ("final_PDF_extended_"+year).c_str(),
					    *final_PDF, *signal_yield);
      }
      else if(dat == 1){
        lambda = new RooRealVar(Form("lambda^{%i}",years[iy]), "lambda", 1., -15., 15.);
        slope = new RooRealVar(Form("slope^{%i}",years[iy]), "slope", 1., -5., 5.);
        mean_cb = new RooRealVar(Form("mean_cb^{%i}",years[iy]), "mean_cb", 5.05, 5., 5.15);
        sigma_cb = new RooRealVar(Form("sigma_cb^{%i}",years[iy]), "sigma_cb", 1., 0., 5.);
        cofs = new RooRealVar(Form("cofs^{%i}",years[iy]), "cofs", 0.5, 0., 1.);
        scale = new RooRealVar("scale", "scale", 12., 1., 20.);
        shift = new RooRealVar("shift", "shift", 5.13, 4.7, 5.2);
        f_erf = new RooRealVar("f_erf", "f_erf", 0.5, 0., 1.);

        /*if( (pdf_model == 3) && (q2Bin != 6) && (q2Bin != 4) ){
          RooExponential* exp_bkg = new RooExponential(("exp_bkg_PDF_"+year).c_str(),
                                                       ("exp_bkg_PDF_"+year).c_str(),
                                                        *mass, *lambda);
          RooPolynomial* poly_bkg = new RooPolynomial(("poly_bkg_PDF_"+year).c_str(),
                                                      ("poly_bkg_PDF_"+year).c_str(),
                                                       *mass, *slope);
        
          RooAddPdf* CB_bkg = new RooAddPdf(("CB_bkg_PDF_"+year).c_str(),
                                            ("CB_bkg_PDF_"+year).c_str(),
                                             RooArgList(*exp_bkg,*poly_bkg), *cofs);
          final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
                                             ("final_PDF_extended_"+year).c_str(),
                                              RooArgList(*final_PDF, *CB_bkg), RooArgList(*signal_yield, *CB_yield));
        }*/
        if((pdf_model == 3) && ((q2Bin == 6) || (q2Bin == 4))){
          RooExponential* exp_bkg = new RooExponential(("exp_bkg_PDF_"+year).c_str(),
                                                       ("exp_bkg_PDF_"+year).c_str(),
                                                        *mass, *lambda);

          RooGenericPdf* erf_bkg = new RooGenericPdf("erf_bkg","erf_bkg","TMath::Erfc((mass-shift)*scale)",RooArgList(*mass,*shift,*scale));

          RooAddPdf* CB_bkg = new RooAddPdf(("CB_bkg_PDF_"+year).c_str(),
                                            ("CB_bkg_PDF_"+year).c_str(),
                                             RooArgList(*exp_bkg,*erf_bkg), *cofs);
          final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
                                             ("final_PDF_extended_"+year).c_str(),
                                              RooArgList(*final_PDF, *CB_bkg), RooArgList(*signal_yield, *CB_yield));
        }
        else{
          RooExponential* CB_bkg = new RooExponential(("CB_bkg_PDF_"+year).c_str(),
                                      ("CB_bkg_PDF_"+year).c_str(),
                                       *mass, *lambda);
          final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
                                             ("final_PDF_extended_"+year).c_str(),
                                             RooArgList(*final_PDF, *CB_bkg), RooArgList(*signal_yield, *CB_yield));
        }
      }
    }
    else if(constrain == 0){ // fit without Gaussian constraints
      if(comp == 0){
        final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
                                           ("final_PDF_extended_"+year).c_str(),
                                            *dcb_rt, *signal_yield);    
      }
      else if(comp == 1){
        final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
                                           ("final_PDF_extended_"+year).c_str(),
                                            *dcb_wt, *signal_yield); 
      } 
    }

    else if(constrain == 2){ // fit with values fixed to MC
      lambda = new RooRealVar(Form("lambda^{%i}",years[iy]), "lambda", -2., -10., 1.);
      cofs = new RooRealVar(Form("cofs^{%i}",years[iy]), "cofs", 0.5, 0., 1.);

      RooRealVar* constrained_yield_RT = (RooRealVar*)fitresult_RT->floatParsFinal().find(Form("sig_yield^{%i}",years[iy]));
      RooRealVar* constrained_yield_WT = (RooRealVar*)fitresult_WT->floatParsFinal().find(Form("sig_yield^{%i}",years[iy]));

      /// create constraint on mFrac (here there is no efficiency, therefore value set to measured value on MC)
      double nrt_mc = constrained_yield_RT->getVal();
      double nwt_mc = constrained_yield_WT->getVal();
      double fraction = nwt_mc / (nrt_mc + nwt_mc);

      RooRealVar* mFrac = new RooRealVar(Form("mFrac^{%i}",years[iy]),"mistag fraction",fraction, 0, 1);
      mFrac->setConstant();

      signal_PDF = new RooAddPdf(("signal_PDF_"+year).c_str(),
                                 ("signal_PDF_"+year).c_str(),
                                  RooArgList(*dcb_wt,*dcb_rt), *mFrac);

      RooExponential* CB_bkg = new RooExponential(("CB_bkg_PDF_"+year).c_str(),
                                                  ("CB_bkg_PDF_"+year).c_str(),
                                                   *mass, *lambda);

      final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
                                         ("final_PDF_extended_"+year).c_str(),
                                          RooArgList(*signal_PDF, *CB_bkg), RooArgList(*signal_yield, *CB_yield));
    }

    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    if(multiSample) for (uint is = firstSample; is <= lastSample; is++){
        if ( !data[iy][is] || data[iy][is]->IsZombie() ){
          cout<<"Dataset " << is  << " not found in file: "<<filename_data<<endl;
          return;
        }
        map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
        simPdf->addPdf(*final_PDF_extended, ("data"+year+Form("_subs%d",is)).c_str());
    }
    else {
      if( !data[iy][0] || data[iy][0]->IsZombie() ){
        cout<<"Dataset " << firstSample  << " not found in file: "<<filename_data<<endl;
        return;
      }
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",firstSample)).c_str(), data[iy][0]) );
      simPdf->addPdf(*final_PDF_extended, ("data"+year+Form("_subs%d",firstSample)).c_str());
    }
  }

  // save initial par values into a workspace 
  // The kTRUE flag imports the values of the objects in (*params) into the workspace
  // If not set, the present values of the workspace parameters objects are stored
  ws_pars->import(*simPdf);
  RooArgSet* params = (RooArgSet*)simPdf->getParameters(*mass);
  ws_pars->saveSnapshot("initial_pars", *params, kTRUE);

  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data",
                          reco_vars,
                          Index(sample),
                          Import(map));
  RooDataSet* combData = 0;
  RooAbsReal* nll = 0;

  for(uint is = firstSample; is <= lastSample; is++){

    TFile* fout = new TFile(("simFitMassResults/simFitResult_recoMC_fullMass" + all_years + stat + Form("_b%ip%ic%im%i_subs%i", q2Bin, parity,constrain,pdf_model,is) + component + ".root").c_str(),"RECREATE");
 
    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }
    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
    if(dat == 0){
      if (nSample>0) cout<<"Fitting subsample "<<is+1<<" with "<<combData->numEntries()<<" entries"<<endl;
      else cout<<"Fitting full MC sample with "<<combData->numEntries()<<" entries"<<endl;
    }
    else if(dat == 1){
      cout << "Fitting data sample with " << combData->numEntries() << " entries" << endl;
    }
    ws_pars->loadSnapshot("initial_pars");

    TStopwatch subTime;
    if(constrain == 1){
      nll = ws_pars->pdf("simPdf")->createNLL(*combData,
                                               RooFit::Extended(kTRUE),
                                               RooFit::Constrain(c_vars),
                                               RooFit::NumCPU(1)
                                               );
    }
    else if( (constrain == 0) || (constrain == 2) ){
      nll = ws_pars->pdf("simPdf")->createNLL(*combData,
                                              RooFit::Extended(kTRUE),
                                              RooFit::NumCPU(1)
                                              );
    }

    RooMinimizer m(*nll) ;
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    m.setPrintLevel(3);//-1
    m.setPrintEvalErrors(-1);
    m.setMinimizerType("Minuit2");
    //m.setVerbose(kTRUE); 
   
    subTime.Start(true);
    m.setStrategy(0);
    m.migrad();
    m.hesse();
    m.setStrategy(2);
    m.migrad();
    m.hesse();
    subTime.Stop();
    cout << "fitting done  " << subTime.CpuTime() << endl;
    
    RooFitResult* fitResult = m.save(("result_" + shortString + "_" + Form("subs%d",is)).c_str()) ;
    fitResult->Print("v");
    double n_float_params = fitResult->floatParsFinal().getSize();    
    cout << "n_float_params = " << n_float_params << endl;

    // Save fit results in file
    if (save) {
      fout->cd();
      //fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete); //only works if root file already exists
      fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str());
      ws_pars->Write();
    }

  fout->Close();

  //if (!plot || multiSample) return;

    string plotString = shortString + "_" + all_years;
    if (nSample>0) plotString = plotString + Form("_s%isubs%i",nSample,is);

    int confIndex = 2*nBins*parity  + q2Bin;
    string longString  = Form("Fit to subsample %i - ",is) + component + " -";
    if(parity < 2){longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);}
    else{longString = longString + " both parity datasets" + Form(" - q2-bin %i ",q2Bin);}

    // plot fit projections 
    c[confIndex] = new TCanvas (("c_"+shortString+Form("subs_%i",is)).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,1400);
    c[confIndex]->cd();

    TPad *p1 = new TPad("p1","p1",0.,0.27,0.7,1.);
    p1->SetBorderMode(1);
    p1->SetFrameBorderMode(0);
    p1->SetBorderSize(2);
    p1->Draw();

    TPad *p2 = new TPad("p2","p2",0.,0.065,0.7,0.235);
    p2->SetBorderMode(1);
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);
    p2->Draw();

    TPad *p3 = new TPad("p3","p3",0.65,0.,1.,1.);
    p3->SetTopMargin(0.075);
    p3->Draw();

    RooPlot* mass_frame = mass->frame(Title("    "));
    mass_frame->GetYaxis()->SetLabelSize(0);
    mass_frame->GetYaxis()->SetTickLength(0);
    mass_frame->GetXaxis()->SetLabelSize(0);
    mass_frame->GetXaxis()->SetTickLength(0);
    mass_frame->GetXaxis()->SetTitle(" ");
    mass_frame->GetYaxis()->SetTitle(" ");

    for (unsigned int iy = 0; iy < years.size(); iy++) {
      year.clear(); year.assign(Form("%i",years[iy]));
      std::vector<RooPlot*> frames;
      frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
      TLegend *leg = new TLegend (0.65,0.7,0.8,0.88);

      for (unsigned int fr = 0; fr < frames.size(); fr++){
          // data
          combData->plotOn(frames[fr], MarkerColor(kBlack), LineColor(kBlack), Binning(80),
                           Cut(("sample==sample::data"+year+Form("_subs%d",is)).c_str()),
                           Name(("plData"+year+Form("_subs_%i",is)).c_str()));
          // total PDF
          (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineColor(kBlack), LineStyle(1), LineWidth(2),
                                           Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                         	           ProjWData(RooArgSet(sample), *combData),
                         	           Name(("plPDF"+year+Form("_subs_%i",is)).c_str()));
          if(comp == 2){
            // signal PDF
            if(constrain == 1){ 
              (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineColor(kRed+3), LineStyle(kSolid),
                                               Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                                               ProjWData(RooArgSet(sample), *combData),
                                               Components(("final_PDF_"+year).c_str()),
                                               Name(("signalPDF"+year+Form("_subs_%i",is)).c_str()));
            }
            else if(constrain == 2){
              (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineColor(kRed+3), LineStyle(kSolid),
                                               Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                                               ProjWData(RooArgSet(sample), *combData),
                                               Components(("signal_PDF_"+year).c_str()),
                                               Name(("signalPDF"+year+Form("_subs_%i",is)).c_str()));
           }
           // RT PDF
           (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineColor(8), LineStyle(kSolid),
                                            FillStyle(3002), FillColor(8), VLines(), DrawOption("F"),
                                            Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                                            ProjWData(RooArgSet(sample), *combData),
                                            Components(("dcb_rt_"+year).c_str()),
                                            Name(("RT_PDF"+year+Form("_subs_%i",is)).c_str()));
           // WT PDF
           (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineColor(2), LineStyle(kSolid),
                                            FillStyle(3002), FillColor(2), VLines(), DrawOption("F"),
                                            Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                                            ProjWData(RooArgSet(sample), *combData),
                                            Components(("dcb_wt_"+year).c_str()),
                                            Name(("WT_PDF"+year+Form("_subs_%i",is)).c_str()));
           // CB bkg PDF
           if(dat == 1){
             (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineColor(kBlue), LineStyle(kSolid), 
                                              Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                                              ProjWData(RooArgSet(sample), *combData),
                                              Components(("CB_bkg_PDF_"+year).c_str()),
                                              Name(("CB_bkgPDF"+year+Form("_subs_%i",is)).c_str()));
  
           }
          }

          combData->plotOn(mass_frame,MarkerColor(kWhite),LineColor(kWhite));
          (ws_pars->pdf("simPdf"))->plotOn(mass_frame,ProjWData(RooArgSet(sample), *combData), LineColor(kWhite));
          (ws_pars->pdf("simPdf"))->paramOn(mass_frame,Layout(0.1,0.9,1.));

          double chis = frames[fr]->chiSquare(("plPDF"+year+Form("_subs_%i",is)).c_str(), ("plData"+year+Form("_subs_%i",is)).c_str(),n_float_params);

          TLatex* tex = new TLatex(0.15, 0.6, Form("#chi^{2}/ndf = %f",chis));//%.3lf
          tex->SetNDC(kTRUE);
          tex->SetTextFont(42);
          tex->SetTextSize(0.04);

          TLatex* tex1;
          if(q2Bin == 0){tex1 = new TLatex(0.15, 0.8, "1 < q^{2} < 2 GeV^{2}");}
          else if(q2Bin == 1){tex1 = new TLatex(0.15, 0.8, "2 < q^{2} < 4.3 GeV^{2}");}
          else if(q2Bin == 2){tex1 = new TLatex(0.15, 0.8, "4.3 < q^{2} < 6 GeV^{2}");}
          else if(q2Bin == 3){tex1 = new TLatex(0.15, 0.8, "6 < q^{2} < 8.68 GeV^{2}");}
          else if(q2Bin == 4){tex1 = new TLatex(0.15, 0.8, "8.68 < q^{2} < 10.09 GeV^{2} (J/ #psi)");}
          else if(q2Bin == 5){tex1 = new TLatex(0.15, 0.8, "10.09 < q^{2} < 12.86 GeV^{2}");}
          else if(q2Bin == 6){tex1 = new TLatex(0.15, 0.8, "12.86 < q^{2} < 14.18 GeV^{2} (#psi (2S))");}
          else if(q2Bin == 7){tex1 = new TLatex(0.15, 0.8, "14.18 < q^{2} < 16 GeV^{2}");}
          tex1->SetNDC(kTRUE);
          tex1->SetTextFont(42);
          tex1->SetTextSize(0.04);

          RooRealVar* yield = (RooRealVar*)fitResult->floatParsFinal().find(Form("sig_yield^{%i}",years[iy]));
          TLatex* tex2;
          if(q2Bin == 4){tex2 = new TLatex(0.15, 0.7, Form("Y_{N} = %.0lf #pm %.0lf",yield->getVal(),yield->getError()));}
          else if(q2Bin == 6){tex2 = new TLatex(0.15, 0.7, Form("Y = %.0lf #pm %.0lf",yield->getVal(),yield->getError()));}
          else{tex2 = new TLatex(0.15, 0.7, Form("Y_{S} = %.0lf #pm %.0lf",yield->getVal(),yield->getError()));}
          tex2->SetNDC(kTRUE);
          tex2->SetTextFont(42);
          tex2->SetTextSize(0.04);

          if (fr == 0) {
            leg->AddEntry(frames[fr]->findObject(("plData"+year+Form("_subs_%i",is)).c_str()), "Data" ,"lep");
            leg->AddEntry(frames[fr]->findObject(("plPDF"+year+Form("_subs_%i",is)).c_str()), "Fit","l");
            if(comp == 2){
              leg->AddEntry(frames[fr]->findObject(("signalPDF"+year+Form("_subs_%i",is)).c_str()), "Signal","l");
              leg->AddEntry(frames[fr]->findObject(("RT_PDF"+year+Form("_subs_%i",is)).c_str()), "RT component","l");
              leg->AddEntry(frames[fr]->findObject(("WT_PDF"+year+Form("_subs_%i",is)).c_str()), "WT component","l");
              if(dat == 1){
                leg->AddEntry(frames[fr]->findObject(("CB_bkgPDF"+year+Form("_subs_%i",is)).c_str()), "Background","l");
              }
            }
          }

          //c[confIndex]->cd(iy+1);

          p1->cd();
          frames[fr]->SetTitle(Form("CMS Preliminary            L = 139.5 fb^{-1}             #sqrt{s} = 13 TeV (pp)            Year %i",years[iy]));
          frames[fr]->SetXTitle("m(K^{+} #pi^{-} #mu^{+} #mu^{-}) (GeV)");
          frames[fr]->SetYTitle(TString::Format("Events / (%g)",(mass->getMax()-mass->getMin())/80));
          frames[fr]->GetYaxis()->SetTitleFont(43);
          frames[fr]->GetYaxis()->SetTitleSize(30);
          frames[fr]->GetYaxis()->SetTitleOffset(1.6);
          frames[fr]->GetYaxis()->SetLabelFont(43);
          frames[fr]->GetYaxis()->SetLabelSize(30);
          frames[fr]->Draw();
          leg->SetTextSize(0.04);
          leg->SetBorderSize(0);
          leg->Draw("same");
          tex->Draw("same");     
          tex1->Draw("same");
          tex2->Draw("same");

          // ratio plot
          RooHist* pull_hist = frames[fr]->pullHist(("plData"+year+Form("_subs_%i",is)).c_str(),("plPDF"+year+Form("_subs_%i",is)).c_str());
          RooPlot *pull_plot = mass->frame();

          pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
          pull_plot->SetTitle("");
          pull_plot->GetXaxis()->SetTitle("m(K^{+} #pi^{-} #mu^{+} #mu^{-}) (GeV)");
          pull_plot->GetXaxis()->SetTitleSize(0.15);
          pull_plot->GetXaxis()->SetTitleOffset(0.9);
          pull_plot->GetXaxis()->SetLabelSize(0.15);
          pull_plot->GetXaxis()->SetTickLength(0.1);
          pull_plot->GetYaxis()->SetTitle("Pull");
          pull_plot->GetYaxis()->SetTitleSize(0.13);
          pull_plot->GetYaxis()->SetTitleOffset(0.13);
          pull_plot->GetYaxis()->SetLabelSize(0.1);
          pull_plot->GetYaxis()->SetNdivisions(305);

          gPad->Update();
          TLine *line = new TLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0); 
          line->SetLineStyle(2);
          line->SetLineColor(kBlue);
   
          p2->cd();
          pull_plot->Draw();
          line->Draw("same");

          p3->cd();
          mass_frame->Draw();

          break;
      }// ends loop over frames
    }// ends loop over years

    if(dat == 0){
      if(nSample > 0){
        if(multiSample){
          c[confIndex]->SaveAs( ("plotSimMassFit_dataStat_multiSample/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
          c[confIndex]->SaveAs( ("plotSimMassFit_dataStat_multiSample/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".pdf").c_str());
        }
        else{
          c[confIndex]->SaveAs( ("plotSimMassFit_dataStat/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
          c[confIndex]->SaveAs( ("plotSimMassFit_dataStat/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".pdf").c_str());
        }
      }   
      else{
        if(parity < 2){
          if(constrain == 1){
            c[confIndex]->SaveAs( ("plotSimMassFit_constrain/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
            c[confIndex]->SaveAs( ("plotSimMassFit_constrain/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".pdf").c_str());
          }
          else if(constrain == 0){ 
            c[confIndex]->SaveAs( ("plotSimMassFit/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str()); 
            c[confIndex]->SaveAs( ("plotSimMassFit/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".pdf").c_str());
          }
        }
        else{
          if(constrain == 1){
            c[confIndex]->SaveAs( ("plotSimMassFit_odd+even_constrain/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
            c[confIndex]->SaveAs( ("plotSimMassFit_odd+even_constrain/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".pdf").c_str());
          }
          else if(constrain == 0){
            c[confIndex]->SaveAs( ("plotSimMassFit_odd+even/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
            c[confIndex]->SaveAs( ("plotSimMassFit_odd+even/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".pdf").c_str());
          }
        }
      }
    }
    else if(dat == 1){
      c[confIndex]->SaveAs( ("plotSimMassFit_DATA/simFitResult_recoDATA_fullMass_" + plotString + "_" + component + ".gif").c_str());
      c[confIndex]->SaveAs( ("plotSimMassFit_DATA/simFitResult_recoDATA_fullMass_" + plotString + "_" + component + ".pdf").c_str());
    }

    // validate fit
    //if( (nSample == 1) || ((dat == 1) && (pdf_model == 0) && (constrain < 2)) ){
    //  validate_fit(ws_pars, sample, c_vars, years[0], q2Bin, parity, dat, constrain, nSample);
    //}
  }// ends loop over samples 
}

void simfit_recoMC_fullMassBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, int constrain, int comp, int dat, int pdf_model, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, dat, pdf_model, years);
  else
    simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, dat, pdf_model, years);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficienc
  //                [1] odd efficiency
  //                [2] for merged even and odd datasets
  //                [-1] for each parity recursively
 
  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);
  if ( argc > 2 ) parity  = atoi(argv[2]);

  bool multiSample = false;
  if ( argc > 3 && atoi(argv[3]) > 0 ) multiSample = true;

  uint nSample = 0;
  if ( argc > 4 ) nSample = atoi(argv[4]);

  if (nSample==0) multiSample = false;

  bool plot = true;
  bool save = true;

  if ( argc > 5 && atoi(argv[5]) == 0 ) plot = false;
  if ( argc > 6 && atoi(argv[6]) == 0 ) save = false;

  int constrain = 0;
  if ( argc > 7 ) constrain = atoi(argv[7]);

  int comp = 0;
  if ( argc > 8 ) comp = atoi(argv[8]);

  int dat = 0;
  if ( argc > 9 ) dat = atoi(argv[9]);

  int pdf_model = 0;
  if ( argc > 10 ) pdf_model = atoi(argv[10]);

  std::vector<int> years;
  if ( argc > 11 && atoi(argv[11]) != 0 ) years.push_back(atoi(argv[11]) );
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 12 && atoi(argv[12]) != 0 ) years.push_back(atoi(argv[12]));
  if ( argc > 13 && atoi(argv[13]) != 0 ) years.push_back(atoi(argv[13]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "plot        " << plot         << endl;
  cout <<  "save        " << save         << endl;
  cout <<  "constrain   " << constrain    << endl;
  cout <<  "comp        " << comp         << endl;
  cout <<  "dat         " << dat          << endl;
  cout <<  "pdf_model   " << pdf_model    << endl;
  cout <<  "years[0]    " << years[0]     << endl;
  //cout <<  "years[1]    " << years[1]     << endl;
  //cout <<  "years[2]    " << years[2]     << endl;

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 2      ) return 1;
  if ( (constrain == 1) & (comp < 2) ) return 1; //constrained fit (for RT+WT full MC, data-like MC and data)
  if ( (constrain == 0) & (comp >=2) ) return 1; //unconstrained fit (for RT MC and WT MC)
  if ( (constrain == 2) & (dat==0) ) return 1; //fit with parameters fixed to MC values and scale factor (for data)
  if ( (dat == 1) & (comp < 2) ) return 1;
  if ( (dat == 1) & (nSample > 0) ) return 1;
  if ( (dat == 1) & (multiSample == true) ) return 1;   

  if ( q2Bin == -1 )   cout << "Running all the q2 bins" << endl;
  if ( parity == -1 )  cout << "Running both the parity datasets recursively" << endl;
  else if (parity == 2) cout << "Running both parity datasets merged" << endl;

  if (dat == 0) cout << "Running MC" << endl;
  else if (dat == 1) cout << "Running data" << endl;

  // https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
  if(nSample > 0){
    if(parity < 2){   
      scale_to_data.insert(std::make_pair(2016, 0.006*2 / 2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
      scale_to_data.insert(std::make_pair(2017, 0.005*2 / 2.05 ));
      scale_to_data.insert(std::make_pair(2018, 0.007*2 / 1.9  ));
    }
    else{
      scale_to_data.insert(std::make_pair(2016, 0.006 / 2.5  )); 
      scale_to_data.insert(std::make_pair(2017, 0.005 / 2.05 ));
      scale_to_data.insert(std::make_pair(2018, 0.007 / 1.9  ));
    }
  }

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, dat, pdf_model, years);
  else
    simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, dat, pdf_model, years);

  return 0;

}
