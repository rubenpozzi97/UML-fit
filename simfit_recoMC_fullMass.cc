#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>

#include <RooArgList.h>
#include <RooAbsPdf.h>
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

void simfit_recoMC_fullMassBin(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, bool constrain, int comp, std::vector<int> years)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%ic%i",q2Bin,parity,constrain);
  cout<<"Conf: "<<shortString<<endl;

  string all_years = "";
  string year = ""; 
  string isample = ""; 
  string stat = nSample > 0 ? "_dataStat":"_MCStat";
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;

  std::vector<RooWorkspace*> wsp, wsp_mcmass, wsp_even, wsp_odd;
  std::vector<std::vector<RooDataSet*>> data;

  std::vector<RooAbsPdf*> PDF_sig_mass(0);
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
  else{
    mass = new RooRealVar("mass","mass", 5.,5.6);
  }
  ws_pars->import(*mass);

  RooRealVar* rand = new RooRealVar("rand","rand", 0,1);
  RooArgSet reco_vars (*mass, *rand);

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    all_years += year;
    for (uint is = firstSample; is <= lastSample; is++) {
      isample.clear(); isample.assign( Form("%i",is) );
      sample.defineType(("data"+year+"_subs"+isample).c_str());
    }
  }

  // Construct a simultaneous pdf using category sample as index (simultaneous over the 3 years)
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("recoMCDataset_b%i_%i.root", q2Bin, years[iy]);
    filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/newphi/", years[iy]) + filename_data;

    // import data (or MC as data proxy)
    if(parity < 2){if (!retrieveWorkspace(filename_data, wsp, Form("ws_b%ip%i", q2Bin, 1-parity ), wsp, Form("ws_b%ip%i", q2Bin, 1-parity ), parity))  return;}
    else{if (!retrieveWorkspace(filename_data, wsp_even, Form("ws_b%ip%i", q2Bin, 0), wsp_odd, Form("ws_b%ip%i", q2Bin, 1), parity))  return;}

    if(parity < 2){
      data.push_back(createDataset(nSample,  firstSample,  lastSample, wsp[iy], wsp[iy], 
                                   q2Bin,  parity,  years[iy], 
                                   reco_vars,  shortString, comp));
    }
    else{
      data.push_back(createDataset(nSample,  firstSample,  lastSample, wsp_even[iy], wsp_odd[iy],
                                   q2Bin,  parity,  years[iy],
                                   reco_vars,  shortString, comp));
    }
    
    // import mass PDF from fits to the MC 
    string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i_fM_newbdt.root",years[iy]);
    if (!retrieveWorkspace(filename_mc_mass, wsp_mcmass, "w", wsp_mcmass, "w", 0))  return;

    RooRealVar* mean_rt = 0;
    RooRealVar* sigma_rt = 0;
    RooRealVar* alpha_rt1 = 0;
    RooRealVar* alpha_rt2 = 0;
    RooRealVar* n_rt1 = 0;
    RooRealVar* n_rt2 = 0;
    RooRealVar* sigma_rt2 = 0;
    RooRealVar* f1rt = 0;
    RooRealVar* mean_wt = 0;
    RooRealVar* sigma_wt = 0;
    RooRealVar* alpha_wt1 = 0;
    RooRealVar* alpha_wt2 = 0;
    RooRealVar* n_wt1 = 0;
    RooRealVar* n_wt2 = 0;

    RooAbsPdf* dcb_rt;
    RooProdPdf* c_dcb_rt;
    RooAbsPdf* dcb_wt;
    RooProdPdf* c_dcb_wt;
    RooProdPdf* final_PDF;
    RooAddPdf* final_PDF_extended;
    RooRealVar* signal_yield;

    RooRealVar* mFrac = new RooRealVar(Form("mFrac^{%i}",years[iy]),"mistag fraction",0.13, 0, 1);
    // mFrac->setConstant();

    int signal_yield_initial = data[iy][0]->numEntries();
    signal_yield = new RooRealVar(Form("sig_yield^{%i}",years[iy]), "signal_yield", signal_yield_initial, 0, 2*signal_yield_initial);
    cout << "signal yield in data = " << signal_yield_initial << endl;

    // create RT component 
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
    mean_rt       = new RooRealVar (Form("mean_{RT}^{%i}",years[iy])    , "massrt"      , wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal()     ,      5,    6, "GeV");
    sigma_rt      = new RooRealVar (Form("#sigma_{RT1}^{%i}",years[iy] ), "sigmart1"    , wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal()  ,      0,    1, "GeV");
    alpha_rt1     = new RooRealVar (Form("#alpha_{RT1}^{%i}",years[iy] ), "alphart1"    , wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    alpha_rt2     = new RooRealVar (Form("#alpha_{RT2}^{%i}",years[iy] ), "alphart2"    , wsp_mcmass[iy]->var(Form("#alpha_{RT2}^{%i}", q2Bin))->getVal() ,    -10,   10 );
    n_rt1         = new RooRealVar (Form("n_{RT1}^{%i}",years[iy])      , "nrt1"        , wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal()      ,      0.,  200.);
    n_rt2         = new RooRealVar (Form("n_{RT2}^{%i}",years[iy])      , "nrt2"        , wsp_mcmass[iy]->var(Form("n_{RT2}^{%i}", q2Bin))->getVal()      ,      0.,  200.);      

    ws_pars->import(*mean_rt);
    ws_pars->import(*sigma_rt);
    ws_pars->import(*alpha_rt1);
    ws_pars->import(*alpha_rt2);
    ws_pars->import(*n_rt1);
    ws_pars->import(*n_rt2);

    if (q2Bin >= 5){
      sigma_rt2 = new RooRealVar (Form("#sigma_{RT2}^{%i}",years[iy] ), "sigmaRT2"  ,    wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal(), 0,   0.12, "GeV");
      f1rt      = new RooRealVar (Form("f^{RT%i}",years[iy])          , "f1rt"      ,   wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal(), 0,  1.);

      ws_pars->import(*sigma_rt2);
      ws_pars->import(*f1rt);

      alpha_rt2->setRange(-10,0);
     
      if(constrain){ 
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt);
      } 
      else{
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy], false, c_vars_rt, c_pdfs_rt);
      }
    }
    else{
      alpha_rt2->setRange(0,10);
       
      if(constrain){  
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2, wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt);
      }
      else{
        dcb_rt = createRTMassShape(q2Bin, mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2, wsp_mcmass[iy], years[iy], false, c_vars_rt, c_pdfs_rt);
      }
    }

    // create WT component
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));
    mean_wt     = new RooRealVar (Form("mean_{WT}^{%i}",years[iy])      , "masswt"     ,  wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal()    ,      5,    6, "GeV");
    sigma_wt    = new RooRealVar (Form("#sigma_{WT1}^{%i}",years[iy])   , "sigmawt"    ,  wsp_mcmass[iy]->var(Form("#sigma_{WT1}^{%i}", q2Bin))->getVal() ,      0,    1, "GeV");
    alpha_wt1   = new RooRealVar (Form("#alpha_{WT1}^{%i}",years[iy] )  , "alphawt1"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    alpha_wt2   = new RooRealVar (Form("#alpha_{WT2}^{%i}",years[iy] )  , "alphawt2"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT2}^{%i}", q2Bin))->getVal() ,      0,   10 );
    n_wt1       = new RooRealVar (Form("n_{WT1}^{%i}",years[iy])        , "nwt1"       ,  wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal()      ,      0., 100.);
    n_wt2       = new RooRealVar (Form("n_{WT2}^{%i}",years[iy])        , "nwt2"       ,  wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal()      ,      0., 100.);

    ws_pars->import(*mean_wt);
    ws_pars->import(*sigma_wt);
    ws_pars->import(*alpha_wt1);
    ws_pars->import(*alpha_wt2);
    ws_pars->import(*n_wt1);
    ws_pars->import(*n_wt2);

    if(constrain){
      dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2, wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt );
    }
    else{
      dcb_wt = createWTMassShape(q2Bin, mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2, wsp_mcmass[iy], years[iy], false, c_vars_wt, c_pdfs_wt );
    }

    if(constrain){
      /// create constrained PDF for RT mass
      RooArgList constr_rt_list = RooArgList(c_pdfs_rt);
      constr_rt_list.add(*dcb_rt);
      c_dcb_rt = new RooProdPdf(("c_dcb_rt_"+year).c_str(), ("c_dcb_rt_"+year).c_str(), constr_rt_list );
      c_vars.add(c_vars_rt);

      /// create constrained PDF for WT mass   
      RooArgList constr_wt_list = RooArgList(c_pdfs_wt);
      constr_wt_list.add(*dcb_wt);
      c_dcb_wt = new RooProdPdf(("c_dcb_wt_"+year).c_str(), ("c_dcb_wt_"+year).c_str(), constr_wt_list );
      c_vars.add(c_vars_wt);

      /// create constraint on mFrac (here there is no efficiency, therefore value set to measured value on MC)
      double nrt_mc   =  wsp_mcmass[iy]->var(Form("nRT_%i",q2Bin))->getVal();
      double nwt_mc   =  wsp_mcmass[iy]->var(Form("nWT_%i",q2Bin))->getVal();
      double fraction = nwt_mc / nrt_mc;

      c_fm.push_back(new RooGaussian(Form("c_fm^{%i}",years[iy]) , "c_fm" , *mFrac,
                                      RooConst(fraction),
                                      RooConst(fM_sigmas[years[iy]][q2Bin])
                                      ));

      cout << fraction << "   " << fM_sigmas[years[iy]][q2Bin] << endl;
      c_vars.add(*mFrac);

      if (q2Bin < 5){
          PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                                 ("PDF_sig_mass_"+year).c_str(),
                                                 *mass,
                                                 *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
                                                 *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
                                                 *mFrac, 
                                                 *c_dcb_rt,
                                                 *c_dcb_wt
                                                 ));
      }
      else{
          PDF_sig_mass.push_back( new PdfSigMass(("PDF_sig_mass_"+shortString+"_"+year).c_str(),
                                                 ("PDF_sig_mass_"+year).c_str(),
                                                 *mass,
                                                 *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
                                                 *mean_wt, *sigma_wt,             *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
                                                 *mFrac, 
                                                 *c_dcb_rt,
                                                 *c_dcb_wt
                                                 ));
      }

      final_PDF = new RooProdPdf(("final_PDF_"+year).c_str(),
                                 ("final_PDF_"+year).c_str(),
                                 RooArgList(*PDF_sig_mass[iy], *c_fm[iy]));                         

      final_PDF_extended = new RooAddPdf(("final_PDF_extended_"+year).c_str(),
					 ("final_PDF_extended_"+year).c_str(),
					  *final_PDF, *signal_yield);

    }

    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    if (multiSample) for (uint is = firstSample; is <= lastSample; is++) {
        if ( !data[iy][is] || data[iy][is]->IsZombie() ) {
          cout<<"Dataset " << is  << " not found in file: "<<filename_data<<endl;
          return;
        }
        map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
        if(comp == 0){simPdf->addPdf(*dcb_rt, ("data"+year+Form("_subs%d",is)).c_str());}
        else if(comp == 1){simPdf->addPdf(*dcb_wt, ("data"+year+Form("_subs%d",is)).c_str());}
        else{simPdf->addPdf(*final_PDF_extended, ("data"+year+Form("_subs%d",is)).c_str());}
    }
    else {
      if ( !data[iy][0] || data[iy][0]->IsZombie() ) {
        cout<<"Dataset " << firstSample  << " not found in file: "<<filename_data<<endl;
        return;
      }
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",firstSample)).c_str(), data[iy][0]) );
      if(comp == 0){simPdf->addPdf(*dcb_rt, ("data"+year+Form("_subs%d",firstSample)).c_str());} 
      else if(comp == 1){simPdf->addPdf(*dcb_wt, ("data"+year+Form("_subs%d",firstSample)).c_str());}
      else{simPdf->addPdf(*final_PDF_extended, ("data"+year+Form("_subs%d",firstSample)).c_str());}
    }
  }

  string component = "";
  if(comp == 0){component = "CT";}
  else if(comp == 1){component = "WT";}
  else{component = "CT+WT";}

  // save initial par values into a workspace 
  // The kTRUE flag imports the values of the objects in (*params) into the workspace
  // If not set, the present values of the workspace parameters objects are stored
  ws_pars->import(*simPdf);
  RooArgSet *params = (RooArgSet*)simPdf->getParameters(*mass);
  ws_pars->saveSnapshot("initial_pars", *params, kTRUE);

  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data",
                          reco_vars,
                          Index(sample),
                          Import(map));
  RooDataSet* combData = 0;
  RooAbsReal* nll = 0;

  for (uint is = firstSample; is <= lastSample; is++) {

    TFile* fout = new TFile(("simFitMassResults/simFitResult_recoMC_fullMass" + all_years + stat + Form("_b%ip%ic%i_subs%i", q2Bin, parity,constrain,is) + component + ".root").c_str(),"UPDATE");
 
    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }
    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
    if (nSample>0) cout<<"Fitting subsample "<<is+1<<" with "<<combData->numEntries()<<" entries"<<endl;
    else cout<<"Fitting full MC sample with "<<combData->numEntries()<<" entries"<<endl;

    ws_pars->loadSnapshot("initial_pars");

    TStopwatch subTime;
    if(constrain){
      nll = ws_pars->pdf("simPdf")->createNLL(*combData,
                                              RooFit::Extended(kTRUE),
                                              RooFit::Constrain(c_vars),
                                              RooFit::NumCPU(1)
                                              );
    }
    else{
      nll = ws_pars->pdf("simPdf")->createNLL(*combData,
                                              RooFit::Extended(kFALSE),
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
  
    // Save fit results in file
    if (save) {
       fout->cd();
       //fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete); //only works if root file already exists
       fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str()); 
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
    c[confIndex] = new TCanvas (("c_"+shortString+Form("subs_%i",is)).c_str(),("Fit to RECO-level MC - "+longString).c_str(),3000,1400);
    c[confIndex]->Divide(years.size());
    cout<<"plotting 4d canvas: "<< years.size() << endl;

    c[confIndex]->cd();

    TPad *p1 = new TPad("p1","p1",0.,0.27,1.,1.);
    p1->SetBorderMode(1);
    p1->SetFrameBorderMode(0);
    p1->SetBorderSize(2);
    //p1->SetBottomMargin(0.40);
    p1->Draw();

    TPad *p2 = new TPad("p2","p2",0.,0.065,1.,0.24);
    p2->SetBorderMode(1);
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);
    p2->Draw();

    for (unsigned int iy = 0; iy < years.size(); iy++) {
      year.clear(); year.assign(Form("%i",years[iy]));
      std::vector<RooPlot*> frames;
      frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
      TLegend *leg = new TLegend (0.6,0.7,0.9,0.9);

      for (unsigned int fr = 0; fr < frames.size(); fr++){
          combData->plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(80),
                           Cut(("sample==sample::data"+year+Form("_subs%d",is)).c_str()),
                           Name(("plData"+year+Form("_subs_%i",is)).c_str()));
          (ws_pars->pdf("simPdf"))->plotOn(frames[fr], LineWidth(1),
                        Slice(sample, ("data"+year+Form("_subs%i",is)).c_str()),
                        ProjWData(RooArgSet(sample), *combData),
                        Name(("plPDF"+year+Form("_subs_%i",is)).c_str()));
          (ws_pars->pdf("simPdf"))->paramOn(frames[fr],Layout(0.65,0.9,0.95));

          if (fr == 0) {
            leg->AddEntry(frames[fr]->findObject(("plData"+year+Form("_subs_%i",is)).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
            leg->AddEntry(frames[fr]->findObject(("plPDF"+year+Form("_subs_%i",is)).c_str()),("Mass component of the pdf "+year).c_str(),"l");
          }

          c[confIndex]->cd(iy+1);

          p1->cd();
          frames[fr]->SetXTitle("mass (GeV)");
          frames[fr]->GetYaxis()->SetTitleFont(43);
          frames[fr]->GetYaxis()->SetTitleSize(30);
          frames[fr]->GetYaxis()->SetTitleOffset(1.);
          frames[fr]->GetYaxis()->SetLabelFont(43);
          frames[fr]->GetYaxis()->SetLabelSize(30);
          frames[fr]->Draw();
          //leg->Draw("same");
     
          // ratio plot
          RooHist* pull_hist = frames[fr]->pullHist(("plData"+year+Form("_subs_%i",is)).c_str(),("plPDF"+year+Form("_subs_%i",is)).c_str());
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

          break;
      }// ends loop over frames
    }// ends loop over years

    if(nSample > 0){
      if(multiSample){
        c[confIndex]->SaveAs( ("plotSimMassFit_dataStat_multiSample/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
      }
      else{
        c[confIndex]->SaveAs( ("plotSimMassFit_dataStat/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
      }
    }   
    else{
      if(parity < 2){
        if(constrain){
          c[confIndex]->SaveAs( ("plotSimMassFit_constrain/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
        }
        else{ 
          c[confIndex]->SaveAs( ("plotSimMassFit/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str()); 
        }
      }
      else{
        if(constrain){
          c[confIndex]->SaveAs( ("plotSimMassFit_odd+even_constrain/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
        }
        else{
          c[confIndex]->SaveAs( ("plotSimMassFit_odd+even/simFitResult_recoMC_fullMass_" + plotString + "_" + component + ".gif").c_str());
        }
      }
    }

  // validate fit
  //validate_fit(ws_pars, sample, years[0], q2Bin, parity);
  }// ends loop over samples 
}

void simfit_recoMC_fullMassBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, bool constrain, int comp, std::vector<int> years)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, years);
  else
    simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, years);
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

  bool constrain = true;
  if ( argc > 7 && atoi(argv[7]) == 0 ) constrain = false;

  int comp = 0;
  if ( argc > 8 ) comp = atoi(argv[8]);

  std::vector<int> years;
  if ( argc > 9 && atoi(argv[9]) != 0 ) years.push_back(atoi(argv[9]) );
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 10 && atoi(argv[10]) != 0 ) years.push_back(atoi(argv[10]));
  if ( argc > 11 && atoi(argv[11]) != 0 ) years.push_back(atoi(argv[11]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "plot        " << plot         << endl;
  cout <<  "save        " << save         << endl;
  cout <<  "constrain   " << constrain    << endl;
  cout <<  "comp        " << comp         << endl;
  cout <<  "years[0]    " << years[0]     << endl;
  //cout <<  "years[1]    " << years[1]     << endl;
  //cout <<  "years[2]    " << years[2]     << endl;

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 2      ) return 1;
  if ( (constrain == true) & (comp < 2) ) return 1;
  if ( (constrain == false) & (comp >=2)) return 1;

  if ( q2Bin == -1 )   cout << "Running all the q2 bins" << endl;
  if ( parity == -1 )  cout << "Running both the parity datasets recursively" << endl;
  else if (parity == 2) cout << "Running both parity datasets merged" << endl;

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
      simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, years);
  else
    simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, constrain, comp, years);

  return 0;

}
