#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <list>
#include <map>

// #include <RooRealVar.h>
#include <RooAbsPdf.h>
// #include <RooWorkspace.h>
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
// #include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooConstVar.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"


// #include "PdfCBShape.h"
#include "PdfSigMass.h"
#include "BoundCheck.h"
#include "BoundDist.h"
#include "Penalty.h"
#include "utils.h"
#include "PdfSigRTMass.h"
#include "PdfSigWTMass.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* cPen;
TCanvas* c [4*nBins];

double power = 1.0;

double maxCoeff = 1e8;

double min_base = 1.05;

void simfit_recoMC_fullMassBin(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, std::vector<int> years, std::map<int,float> scale_to_data, double fac1, double fac4, double fac5, double base1, double base4, double base5, double max1, double max4, double max5)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency

  string all_years = "";
  string year = ""; 
  string isample = ""; 
  string stat = nSample > 0 ? "_dataStat":"_MCStat";
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;
  
  std::vector<RooWorkspace*> wsp, wsp_mcmass;
  std::vector<std::vector<RooDataSet*>> data;
  std::vector<RooAbsReal*> effC, effW;
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular (0);
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular_penalty (0);
  std::vector<RooAbsPdf*> PDF_sig_mass(0);
  std::vector<RooAbsPdf*> PDF_sig_ang_mass(0);
//   std::vector<RooAbsPdf*> PDF_sig_ang_mass_penalty(0);
  std::vector<RooGaussian*> c_sigma_rt, c_sigma_rt2, c_mean_rt, c_alpha_rt1, c_alpha_rt2, c_n_rt1, c_n_rt2, c_f1rt;
  std::vector<RooGaussian*> c_sigma_wt,              c_mean_wt, c_alpha_wt1, c_alpha_wt2, c_n_wt1, c_n_wt2;
  std::vector<RooGaussian*> c_deltaPeaks, c_fm;
  RooArgSet c_vars_rt, c_pdfs_rt;
  RooArgSet c_vars_wt, c_pdfs_wt;
  RooArgSet c_vars; 

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

//   RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
//   RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
//   RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
//   RooArgList vars (* ctK,* ctL,* phi);
  RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6);
  RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
  RooRealVar* genSignal = new RooRealVar("genSignal", "genSignal", 0,10);
  RooRealVar* tagB0     = new RooRealVar("tagB0", "tagB0", 0,10);
  RooArgSet reco_vars ( *mass, *rand, *genSignal, *tagB0);
//   RooArgSet reco_vars (*ctK, *ctL, *phi, *mass, *rand);

  // define angular parameters with ranges from positiveness requirements on the decay rate
//   RooRealVar* Fl    = new RooRealVar("Fl","F_{L}",0.5,0,1);
//   RooRealVar* P1    = new RooRealVar("P1","P_{1}",0,-1,1);   
//   RooRealVar* P2    = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
//   RooRealVar* P3    = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
//   RooRealVar* P4p   = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
//   RooRealVar* P5p   = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
//   RooRealVar* P6p   = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
//   RooRealVar* P8p   = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* mFrac = new RooRealVar("mFrac","mistag fraction",0.87, 0, 2);
  mFrac->setConstant();

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    all_years += year;
    for (uint is = firstSample; is <= lastSample; is++) {
      isample.clear(); isample.assign( Form("%i",is) );
      sample.defineType(("data"+year+"_subs"+isample).c_str());
    }
  } 

  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/recoMCDataset_b%i_%i.root", years[iy], q2Bin, years[iy]); 

    // import data (or MC as data proxy)
    if (!retrieveWorkspace( filename_data, wsp, Form("ws_b%ip%i", q2Bin, 1-parity )))  return;

    // create roodataset (in case data-like option is selected, only import the correct % of data)
    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> data_isample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {
        dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[years[iy]], (is+1)*scale_to_data[years[iy]] )) ;
        dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[years[iy]], (is+1)*scale_to_data[years[iy]] )) ;

        RooDataSet* datatmp = new RooDataSet(*dataCT,("data_"+shortString + Form("_subs%i", is)).c_str());
        datatmp->append(*dataWT);
        datatmp->Print();
        data_isample.push_back (datatmp);
      }
    }
    else{
      dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
    
      RooDataSet* datatmp = new RooDataSet(*dataCT,("data_"+shortString + "_subs0").c_str());
      datatmp->append(*dataWT);
      data_isample.push_back (datatmp);
    }

    
    data.push_back(data_isample) ;

    // Mass Component
    // import mass PDF from fits to the MC
    string filename_mc_mass = Form("/eos/cms/store/user/fiorendi/p5prime/massFits/results_fits_%i.root",years[iy]);
    if (!retrieveWorkspace( filename_mc_mass, wsp_mcmass, "w"))  return;

    RooRealVar* mean_rt       = new RooRealVar (Form("mean_{RT}^{%i}",years[iy])    , "massrt"      , wsp_mcmass[iy]->var(Form("mean_{RT}^{%i}",q2Bin))->getVal()     ,      5,    6, "GeV");
    RooRealVar* sigma_rt      = new RooRealVar (Form("#sigma_{RT1}^{%i}",years[iy] ), "sigmart1"    , wsp_mcmass[iy]->var(Form("#sigma_{RT1}^{%i}",q2Bin))->getVal()  ,      0,    1, "GeV");
    RooRealVar* alpha_rt1     = new RooRealVar (Form("#alpha_{RT1}^{%i}",years[iy] ), "alphart1"    , wsp_mcmass[iy]->var(Form("#alpha_{RT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* alpha_rt2     = new RooRealVar (Form("#alpha_{RT2}^{%i}",years[iy] ), "alphart2"    , wsp_mcmass[iy]->var(Form("#alpha_{RT2}^{%i}", q2Bin))->getVal() ,    -10,   10 );
    RooRealVar* n_rt1         = new RooRealVar (Form("n_{RT1}^{%i}",years[iy])      , "nrt1"        , wsp_mcmass[iy]->var(Form("n_{RT1}^{%i}", q2Bin))->getVal()      ,      0.,  100.);
    RooRealVar* n_rt2         = new RooRealVar (Form("n_{RT2}^{%i}",years[iy])      , "nrt2"        , wsp_mcmass[iy]->var(Form("n_{RT2}^{%i}", q2Bin))->getVal()      ,      0.,  100.);
    /// create constrain RT and add them to list of constraining pdf and vars
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_RT_%i",q2Bin));
//     c_sigma_rt.push_back(  constrainVar(sigma_rt , Form("#sigma_{RT1}^{%i}",q2Bin) , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));
//     c_alpha_rt1.push_back( constrainVar(alpha_rt1, Form("#alpha_{RT1}^{%i}",q2Bin) , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));
//     c_alpha_rt2.push_back( constrainVar(alpha_rt2, Form("#alpha_{RT2}^{%i}",q2Bin) , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));
//     c_n_rt1.push_back(     constrainVar(n_rt1    , Form("n_{RT1}^{%i}",q2Bin)      , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));
//     c_n_rt2.push_back(     constrainVar(n_rt2    , Form("n_{RT2}^{%i}",q2Bin)      , wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));

    RooAbsPdf* dcb_rt;
    RooRealVar* sigma_rt2 = new RooRealVar (Form("#sigma_{RT2}^{%i}",years[iy] ), "sigmaRT2"  ,   0 , 0,   0.12, "GeV");
    RooRealVar* f1rt      = new RooRealVar (Form("f^{RT%i}",years[iy])          , "f1rt"      ,   0 , 0.,  1.);
    if (q2Bin >= 5){
      sigma_rt2-> setVal(wsp_mcmass[iy]->var(Form("#sigma_{RT2}^{%i}",q2Bin))->getVal() );
      f1rt     -> setVal(wsp_mcmass[iy]->var(Form("f^{RT%i}", q2Bin))->getVal() );
//       c_sigma_rt2.push_back(  constrainVar(sigma_rt2 , Form("#sigma_{RT2}^{%i}",q2Bin), wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));
//       c_f1rt     .push_back(  constrainVar(f1rt      , Form("f^{RT%i}"         ,q2Bin), wsp_mcmass[iy], years[iy], true, c_vars_rt, c_pdfs_rt));
      dcb_rt = createRTMassShape(mass, mean_rt, sigma_rt, sigma_rt2, alpha_rt1, alpha_rt2, n_rt1, n_rt2 ,f1rt, wsp_mcmass[iy], years[iy] );
    } 
    else    
        dcb_rt = createRTMassShape(mass, mean_rt, sigma_rt, alpha_rt1, alpha_rt2, n_rt1, n_rt2 , wsp_mcmass[iy], years[iy] );
   
    /// create WT component
    RooRealVar* mean_wt     = new RooRealVar (Form("mean_{WT}^{%i}",years[iy])      , "masswt"     ,  wsp_mcmass[iy]->var(Form("mean_{WT}^{%i}", q2Bin))->getVal()    ,      5,    6, "GeV");
    RooRealVar* sigma_wt    = new RooRealVar (Form("#sigma_{WT1}^{%i}",years[iy])   , "sigmawt"    ,  wsp_mcmass[iy]->var(Form("#sigma_{WT1}^{%i}", q2Bin))->getVal() ,      0,    1, "GeV");
    RooRealVar* alpha_wt1   = new RooRealVar (Form("#alpha_{WT1}^{%i}",years[iy] )  , "alphawt1"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT1}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* alpha_wt2   = new RooRealVar (Form("#alpha_{WT2}^{%i}",years[iy] )  , "alphawt2"   ,  wsp_mcmass[iy]->var(Form("#alpha_{WT2}^{%i}", q2Bin))->getVal() ,      0,   10 );
    RooRealVar* n_wt1       = new RooRealVar (Form("n_{WT1}^{%i}",years[iy])        , "nwt1"       ,  wsp_mcmass[iy]->var(Form("n_{WT1}^{%i}", q2Bin))->getVal()      ,      0., 100.);
    RooRealVar* n_wt2       = new RooRealVar (Form("n_{WT2}^{%i}",years[iy])        , "nwt2"       ,  wsp_mcmass[iy]->var(Form("n_{WT2}^{%i}", q2Bin))->getVal()      ,      0., 100.);
    /// create constrain WT 
    wsp_mcmass[iy]->loadSnapshot(Form("reference_fit_WT_%i",q2Bin));
//     c_sigma_wt.push_back(  constrainVar(sigma_wt , Form("#sigma_{WT1}^{%i}",q2Bin) , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt));
//     c_alpha_wt1.push_back( constrainVar(alpha_wt1, Form("#alpha_{WT1}^{%i}",q2Bin) , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt));
//     c_alpha_wt2.push_back( constrainVar(alpha_wt2, Form("#alpha_{WT2}^{%i}",q2Bin) , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt));
//     c_n_wt1.push_back(     constrainVar(n_wt1    , Form("n_{WT1}^{%i}",q2Bin)      , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt));
//     c_n_wt2.push_back(     constrainVar(n_wt2    , Form("n_{WT2}^{%i}",q2Bin)      , wsp_mcmass[iy], years[iy], true, c_vars_wt, c_pdfs_wt));
// 
    RooAbsPdf* dcb_wt = createWTMassShape(mass, mean_wt, sigma_wt, alpha_wt1, alpha_wt2, n_wt1, n_wt2 , wsp_mcmass[iy], years[iy] );
      

//     PdfCBShape* cbs = new PdfCBShape("cbs", "Crystal Ball shape", 
//                                       *mass, 
//                                       *mean_rt, *sigma_rt, *alpha_rt1, *n_rt1, *alpha_rt2, *n_rt2);
//     std::cout << "[yVal of CBS]: \t" << cbs->eval() << std::endl;
    std::cout << "[yVal of CBS]: \t" << dcb_rt->getVal() << std::endl;
      
    PdfSigMass* PDF_sig_ang_mass_unc;
    if (q2Bin < 5)  
        PDF_sig_ang_mass_unc = new PdfSigMass( ("PDF_sig_ang_mass_unc_"+shortString+"_"+year).c_str(),
                                               ("PDF_sig_ang_mass_unc_"+year).c_str(),
                                               *mass,
                                               *mean_rt, *sigma_rt, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2,
                                               *mean_wt, *sigma_wt, *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
      		                               *mFrac,
      		                               *dcb_rt,
      		                               *dcb_wt  
      		                              );
    else  
        PDF_sig_ang_mass_unc = new PdfSigMass( ("PDF_sig_ang_mass_unc_"+shortString+"_"+year).c_str(),
                                               ("PDF_sig_ang_mass_unc_"+year).c_str(),
                                               *mass,
                                               *mean_rt, *sigma_rt, *sigma_rt2, *alpha_rt1, *alpha_rt2, *n_rt1, *n_rt2, *f1rt,
                                               *mean_wt, *sigma_wt,             *alpha_wt1, *alpha_wt2, *n_wt1, *n_wt2,
      		                               *mFrac,
      		                               *dcb_rt ,
      		                               *dcb_wt
      		                               );

    PDF_sig_ang_mass.push_back(PDF_sig_ang_mass_unc); 

    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    if ( !data[iy][0] || data[iy][0]->IsZombie() ) {
      cout<<"Dataset " << firstSample  << " not found in file: "<<filename_data<<endl;
      return;
     }
    map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",firstSample)).c_str(), data[iy][0]) );
//     simPdf        ->addPdf(*dcb_rt, ("data"+year+Form("_subs%d",firstSample)).c_str());
    std::cout << "[yVal of simFit]: \t" << PDF_sig_ang_mass[iy]->getVal() << std::endl;
    simPdf        ->addPdf(*PDF_sig_ang_mass[iy], ("data"+year+Form("_subs%d",firstSample)).c_str());
  }
   
  TFile* fout = new TFile(("simFitResults4d/simFitResult_recoMC_fullMass" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"RECREATE");
  
  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data", 
                            reco_vars,
                            Index(sample), 
                            Import(map)); 
  RooDataSet* combData = 0;
  RooAbsReal* nll = 0;

  for (uint is = firstSample; is <= lastSample; is++) {

    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }
    
    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
    if (nSample>0) cout<<"Fitting subsample "<<is+1<<" with "<<combData->numEntries()<<" entries"<<endl;
    else cout<<"Fitting full MC sample with "<<combData->numEntries()<<" entries"<<endl;

    TStopwatch subTime;
    cout << "create nll" << endl;
    nll = simPdf->createNLL(*combData,
                            RooFit::Extended(kFALSE),
                            RooFit::NumCPU(1)
//                             RooFit::Constrain(c_vars)
                            );
         
    RooMinimizer m(*nll) ;   
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    m.setPrintLevel(-1);
    m.setPrintEvalErrors(-1);
    m.setMinimizerType("Minuit2");

    subTime.Start(true);
    m.setStrategy(0);
    m.migrad() ;
    m.hesse() ;

    subTime.Stop();
    cout << "fitting done  " << subTime.CpuTime() << endl;
//     m.setStrategy(2);
//     m.migrad() ;
//     m.hesse() ;
    
    RooFitResult* fitResult = m.save(("result_" + shortString + Form("subs%d",is)).c_str()) ; 
    fitResult->Print("v");

    // Save fit results in file
    if (save) {
      cout << fout->GetName() << endl;
      fout->cd();
      fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
    }
    
  }  

  fout->Close();

  string plotString = shortString + "_" + all_years;
  if (nSample>0) plotString = plotString + Form("_s%i",nSample);

  int confIndex = 2*nBins*parity  + q2Bin;
  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections 
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),3000,1400);
  c[confIndex]->Divide(4, years.size());
  
  cout<<"plotting 4d canvas"<<endl;
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    std::vector<RooPlot*> frames;
    frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);

    cout<<"canvas ready"<<endl;
    for (unsigned int fr = 0; fr < frames.size(); fr++){
        cout<<"fr " << fr<<endl;
        combData->plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(40), Cut(("sample==sample::data"+year+"_subs0").c_str()), Name(("plData"+year).c_str()));
         
//         mass->setBins(20) ;
//         RooDataHist* projData;
//         if (fr==0) projData = new RooDataHist ("projData","projData",RooArgSet(*ctK, *ctL, *phi),*combData) ;
        simPdf  ->plotOn(frames[fr], Slice(sample, ("data"+year+"_subs0").c_str()), 
                                     ProjWData(RooArgSet(sample), *combData), 
                                     LineWidth(1), 
                                     Name(("plPDF"+year).c_str()), 
                                     NumCPU(4));
        if (fr == 0) { 
          leg->AddEntry(frames[fr]->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
          leg->AddEntry(frames[fr]->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
        }
        c[confIndex]->cd(iy*4+fr+1);
        gPad->SetLeftMargin(0.19); 
        frames[fr]->Draw();
        leg->Draw("same");
        break;
    }
  }
  c[confIndex]->SaveAs( ("plotSimFit4d_d/simFitResult_recoMC_fullMass_" + plotString +  ".pdf").c_str() );

}



void simfit_recoMC_fullMassBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, std::vector<int> years, std::map<int,float> scale_to_data, double fac1, double fac4, double fac5, double base1, double base4, double base5, double max1, double max4, double max5)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);
  else
    simfit_recoMC_fullMassBin(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);
  if ( argc > 2 ) parity  = atoi(argv[2]);

  double fac1 = 1;
  double fac4 = 1;
  double fac5 = 1;
  double base1 = 3;
  double base4 = 3;
  double base5 = 3;
  double max1 = 0;
  double max4 = 0;
  double max5 = 200;

  if ( argc > 3  ) fac1  = atof(argv[3]) / 1000.0;
  if ( argc > 4  ) fac4  = atof(argv[4]) / 1000.0;
  if ( argc > 5  ) fac5  = atof(argv[5]) / 1000.0;
  if ( argc > 6  ) base1 = atof(argv[6]) / 1000.0;
  if ( argc > 7  ) base4 = atof(argv[7]) / 1000.0;
  if ( argc > 8  ) base5 = atof(argv[8]) / 1000.0;
  if ( argc > 9  ) max1  = atof(argv[9]);
  if ( argc > 10 ) max4  = atof(argv[10]);
  if ( argc > 11 ) max5  = atof(argv[11]);

  bool multiSample = false;
  uint nSample = 0;
  if ( argc > 12 && atoi(argv[12]) > 0 ) multiSample = true;
  if ( argc > 13 ) nSample = atoi(argv[13]);

  if (nSample==0) multiSample = false;

  bool plot = true;
  bool save = true;

  if ( argc > 14 && atoi(argv[14]) == 0 ) plot = false;
  if ( argc > 15 && atoi(argv[15]) == 0 ) save = false;

  std::vector<int> years;
  if ( argc > 16 && atoi(argv[16]) != 0 ) years.push_back(atoi(argv[16]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 17 && atoi(argv[17]) != 0 ) years.push_back(atoi(argv[17]));
  if ( argc > 18 && atoi(argv[18]) != 0 ) years.push_back(atoi(argv[18]));

  cout <<  "q2Bin       " << q2Bin        << endl;
  cout <<  "parity      " << parity       << endl;
  cout <<  "multiSample " << multiSample  << endl;
  cout <<  "nSample     " << nSample      << endl;
  cout <<  "plot        " << plot         << endl;
  cout <<  "save        " << save         << endl;
  cout <<  "years[0]    " << years[0]     << endl;
//   cout << years[1] << endl;
//   cout << years[2] << endl;


  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  std::map<int,float> scale_to_data;
  // https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
  scale_to_data.insert(std::make_pair(2016, 0.01 )); // *2 since we are using only odd/even events, second factor is "data-driven"
  scale_to_data.insert(std::make_pair(2017, 0.01 ));
  scale_to_data.insert(std::make_pair(2018, 0.01 ));
//   scale_to_data.insert(std::make_pair(2016, 0.006*2 /2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
//   scale_to_data.insert(std::make_pair(2017, 0.005*2 /2.05 ));
//   scale_to_data.insert(std::make_pair(2018, 0.007*2 /1.9  ));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);
  else
    simfit_recoMC_fullMassBin1(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);

  return 0;

}
