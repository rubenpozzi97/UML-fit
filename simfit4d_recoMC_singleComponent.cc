#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <TSystem.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooCmdArg.h>
#include <RooConstVar.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>

#include "RooDoubleCBFast.h"
#include "PdfRT.h"
#include "PdfWT.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* c [4*nBins];

RooGaussian* constrainVar(RooRealVar* var, string inVarName,  RooWorkspace *w, int year){
    RooGaussian* gauss_constr = new RooGaussian(  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  *var,  
                                                  RooConst( w->var(inVarName.c_str())->getVal()  ), 
                                                  RooConst( w->var(inVarName.c_str())->getError())
                                                 ); 
    gauss_constr ->Print();
    return gauss_constr;                 
}

RooPlot* prepareFrame(RooPlot* frame){
    frame->GetYaxis()->SetTitleOffset(1.8);
    frame->SetMaximum(frame->GetMaximum()*1.15);
    frame->SetMinimum(0);
    return frame;
}


void simfit4d_recoMC_singleComponentBin(int q2Bin, int parity, int tagFlag, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data)
{
  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency
  string shortString = Form("b%ip%it%i",q2Bin,parity,tagFlag);
  cout<<"Conf: "<<shortString<<endl;
  string datasetString = (tagFlag==1?"data_ctRECO":"data_wtRECO");
  datasetString = datasetString + Form((parity==1?"_ev_b%i":"_od_b%i"),q2Bin);
  string effString = Form((tagFlag==0?"effWHist_b%ip%i":"effCHist_b%ip%i"),q2Bin,parity);
  string intHistString = "MCint_"+shortString;
  string all_years = "";
  string year = ""; 
  string stat = datalike ? "_dataStat":"_MCStat";

  TFile* fin_data, *fin_eff, *fin_MC_mass;
  RooWorkspace* wsp;
  std::vector<RooDataSet*> data;
  std::vector<RooAbsReal*> eff;
  std::vector<TH3D*> effHist;
  std::vector<TH1D*> intHist;
  std::vector< std::vector<double> > intVec(years.size(), std::vector<double>(0));
  std::vector<RooAbsPdf*> PDF_sig_ang_singleComponent (0), PDF_sig_mass_singleComponent (0);
  std::vector<RooGaussian*> c_sigma_rt1, c_sigma_rt2, c_mean_rt, c_f1rt;
  std::vector<RooGaussian*> c_sigma_wt, c_mean_wt, c_alpha_wt1, c_alpha_wt2, c_n_1, c_n_2;
  RooArgSet c_vars;
  

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

  RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);
  RooRealVar* mass = new RooRealVar("tagged_mass","mass", 5.,5.6);
  RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
  RooArgSet reco_vars (*ctK, *ctL, *phi, *mass, *rand);

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl  = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1  = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2  = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3  = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    sample.defineType(("data"+year).c_str());
    all_years += year;
  }
  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));

    // import data (or MC as data proxy)
    string filename = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/recoMCDataset_b%i_%i_tagged.root", years[iy], q2Bin, years[iy]); 
    fin_data = TFile::Open( filename.c_str() ) ;
    if ( !fin_data || !fin_data->IsOpen() ) {
      cout << "File not found: " << filename << endl;
      return;
    }
    wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i", q2Bin, 1-parity ) );
    if ( !wsp || wsp->IsZombie() ) {
      cout<<"Workspace not found in file: "<<filename<<endl; 
      return;
    }
  
    // create roodataset (in case data-like option is selected, only import the correct % of data)
    if (datalike)
      data.push_back( (RooDataSet*)wsp->data(datasetString.c_str())->reduce( RooArgSet(reco_vars), Form("rand < %f", scale_to_data[years[iy]] ))) ;
    else
      data.push_back( (RooDataSet*)wsp->data(datasetString.c_str())) ;
    if ( !data[iy] || data[iy]->IsZombie() ) {
      cout<<"Dataset "<<datasetString<<" not found in file: "<<filename<<endl;
      return;
    }

    // import KDE efficiency histograms and partial integral histograms
    filename = "../eff-KDE/files_forSimFit/KDEeff_b";
    filename = filename + Form((parity==0 ? "%i_ev_%i.root" : "%i_od_%i.root"),q2Bin,years[iy]);
    fin_eff  = new TFile( filename.c_str(), "READ");
    if ( !fin_eff || !fin_eff->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      return;
    }
    effHist.push_back( (TH3D*)fin_eff->Get(effString.c_str()));
    if ( !effHist[iy] || effHist[iy]->IsZombie() ) {
      cout<<"Efficiency histogram "<< effString <<" not found in file: "<< filename <<endl;
      return;
    }

    // create efficiency functions
    RooDataHist* effData = new RooDataHist(("effData_"+shortString+"_"+year).c_str(),"effData",vars,effHist[iy]);
    eff.push_back( new RooHistFunc(("eff_"+shortString+"_"+year).c_str(),
                                   ("eff"+year).c_str() ,
                                   vars,
                                   *effData,
                                   1));

    // import precomputed integrals and fill a std::vector
    intHist.push_back( (TH1D*)fin_eff->Get(intHistString.c_str()));
    intVec.push_back (vector<double> (0));
    if ( !intHist[iy] || intHist[iy]->IsZombie() ) {
      cout << "Integral histogram " << intHistString <<" not found in file: "<< filename << endl << "Abort" << endl;
      return;
    } else if ( strcmp( intHist[iy]->GetTitle(), effHist[iy]->GetTitle() ) ) {
    // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
      cout << "Integral histogram "<< intHistString << " is incoherent with efficiency " << effString << " in file: " << filename << endl;
      cout << "Efficiency conf: "  << effHist[iy]->GetTitle() << endl;
      cout << "Integral conf: "    << intHist[iy]->GetTitle() << endl << "Abort" << endl;
      return;
    } 
    else {
      for (int i=1; i<=intHist[iy]->GetNbinsX(); ++i) {
        intVec[iy].push_back(intHist[iy]->GetBinContent(i));
      }
    }

    // define angular PDF for signal and only one tag component, using the custom class
    // efficiency function and integral values are passed as arguments
    if (tagFlag==1) PDF_sig_ang_singleComponent.push_back( new PdfRT(("PDF_sig_ang_singleComponent_"+shortString+"_"+year).c_str(),
                                                                     ("PDF_sig_ang_singleComponent_"+year).c_str(),
      							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p, *eff[iy], intVec[iy]));
    else            PDF_sig_ang_singleComponent.push_back( new PdfWT(("PDF_sig_ang_singleComponent_"+shortString+"_"+year).c_str(),
                                                                     ("PDF_sig_ang_singleComponent_"+year).c_str(),
  							  *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p, *eff[iy], intVec[iy]));
    

    // import mass parametrisation from fits to the MC
    filename = Form("results_fits_%i_newCB.root",years[iy]);
    fin_MC_mass = new TFile(filename.c_str(), "READ");
    if ( !fin_MC_mass || !fin_MC_mass->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      return;
    }
//     string tagstr = "RT";
//     if (tagFlag==0) tagstr = "WT";
//     string fitresname = Form("results_%s_%i",tagstr.c_str(),q2Bin);
    
    RooWorkspace *w = (RooWorkspace*)fin_MC_mass->Get("w");
    if ( !w || w->IsZombie() ) {
      cout << "Workspace not found in fitResult file: " << filename << endl;
      fin_MC_mass->ls();
      return;
    }
    
    // build mass PDF for the RT component
    if (tagFlag==1){
      RooRealVar* meanrt   = new RooRealVar (Form("mean^{%i}",years[iy])       , "massRT"    ,   w->var(Form("mean^{RT%i}",q2Bin))->getVal() ,      5,    6, "GeV");
      RooRealVar* sigmart1 = new RooRealVar (Form("#sigma_{1}^{%i}",years[iy] ), "sigmaRT1"  ,   w->var(Form("#sigma_{1}^{RT%i}",q2Bin))->getVal() ,      0,    1, "GeV");
      RooRealVar* sigmart2 = new RooRealVar (Form("#sigma_{2}^{%i}",years[iy] ), "sigmaRT2"  ,   w->var(Form("#sigma_{2}^{RT%i}",q2Bin))->getVal() ,      0,   0.12, "GeV");
      RooRealVar* f1rt     = new RooRealVar (Form("f^{%i}",years[iy])          , "f1rt"      ,   w->var(Form("f^{RT%i}"   ,q2Bin))->getVal()       ,      0.,    1.);
  
      RooGaussian* sigGauss1 = new RooGaussian (Form("dg_gaus1_%i", years[iy]) , "gaus1"      , *mass, *meanrt, *sigmart1);
      RooGaussian* sigGauss2 = new RooGaussian (Form("dg_gaus2_%i", years[iy]) , "gaus2"      , *mass, *meanrt, *sigmart2);
      RooAddPdf *theRTgauss  = new RooAddPdf (Form("doublegaus_%i", years[iy]) , "gaus1+gaus2", RooArgList(*sigGauss1,*sigGauss2), RooArgList(*f1rt));
      
      /// create constraints for the RT component
      c_sigma_rt1.push_back( constrainVar(sigmart1, Form("#sigma_{1}^{RT%i}",q2Bin), w, years[iy]));
      c_sigma_rt2.push_back( constrainVar(sigmart2, Form("#sigma_{2}^{RT%i}",q2Bin), w,  years[iy]));
      c_mean_rt.push_back(   constrainVar(meanrt,   Form("mean^{RT%i}",q2Bin),       w,  years[iy]));
      c_f1rt.push_back(      constrainVar(f1rt,     Form("f^{RT%i}"   ,q2Bin),       w,  years[iy]));
      c_vars.add(*c_sigma_rt1[iy]);
      c_vars.add(*c_sigma_rt2[iy]);
      c_vars.add(*c_mean_rt[iy]);
      c_vars.add(*c_f1rt[iy]);
      
      PDF_sig_mass_singleComponent.push_back( new RooProdPdf  (("mass_pdf"+year).c_str() , 
                                                               ("mass_pdf"+year).c_str() , 
                                                               RooArgList(*theRTgauss, 
                                                                          *c_sigma_rt1[iy],
                                                                          *c_sigma_rt2[iy], 
                                                                          *c_mean_rt[iy], 
                                                                          *c_f1rt[iy]  
                                                               )));  
    }  
    // build mass PDF for the WT component
    else{

      RooRealVar* meanwt   = new RooRealVar (Form("mean^{%i}",years[iy])       , "massWT"    ,   w->var(Form("mean^{WT%i}", q2Bin))->getVal()       ,      5,    6, "GeV");
      RooRealVar* sigmaCB  = new RooRealVar (Form("#sigma_{CB}^{%i}",years[iy]), "sigmaCB"   ,   w->var(Form("#sigma_{CB}^{WT%i}", q2Bin))->getVal(),      0,    1, "GeV");
      RooRealVar* alpha1   = new RooRealVar (Form("#alpha_{1}^{%i}",years[iy] ), "alpha1"    ,   w->var(Form("#alpha_{1}^{WT%i}", q2Bin))->getVal() ,      0,   10 );
      RooRealVar* alpha2   = new RooRealVar (Form("#alpha_{2}^{%i}",years[iy] ), "alpha2"    ,   w->var(Form("#alpha_{2}^{WT%i}", q2Bin))->getVal() ,      0,   10 );
      RooRealVar* n1       = new RooRealVar (Form("n_{1}^{%i}",years[iy])      , "n1"        ,   w->var(Form("n_{1}^{WT%i}", q2Bin))->getVal()      ,      0.,  20.);
      RooRealVar* n2       = new RooRealVar (Form("n_{2}^{%i}",years[iy])      , "n2"        ,   w->var(Form("n_{2}^{WT%i}", q2Bin))->getVal()      ,      0.,  20.);
      RooDoubleCBFast* doublecb = new RooDoubleCBFast ( Form("dcb_%i", years[iy]), "dcb"     , *mass, *meanwt, *sigmaCB, *alpha1, *n1, *alpha2, *n2);
      
      /// create constraints for the WT component
      c_mean_wt.push_back(   constrainVar(meanwt, Form("mean^{WT%i}"       ,q2Bin), w, years[iy]));
      c_sigma_wt.push_back(  constrainVar(sigmaCB,Form("#sigma_{CB}^{WT%i}",q2Bin), w, years[iy]));
      c_alpha_wt1.push_back( constrainVar(alpha1, Form("#alpha_{1}^{WT%i}" ,q2Bin), w, years[iy]));
      c_alpha_wt2.push_back( constrainVar(alpha2, Form("#alpha_{2}^{WT%i}" ,q2Bin), w, years[iy]));
      c_n_1.push_back(       constrainVar(n1,     Form("n_{1}^{WT%i}"      ,q2Bin), w, years[iy]));
      c_n_2.push_back(       constrainVar(n2,     Form("n_{2}^{WT%i}"      ,q2Bin), w, years[iy]));
      c_vars.add(*c_mean_wt[iy]);
      c_vars.add(*c_sigma_wt[iy]);
      c_vars.add(*c_alpha_wt1[iy]);
      c_vars.add(*c_alpha_wt2[iy]);
      c_vars.add(*c_n_1[iy]);
      c_vars.add(*c_n_2[iy]);
      
      PDF_sig_mass_singleComponent.push_back( new RooProdPdf  (("mass_pdf"+year).c_str() , 
                                                               ("mass_pdf"+year).c_str() , 
                                                               RooArgList( *doublecb, 
                                                                           *c_mean_wt[iy], 
                                                                           *c_sigma_wt[iy],
                                                                           *c_alpha_wt1[iy], 
                                                                           *c_alpha_wt2[iy], 
                                                                           *c_n_1[iy], 
                                                                           *c_n_2[iy]  
                                                               )));  
    }                                                           
    // insert sample in the category map, to be imported in the combined dataset
    map.insert(map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year).c_str(), data[iy]));
    // associate model with the data
    RooProdPdf* mass_ang_pdf = new RooProdPdf(("mass_ang_pdf_"+year).c_str(), ("mass_ang_pdf_"+year).c_str(), RooArgList(*PDF_sig_ang_singleComponent[iy], *PDF_sig_mass_singleComponent[iy]));
    simPdf->addPdf(*mass_ang_pdf, ("data"+year).c_str());
    
    fin_data    -> Close();
    fin_eff     -> Close();
    fin_MC_mass -> Close();
  }

  for(auto it = map.cbegin(); it != map.cend(); ++it)
    std::cout << "dataset: " << it->first << ", with n entries: " << it->second->sumEntries() << "\n";

  // to start the fit, parameters are restored to the center of the parameter space
  Fl ->setVal(0.5);
  P1 ->setVal(0);
  P2 ->setVal(0);
  P3 ->setVal(0);
  P4p->setVal(0);
  P5p->setVal(0);
  P6p->setVal(0);
  P8p->setVal(0);

  // Construct combined dataset in (x,sample)
  RooDataSet combData ("combData", "combined data", 
                         reco_vars,
                         Index(sample),
                         Import(map)); 


  // perform fit in two steps:
  // first with strategy=0 and no MINOS, to get close to the likilihood maximum
  simPdf->fitTo( combData,
                 Minimizer("Minuit2","migrad"), 
                 Extended(false), 
//                  Timer(true),
                 NumCPU(1),
                 Hesse(true),
                 Strategy(0),
//                  Offset(true),
                 Constrain(c_vars)
                 );  
  // second with full accuracy (strategy=2) and MINOS, to get the best result
  RooFitResult* fitResult = simPdf->fitTo( combData,
                                           Minimizer("Minuit2","migrad"),
//                                            Extended(false), 
                                           Save(true),
                                           NumCPU(1),
                                           Hesse(true),
                                           Strategy(2),
                                           Minos(true),
                                           Offset(true),
                                           Constrain(c_vars)
                                         );
                                         
  // The two step fit seemed necessary to get convergence in all q2 bins
  // if one observes that the convergence is obtained just running the second step, it is better to avoid the first step (to gain computing time)
  // Offset(true) parameter make the FCN value to be defined at 0 at the begin of the minimization
  // it is needed to get convergence in all q2 bins, and avoid the "machine precision reached" errors
  // which are due to too high FCN absolute values compared to the minimisation steps

  fitResult->Print("v");

  if (save) {
    // Save fit results in file
    TFile* fout = new TFile(("simFit4dResults/simFit4dResult_recoMC_singleComponent" + all_years + stat + ".root").c_str(),"UPDATE");
    fitResult->Write(("simFitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close();
  }

  if (!plot) return;

  int confIndex = 2*nBins*parity + nBins*tagFlag + q2Bin;
  string longString  = (tagFlag==1?"Fit to correctly tagged events":"Fit to wrongly tagged events");
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,1400);
  c[confIndex]->Divide(4, years.size());
  
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    std::vector<RooPlot*> frames;
    frames.push_back( prepareFrame( ctK ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( ctL ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( phi ->frame(Title((longString+year).c_str())) ));
    frames.push_back( prepareFrame( mass ->frame(Title((longString+year).c_str()))));
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
    leg->SetTextSize(0.03);

    for (unsigned int fr = 0; fr < frames.size(); fr++){
        combData.plotOn(frames[fr], MarkerColor(kRed+1), LineColor(kRed+1), Binning(40), Cut(("sample==sample::data"+year).c_str()), Name(("plData"+year).c_str()));
        simPdf ->plotOn(frames[fr], Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1), Name(("plPDF"+year).c_str()));
        if (fr == 0) {
          leg->AddEntry(frames[fr]->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
          leg->AddEntry(frames[fr]->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
        }
        c[confIndex]->cd(iy*4+fr+1);
        gPad->SetLeftMargin(0.19); 
        frames[fr]->Draw();
        leg->Draw("same");
    }
  }

  c[confIndex]->SaveAs( ("plotSimFit4d_d/simFit4dResult_recoMC_singleComponent_" + shortString + "_" + all_years + stat + ".pdf").c_str() );
}


void simfit4d_recoMC_singleComponentBin2(int q2Bin, int parity, int tagFlag, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data)
{
  if ( tagFlag==-1 )
    for (tagFlag=0; tagFlag<2; ++tagFlag)
      simfit4d_recoMC_singleComponentBin(q2Bin, parity, tagFlag, plot, save, datalike,  years, scale_to_data);
  else
    simfit4d_recoMC_singleComponentBin(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);
}

void simfit4d_recoMC_singleComponentBin1(int q2Bin, int parity, int tagFlag, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit4d_recoMC_singleComponentBin2(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);
  else
    simfit4d_recoMC_singleComponentBin2(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively
  // tag format: [0] mistagged
  //             [1] correctly-tagged
  //             [-1] for each tag recursively

//   gSystem->Load("doubleCB/libRooDoubleCBFast.so");

  int q2Bin   = -1;
  int parity  = -1; 
  int tagFlag = -1;

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);
  if ( argc >= 4 ) tagFlag = atoi(argv[3]);

  bool plot = true;
  bool save = true;
  bool datalike = true;

  if ( argc >= 5 && atoi(argv[4]) == 0 ) plot     = false;
  if ( argc >= 6 && atoi(argv[5]) == 0 ) save     = false;
  if ( argc >= 7 && atoi(argv[6]) == 0 ) datalike = false;

  std::vector<int> years;
  if ( argc >= 8 && atoi(argv[7]) != 0 ) years.push_back(atoi(argv[7]));
  else {
    cout << " no specific years selected, using default: 2016";
    years.push_back(2016);
  }
  if ( argc >= 9  && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));
  if ( argc >= 10 && atoi(argv[9]) != 0 ) years.push_back(atoi(argv[9]));

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;
  if ( tagFlag < -1 || tagFlag > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;
  if ( tagFlag==-1 ) cout << "Running both the tag conditions"  << endl;
  if ( datalike )    cout << "Considering data-like statistics" << endl;

  std::map<int,float> scale_to_data;
  scale_to_data.insert(std::make_pair(2016, 0.01*2)); // *2 since we are using only odd/even events 
  scale_to_data.insert(std::make_pair(2017, 0.01*2));
  scale_to_data.insert(std::make_pair(2018, 0.015*2));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit4d_recoMC_singleComponentBin1(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);
  else
    simfit4d_recoMC_singleComponentBin1(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);

  return 0;

}
