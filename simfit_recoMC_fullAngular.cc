#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
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

#include "ShapeSigAng.h"

#include "PdfSigAng.h"
#include "BoundCheck.h"
#include "BoundDist.h"
#include "Penalty.h"

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

void simfit_recoMC_fullAngularBin(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, std::vector<int> years, std::map<int,float> scale_to_data, double fac1, double fac4, double fac5, double base1, double base4, double base5, double max1, double max4, double max5)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency

  string effCString = Form("effCHist_b%ip%i",q2Bin,parity);
  string effWString = Form("effWHist_b%ip%i",q2Bin,parity);
  string intCHistString = "MCint_"+shortString + "t1";
  string intWHistString = "MCint_"+shortString + "t0";
  string all_years = "";
  string year = ""; 
  string isample = ""; 
  string stat = nSample > 0 ? "_dataStat":"_MCStat";
  uint firstSample = ( multiSample || nSample==0 ) ? 0 : nSample-1;
  uint lastSample = nSample > 0 ? nSample-1 : 0;
  
  std::vector<TFile*> fin_data, fin_eff;
  std::vector<RooWorkspace*> wsp;
  std::vector<std::vector<RooDataSet*>> data;
  std::vector<RooAbsReal*> effC, effW;
  std::vector<TH3D*> effCHist, effWHist;
  std::vector<TH1D*> intCHist, intWHist;
  std::vector< std::vector<double> > intCVec(years.size(), std::vector<double>(0));
  std::vector< std::vector<double> > intWVec(years.size(), std::vector<double>(0));
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular (0);
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular_penalty (0);

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

  RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);
  RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
  RooRealVar* mass = new RooRealVar("mass","mass", 5.,5.6);
  RooArgSet reco_vars (*ctK, *ctL, *phi, *rand, *mass);
//   RooArgSet reco_vars (*ctK, *ctL, *phi, *rand);

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl    = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1    = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2    = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3    = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* P4p   = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p   = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p   = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p   = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* mFrac = new RooRealVar("mFrac","mistag fraction",1, 0, 2);
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
  RooSimultaneous* simPdf_penalty = new RooSimultaneous("simPdf_penalty", "simultaneous pdf with penalty term", sample);

  // Define boundary check (returning 0 in physical region and 1 outside)
  BoundCheck* boundary = new BoundCheck("bound","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // Define boundary distance calculator
  BoundDist* bound_dist = new BoundDist("bound","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,true,0,false);

  // Define penalty term (parameters set to zero and will be set sample-by-sample)
  Penalty* penTerm = new Penalty("penTerm","Penalty term",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,0,0,0,0);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
//     string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/newphi/effDataset_b%i_%i.root", years[iy], q2Bin, years[iy]); 
    string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/recoMCDataset_b%i_%i.root", years[iy], q2Bin, years[iy]); 
//     string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/newphi/recoMCDataset_b%i_%i_tagged_newphi.root", years[iy], q2Bin, years[iy]); 

    // import data (or MC as data proxy)
    fin_data.push_back( TFile::Open( filename_data.c_str() ) );
    if ( !fin_data[iy] || !fin_data[iy]->IsOpen() ) {
      cout << "File not found: " << filename_data << endl;
      return;
    }
    wsp.push_back( (RooWorkspace*)fin_data[iy]->Get(Form("ws_b%ip%i", q2Bin, 1-parity ) ) );
    if ( !wsp[iy] || wsp[iy]->IsZombie() ) {
      cout<<"Workspace not found in file: "<<filename_data<<endl;
      return;
    }
  

    // import KDE efficiency histograms and partial integr histograms
    string filename = "/eos/cms/store/user/fiorendi/p5prime/effKDE/";
    filename = filename + Form((parity==0 ? "%i/lmnr/KDEeff_b%i_ev_%i.root" : "%i/lmnr/KDEeff_b%i_od_%i.root"),years[iy],q2Bin,years[iy]);
//     string filename = "/eos/cms/store/user/fiorendi/p5prime/effKDE/";
//     filename = filename + Form((parity==0 ? "%i/lmnr/newphi/KDEeff_b%i_ev_%i.root" : "%i/lmnr/newphi/KDEeff_b%i_od_%i.root"),years[iy],q2Bin,years[iy]);
    fin_eff.push_back( new TFile( filename.c_str(), "READ" ));
    if ( !fin_eff[iy] || !fin_eff[iy]->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      return;
    }

    effCHist.push_back( (TH3D*)fin_eff[iy]->Get(effCString.c_str()));
    effWHist.push_back( (TH3D*)fin_eff[iy]->Get(effWString.c_str()));
    if ( !effCHist[iy] || effCHist[iy]->IsZombie() || !effWHist[iy] || effWHist[iy]->IsZombie() ) {
      cout<<"Efficiency histogram "<< effCString <<" or " << effWString << " not found in file: "<< filename <<endl;
      return;
    }

    // create efficiency functions
    RooDataHist* effCData = new RooDataHist(("effCData_"+shortString+"_"+year).c_str(),"effCData",vars,effCHist[iy]);
    RooDataHist* effWData = new RooDataHist(("effWData_"+shortString+"_"+year).c_str(),"effWData",vars,effWHist[iy]);
    effC.push_back( new RooHistFunc(("effC_"+shortString+"_"+year).c_str(),
                                    ("effC"+year).c_str() ,
                                    vars,
                                    *effCData,
                                    1));
    effW.push_back( new RooHistFunc(("effW_"+shortString+"_"+year).c_str(),
                                    ("effW"+year).c_str() ,
                                    vars,
                                    *effWData,
                                    1));

    // import precomputed integrals and fill a std::vector
    intCHist.push_back( (TH1D*)fin_eff[iy]->Get(intCHistString.c_str()));
    intWHist.push_back( (TH1D*)fin_eff[iy]->Get(intWHistString.c_str()));
    intCVec.push_back (vector<double> (0));
    intWVec.push_back (vector<double> (0));
    if ( !intCHist[iy] || intCHist[iy]->IsZombie() || !intWHist[iy] || intWHist[iy]->IsZombie() ) {
      cout << "Integral histogram " << intCHistString <<" or " << intWHistString << " not found in file: "<< filename << endl << "Abort" << endl;
      return;
    } else if ( strcmp( intCHist[iy]->GetTitle(), effCHist[iy]->GetTitle() ) || strcmp( intWHist[iy]->GetTitle(), effWHist[iy]->GetTitle() )) {
    // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
      cout << "Integral histograms are incoherent with efficiency in file: " << filename << endl;
      cout << "Efficiency (CT) conf: " << effCHist[iy]->GetTitle() <<endl;
      cout << "Integral (CT) conf: "   << intCHist[iy]->GetTitle() <<endl;
      cout << "Efficiency (WT) conf: " << effWHist[iy]->GetTitle() <<endl;
      cout << "Integral (WT) conf: "   << intWHist[iy]->GetTitle() <<endl;
      cout << "Abort"<<endl;
      return;
    } 
    else {
      for (int i=1; i<=intCHist[iy]->GetNbinsX(); ++i) {
        intCVec[iy].push_back(intCHist[iy]->GetBinContent(i));
      }
      for (int i=1; i<=intWHist[iy]->GetNbinsX(); ++i) {
        intWVec[iy].push_back(intWHist[iy]->GetBinContent(i));
      }
    }


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

    // define angular PDF for signal, using the custom class
    // efficiency function and integral values are passed as arguments

    RooAbsReal* ang_rt = new ShapeSigAng( ("PDF_sig_ang_rt_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_rt_"+year).c_str(),
         		                 *ctK,*ctL,*phi,
         		                 *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
         		                 *effC[iy], *effW[iy], intCVec[iy],intWVec[iy],
         		                 true
         		                 );
    
    RooAbsReal* ang_wt = new ShapeSigAng( ("PDF_sig_ang_wt_"+shortString+"_"+year).c_str(),
                                         ("PDF_sig_ang_wt_"+year).c_str(),
         		                 *ctK,*ctL,*phi,
         		                 *Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,
         		                 *effC[iy], *effW[iy], intCVec[iy],intWVec[iy],
         		                 false
         		                 );


    PDF_sig_ang_fullAngular.push_back( new PdfSigAng(("PDF_sig_ang_fullAngular_"+shortString+"_"+year).c_str(),
                                                     ("PDF_sig_ang_fullAngular_"+year).c_str(),
      		                                     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
         		                             *ang_rt, *ang_wt
      		                                     ) );
    // define PDF with penalty term
    PDF_sig_ang_fullAngular_penalty.push_back( new PdfSigAng(("PDF_sig_ang_fullAngular_penalty_"+shortString+"_"+year).c_str(),
							     ("PDF_sig_ang_fullAngular_penalty_"+year).c_str(),
							     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
         		                                     *ang_rt, *ang_wt,
							     *penTerm));

    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    if (multiSample) for (uint is = firstSample; is <= lastSample; is++) {
	if ( !data[iy][is] || data[iy][is]->IsZombie() ) {
	  cout<<"Dataset " << is  << " not found in file: "<<filename_data<<endl;
	  return;
	}
	map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
	simPdf->addPdf(*PDF_sig_ang_fullAngular[iy], ("data"+year+Form("_subs%d",is)).c_str());
	simPdf_penalty->addPdf(*PDF_sig_ang_fullAngular_penalty[iy], ("data"+year+Form("_subs%d",is)).c_str());
      }
    else {
      if ( !data[iy][0] || data[iy][0]->IsZombie() ) {
	cout<<"Dataset " << firstSample  << " not found in file: "<<filename_data<<endl;
	return;
      }
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",firstSample)).c_str(), data[iy][0]) );
      simPdf->addPdf(*PDF_sig_ang_fullAngular[iy], ("data"+year+Form("_subs%d",firstSample)).c_str());
      simPdf_penalty->addPdf(*PDF_sig_ang_fullAngular_penalty[iy], ("data"+year+Form("_subs%d",firstSample)).c_str());
    }

  }


  TFile* fout = new TFile(("simFitResults/newphi/simFitResult_recoMC_fullAngular" + all_years + stat + Form("_b%i_newShape.root", q2Bin)).c_str(),"UPDATE");

  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data", 
                            vars,
                            Index(sample), 
                            Import(map)); 
  RooDataSet* combData = 0;
  RooAbsReal* nll = 0;
  RooAbsReal* nll_penalty = 0;

  // Results' containers
  RooRealVar* fitTime = new RooRealVar("fitTime","fit time",0,"s");
  RooRealVar* co1 = new RooRealVar("co1","Coefficient 1",0);
  RooRealVar* co4 = new RooRealVar("co4","Coefficient 4",0);
  RooRealVar* co5 = new RooRealVar("co5","Coefficient 5",0);
  RooRealVar* boundDist = new RooRealVar("boundDist","Distance from boundary",0);
  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  RooArgSet savePars (*co1,*co4,*co5,*fitTime,*boundDist);
  savePars.add(pars);
  RooCategory resStatus ("resStatus","Status of the fit result");
  resStatus.defineType("convergent-positive-noPenalty",0);
  resStatus.defineType("convergent-positive",1);
  resStatus.defineType("convergent-negative",2);
  resStatus.defineType("notconvergent-positive",3);
  resStatus.defineType("notconvergent-negative",4);
  RooDataSet* subResults = 0;
  RooDataSet* subNoPen = new RooDataSet("subNoPen","subNoPen",savePars);
  RooDataSet* subPosConv = new RooDataSet("subPosConv","subPosConv",savePars);
  RooDataSet* subPosNotc = new RooDataSet("subPosNotc","subPosNotc",savePars);
  RooDataSet* subNegConv = new RooDataSet("subNegConv","subNegConv",savePars);
  RooDataSet* subNegNotc = new RooDataSet("subNegNotc","subNegNotc",savePars);

  // Timer for fitting time
  TStopwatch subTime;

  // counters to monitor results' status
  int cnt[9];
  for (int iCnt=0; iCnt<9; ++iCnt) cnt[iCnt] = 0;

  bool usedPenalty = false;
                         
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

//     for(auto it = map.cbegin(); it != map.cend(); ++it)
//       std::cout << "dataset: " << it->first << ", with n entries: " << it->second->sumEntries() << "\n";

    // to start the fit, parameters are restored to the center of the parameter space
    Fl ->setVal(0.5);
    P1 ->setVal(0);
    P2 ->setVal(0);
    P3 ->setVal(0);
    P4p->setVal(0);
    P5p->setVal(0);
    P6p->setVal(0);
    P8p->setVal(0);

    // set penalty term power parameter
    int combEntries = combData->numEntries();
    penTerm->setPower(power/combEntries);

    // define variables needed for adaptive procedure
    bool isPhysical = false;
    double coeff1 = 0;
    double coeff4 = 0;
    double coeff5 = 0;
    int totCoeff, iCoeff1, iCoeff4;
    double base1_corr = base1*sqrt(combEntries);
    double base4_corr = base4*sqrt(combEntries);
    double base5_corr = base5*sqrt(combEntries);
    if (base1_corr<min_base) base1_corr = min_base;
    if (base4_corr<min_base) base4_corr = min_base;
    if (base5_corr<min_base) base5_corr = min_base;

    bool inCTL4  = true;
    bool inCTL15 = true;

    nll = simPdf->createNLL(*combData,
                            RooFit::Extended(kFALSE),
                            RooFit::NumCPU(1)
                            );
         
    RooMinimizer m(*nll) ;
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    // m.setVerbose(kTRUE);
    m.setPrintLevel(-1);  
    m.setPrintEvalErrors(-1);
    //  Minuit2.setEps(1e-16) ;
    m.setMinimizerType("Minuit2");
     
    subTime.Start(true); 

    m.setStrategy(0);   
  //   m.setEvalErrorWall(false);
    m.migrad() ;
    m.hesse() ;
    // std::cout << std::endl;
    // std::cout << "######################### now strategy 2 #########################"<< std::endl;
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;   
    // m.minos() ;
    
    RooFitResult* fitResult = m.save(("result_" + shortString + Form("subs%d",is)).c_str()) ; 
    fitResult->Print("v");

    RooFitResult* fitResult_penalty = 0;
    usedPenalty = false;

 //    if ( fitResult->status()!=0 || fitResult->covQual()!=3 || boundary->getValV() > 0 ) {
//       usedPenalty = true;
// 
//       if ( !boundary->isInCTL4()  ) inCTL4  = false;
//       if ( !boundary->isInCTL15() ) inCTL15 = false;
// 
//       for (totCoeff=0; fac1*pow(base1_corr,totCoeff)<=maxCoeff; ++totCoeff) {
// 
// 	for (iCoeff1=totCoeff; iCoeff1>=0; --iCoeff1) {
// 	  coeff1 = fac1 * pow(base1_corr,iCoeff1);
// 	  if (max1>0 && coeff1>max1) continue;
// 	  // if ( inCTL15 ) {
// 	  //   if ( iCoeff1>0 ) continue;
// 	  //   coeff1 = 0;
// 	  // }
// 	  penTerm->setCoefficient(1,coeff1);
// 
// 	  // for (iCoeff4=totCoeff-iCoeff1; iCoeff4>=0; --iCoeff4) {
// 	  iCoeff4=totCoeff-iCoeff1; // new
// 	  {                         // new
// 
// 	    coeff4 = fac4 * pow(base4_corr,iCoeff4);
// 	    if (max4>0 && coeff4>max4) continue;
// 	    // if ( inCTL4 ) {
// 	    //   if ( iCoeff4>0 ) continue;
// 	    //   coeff4 = 0;
// 	    // }
// 
// 	    // coeff5 = fac5 * pow(base5_corr,totCoeff-iCoeff1-iCoeff4);
// 	    coeff5 = pow(coeff1,1.5) / 316.2; // new
// 
// 	    if (max5>0 && coeff5>max5) continue;
// 	    // if ( inCTL15 ) {
// 	    //   if ( totCoeff-iCoeff1-iCoeff4>0 ) continue;
// 	    //   coeff5 = 0;
// 	    // }
// 
// 	    penTerm->setCoefficient(4,coeff4);
// 	    penTerm->setCoefficient(5,coeff5);
// 
// 	    nll_penalty = simPdf_penalty->createNLL(*combData,
// 						    RooFit::Extended(kFALSE),
// 						    RooFit::NumCPU(1)
// 						    );
// 
// 	    RooMinimizer m_penalty (*nll_penalty) ;
// 	    m_penalty.optimizeConst(kTRUE);
// 	    m_penalty.setOffsetting(kTRUE);
// 	    // m_penalty.setVerbose(kTRUE);
// 	    m_penalty.setMinimizerType("Minuit2");
// 	    // m_penalty.setProfile(kTRUE);
// 	    m_penalty.setPrintLevel(-1);
// 	    m_penalty.setPrintEvalErrors(-1);
// 	    m_penalty.setStrategy(2);
//     
// 	    m_penalty.migrad() ;
// 	    m_penalty.hesse() ;
// 	    fitResult_penalty = m_penalty.save(Form("subRes_%s_%i",shortString.c_str(),is),Form("subRes_%s_%i",shortString.c_str(),is));
// 	    
// 	    // cout<<penTerm->getCoefficient(1)<<"\t"<<penTerm->getCoefficient(5)<<"\t"<<P5p->getValV()<<endl;
// 	    // fitResult_penalty->Print("v");
// 	    
// 	    if ( fitResult_penalty->status()==0 && fitResult_penalty->covQual()==3 ) {
// 	      if ( boundary->getValV()==0 ) {
// 		isPhysical = true;
// 		cout<<"P "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
// 		break;
// 	      } else cout<<"O "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
// 	    } else cout<<"N "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
// 	  }
// 	  if (isPhysical) break;
// 	}
// 	if (isPhysical) break;
//       }
//     
//     }
// 
//     subTime.Stop();
//     fitTime->setVal(subTime.CpuTime());
    // fitTime->setVal(subTime.RealTime());

//     co1->setVal(0);
//     co4->setVal(0);
//     co5->setVal(0);
// 
//     if (usedPenalty) {
//       if (isPhysical) {
// 	co1->setVal(coeff1);
// 	co4->setVal(coeff4);
// 	co5->setVal(coeff5);
//       }
//     }
// 
//     double boundCheck = boundary->getValV();
//     bool convCheck = false;
//     if (usedPenalty && fitResult_penalty->status()==0 && fitResult_penalty->covQual()==3) convCheck = true;
//     if (!usedPenalty && fitResult->status()==0 && fitResult->covQual()==3) convCheck = true;

//     TStopwatch distTime;
//     distTime.Start(true);
//     double boundDistVal = bound_dist->getValV();
//     distTime.Stop();
//     cout<<"Distance from boundary: "<<boundDistVal<<" (computed in "<<distTime.CpuTime()<<" s)"<<endl;
//     boundDist->setVal(boundDistVal);
// 
//     if (boundDistVal>0.02 && usedPenalty)
//       cout<<"WARNING high distance: "<<boundDistVal<<" with coeff1 "<<coeff1<<" coeff4 "<<coeff4<<" coeff5 "<<coeff5<<endl;
// 
//     ++cnt[8];
//     int iCnt = 0;
//     if (!convCheck) iCnt += 4;
//     if (boundCheck>0) iCnt += 2;
//     if (usedPenalty) iCnt += 1;
//     ++cnt[iCnt];
// 
//     if (boundCheck>0) {
//       if (convCheck) {
// 	subNegConv->add(savePars);
// 	cout<<"Converged in unphysical region"<<endl;
//       } else {
// 	subNegNotc->add(savePars);
// 	cout<<"Not converged (result in unphysical region)"<<endl;
//       }
//     } else {
//       if (convCheck) {
// 	if (usedPenalty) {
// 	  subPosConv->add(savePars);
// 	  cout<<"Converged with penalty term with coeff: "<<coeff1<<" "<<coeff4<<" "<<coeff5<<endl;
// 	} else {
// 	  subNoPen->add(savePars);
// 	  cout<<"Converged without penalty"<<endl;
// 	}
//       } else {
// 	subPosNotc->add(savePars);
// 	cout<<"Not converged (result in physical region)"<<endl;
//       }
//     }
// 
//     // Save fit results in file
//     if (save) {
//       fout->cd();
//       if (usedPenalty) fitResult_penalty->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
//       else fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
//     }
//   }  

//   if (multiSample) {
    subResults = new RooDataSet("subResults",
				"Results of RECO sub-sample fitting",
				savePars,Index(resStatus),
				Import("convergent-positive-noPenalty",*subNoPen),
				Import("convergent-positive",*subPosConv),
				Import("convergent-negative",*subNegConv),
				Import("notconvergent-positive",*subPosNotc),
				Import("notconvergent-negative",*subNegNotc));

    double time90quant = 0;
    double quant = 0;
    double totEntries = subResults->sumEntries();
    for (time90quant = 0; quant<0.9; time90quant += 0.1)
      quant = subResults->sumEntries(Form("fitTime<%.2f",time90quant))/totEntries;
    cout<<"Average fit time: "<<subResults->mean(*fitTime)<<" sec (90% quantile: "<<time90quant<<" sec)"<<endl;

    cout<<"Fitted subsamples: "<<cnt[8]<<" of which good: "<<cnt[0]+cnt[1]<<" ("<<cnt[1]<<" with the use of the penalty term)"<<endl;
    cout<<"Bad fits: "<<cnt[3]<<" converging outside physical region, "<<cnt[5]+cnt[7]<<" not converged ("<<cnt[5]<<" in ph region)"<<endl;
  }

//   if (save) {
//     RooWorkspace* wksp = new RooWorkspace(((multiSample?"wsMulti_":"ws_")+shortString+Form("_s%i_pow%.1f",nSample,power)).c_str(),
// 					  (multiSample?"Workspace with set of RECO subsample fit results":
// 					   (nSample>0?"Workspace with RECO subsample fit result":
// 					    "Workspace with full RECO fit result")));
// 
//     if (multiSample) {
//       wksp->import(*subResults);
//     } else {
//       wksp->import(*combData,Rename("data"));
//       wksp->import(*simPdf,RenameVariable(simPdf->GetName(),"pdf"),Silence());
//       if (usedPenalty) {
// 	wksp->import(*simPdf_penalty,RenameVariable(simPdf_penalty->GetName(),"pdfPen"),Silence(),RecycleConflictNodes());
// 	wksp->import(*penTerm,Silence(),RecycleConflictNodes());
//       }
//     }
// 
//     fout->cd();
//     wksp->Write();
//   }

//   fout->Close();

//   if (multiSample) {
//     TCanvas* cDist = new TCanvas (("cDist_"+shortString).c_str(),("cDist_"+shortString).c_str(),1800,1800);
//     RooPlot* fDist = boundDist->frame(Name("fDist"),Title("Distribution of results' distance fram boundary"),Range(0,0.1));
//     subNoPen->plotOn(fDist,Binning(50,0,0.1),LineColor(kBlue),MarkerColor(kBlue),MarkerStyle(19),DrawOption("XL"));
//     subPosConv->plotOn(fDist,Binning(50,0,0.1),LineColor(kRed),MarkerColor(kRed),MarkerStyle(19),DrawOption("XL"));
//     cDist->cd();
//     fDist->Draw();
// //     cDist->SaveAs( ("plotSimFit_d/newphi/recoBoundDist_" + shortString + "_" + all_years + Form("_f-%.3f-%.3f-%.3f_b-%.3f-%.3f-%.3f_m-%.0f-%.0f-%.0f.pdf",fac1,fac4,fac5,base1,base4,base5,max1,max4,max5)).c_str() );
//   }

//   if (!plot || multiSample) return;

  // For plotting the effective penalty term is used
//   Penalty* penTerm_eff = new Penalty(*penTerm,"penTerm_eff");
//   penTerm_eff->setPower(power);
//   RooFormulaVar* penLog = new RooFormulaVar("penLog","penLog","-1.0 * log(penTerm_eff)",RooArgList(*penTerm_eff));
// 
//   RooAbsReal* nll_penalty_plot = new RooAddition("nll_penalty_plot", "nll_penalty_plot", RooArgList(*nll,*penLog));
// 
//   double xZoom = 200.0;
//   if (nSample>0) xZoom = 2.0;

 //  cnll  = new TCanvas (("cnll_"+shortString).c_str(),("cnll_"+shortString).c_str(),1800,1800);
//   cZoom = new TCanvas (("cZoom_"+shortString).c_str(),("cZoom_"+shortString).c_str(),1800,1800);
//   cPen = new TCanvas (("cPen_"+shortString).c_str(),("cPen_"+shortString).c_str(),1800,1800);
//   cnll->Divide(3,3);
//   cZoom->Divide(3,3);
//   cPen->Divide(3,3);
// 
//   RooPlot* frame [8];
//   RooPlot* fZoom [8];
//   RooPlot* fPenTerm [8];

//   for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
// 
//     RooRealVar* par = (RooRealVar*)pars.at(iPar);
// 
//     frame[iPar] = par->frame(Name(Form("f1%s",par->GetName())),Title(Form("-log(L) scan vs %s",par->GetTitle()))) ;
//     fZoom[iPar] = par->frame(Name(Form("f2%s",par->GetName())),Title(Form("zoom on -log(L) scan vs %s",par->GetTitle())),
// 			     Range(TMath::Max(par->getMin(),par->getValV()+xZoom*par->getErrorLo()),
// 				   TMath::Min(par->getMax(),par->getValV()+xZoom*par->getErrorHi()) )) ;
// 
//     nll->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;
//     nll->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;
// 
//     if (iPar>0) {
// 
//       double hMax = frame[iPar]->GetMaximum();
// 
//       boundary->plotOn(frame[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
//       boundary->plotOn(fZoom[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
// 
//       if (usedPenalty) {
// 
// 	nll_penalty_plot->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll_penalty_plot->getVal()+10),LineColor(kBlue),LineWidth(2));
// 	nll_penalty_plot->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll_penalty_plot->getVal()+10),LineColor(kBlue),LineWidth(2));
// 
// 	penLog->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));
// 	penLog->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));
// 
//       }
// 
//       frame[iPar]->SetMaximum(hMax);
// 
//       fPenTerm[iPar] = par->frame(Name(Form("f3%s",par->GetName())),Title(Form("Penalty term vs %s",par->GetTitle()))) ;
//       penTerm_eff->plotOn(fPenTerm[iPar],LineColor(4),LineWidth(2)) ;
//       double hMaxP = fPenTerm[iPar]->GetMaximum();
//       boundary->plotOn(fPenTerm[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMaxP,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
//       fPenTerm[iPar]->SetMaximum(hMaxP);
//       cPen->cd(iPar+1);
//       fPenTerm[iPar]->Draw();
// 
//     }
// 
//     fZoom[iPar]->SetMaximum(0.5*xZoom*xZoom);
// 
//     cnll->cd(iPar+1);
//     frame[iPar]->Draw();
// 
//     cZoom->cd(iPar+1);
//     fZoom[iPar]->Draw();
// 
// 
//   }
// 
//   string plotString = shortString + "_" + all_years;
//   if (nSample>0) plotString = plotString + Form("_s%i",nSample);
// 
//   cnll->Update();
//   cnll->SaveAs( ("plotSimFit_d/newphi/recoNLL_scan_" + plotString + ".pdf").c_str() );
// 
//   cZoom->Update();
//   cZoom->SaveAs( ("plotSimFit_d/newphi/recoNLL_scan_" + plotString + "_zoom.pdf").c_str() );
// 
//   cPen->Update();
//   cPen->SaveAs( ("plotSimFit_d/newphi/recoPenTerm_" + plotString + ".pdf").c_str() );
//   return;

  int confIndex = 2*nBins*parity  + q2Bin;
  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,1400);
  c[confIndex]->Divide(3, years.size());  
  
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    RooPlot* xframe = ctK->frame(Title((longString+year).c_str()));
    RooPlot* yframe = ctL->frame(Title((longString+year).c_str()));
    RooPlot* zframe = phi->frame(Title((longString+year).c_str()));
    xframe->GetYaxis()->SetTitleOffset(1.8);
    yframe->GetYaxis()->SetTitleOffset(1.8);
    zframe->GetYaxis()->SetTitleOffset(1.8);
    xframe->SetMaximum(xframe->GetMaximum()*1.15);
    yframe->SetMaximum(yframe->GetMaximum()*1.15);
    zframe->SetMaximum(zframe->GetMaximum()*1.15);
    xframe->SetMinimum(0);
    yframe->SetMinimum(0);
    zframe->SetMinimum(0);
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);

    combData->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year+Form("_subs%d",firstSample)).c_str()), Name(("plData"+year).c_str()));
    combData->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year+Form("_subs%d",firstSample)).c_str()));
    combData->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year+Form("_subs%d",firstSample)).c_str()));

    simPdf->plotOn(xframe,Slice(sample, ("data"+year+Form("_subs%d",firstSample)).c_str()), ProjWData(RooArgSet(sample), *combData), LineWidth(1),Name(("plPDF"+year).c_str()));
    simPdf->plotOn(yframe,Slice(sample, ("data"+year+Form("_subs%d",firstSample)).c_str()), ProjWData(RooArgSet(sample), *combData), LineWidth(1));
    simPdf->plotOn(zframe,Slice(sample, ("data"+year+Form("_subs%d",firstSample)).c_str()), ProjWData(RooArgSet(sample), *combData), LineWidth(1));

    c[confIndex]->cd(iy*3+1); 
    gPad->SetLeftMargin(0.19);        
    xframe->Draw(); 
    leg->Draw("same");
    c[confIndex]->cd(iy*3+2);   
    gPad->SetLeftMargin(0.19); 
    yframe->Draw(); 
    leg->Draw("same");
    c[confIndex]->cd(iy*3+3);    
    gPad->SetLeftMargin(0.19); 
    zframe->Draw(); 
    leg->SetTextSize(0.03);
    leg->AddEntry(xframe->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
    leg->AddEntry(xframe->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
    leg->Draw("same");
  }

  
  c[confIndex]->SaveAs( ("plotSimFit_d/newphi/simFitResult_recoMC_fullAngular_" + shortString + "_" + all_years + stat + "_checkProj_newShapeAnInt.pdf").c_str() );

}


void simfit_recoMC_fullAngularBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool plot, bool save, std::vector<int> years, std::map<int,float> scale_to_data, double fac1, double fac4, double fac5, double base1, double base4, double base5, double max1, double max4, double max5)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullAngularBin(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);
  else
    simfit_recoMC_fullAngularBin(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);
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

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;

  std::map<int,float> scale_to_data;
  // https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
  scale_to_data.insert(std::make_pair(2016, 0.006*2 /2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
  scale_to_data.insert(std::make_pair(2017, 0.005*2 /2.05 ));
  scale_to_data.insert(std::make_pair(2018, 0.007*2 /1.9  ));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_recoMC_fullAngularBin1(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);
  else
    simfit_recoMC_fullAngularBin1(q2Bin, parity, multiSample, nSample, plot, save, years, scale_to_data, fac1, fac4, fac5, base1, base4, base5, max1, max4, max5);

  return 0;

}
