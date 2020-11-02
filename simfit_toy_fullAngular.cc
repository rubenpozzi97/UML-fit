#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <TLine.h>
#include <TRandom3.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooSuperCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>
#include <RooRandom.h>

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

// lower threshold to parameters' uncertainties
// to build the randomisation models (too small leads to many useless points)
double minParError = 0.01;

// MINOS parameters
int nGenMINOS = 1e5;
double widthScale = 0.05;

// Variables to be used both in the main function and the fit subfunc
double coeff1 = 0;
double coeff4 = 0;
double coeff5 = 0;
bool usedPenalty = false;

RooFitResult* fit (RooDataSet* combData, RooAbsPdf* simPdf, RooAbsPdf* simPdf_penalty, RooAbsReal* & nll, RooAbsReal* & nll_penalty, BoundCheck* boundary, Penalty* penTerm, double fac1, double fac4, double base1, double base4, double max1, double max4);
                         
void simfit_toy_fullAngularBin(int q2Bin, vector<double> genPars, uint seed, uint nSample, bool localFiles, bool save, std::vector<int> years, std::map<int,float> scale_to_data, double fac1, double fac4, double base1, double base4, double max1, double max4)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%i",q2Bin);
  cout<<"Conf: "<<shortString<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency

  string effCString = Form("effCHist_b%ip1",q2Bin);
  string effWString = Form("effWHist_b%ip1",q2Bin);
  string intCHistString = "MCint_"+shortString + "p1t1";
  string intWHistString = "MCint_"+shortString + "p1t0";
  string all_years = "";
  string year = ""; 
  string isample = ""; 
  
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
  RooArgSet reco_vars (*ctK, *ctL, *phi, *rand);

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

  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  for (int iPar=0; iPar<pars.getSize(); ++iPar)
    ((RooRealVar*)pars.at(iPar))->setVal(genPars[iPar]);

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    all_years += year;
    for (uint is = 0; is < nSample; is++) {
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

  // check if generation point is physical
  if ( boundary->getValV() > 0 ) {
    cout<<"GEN point is outside physical region. Abort!"<<endl;
    return;
  }

  // Random generators
  RooRandom::randomGenerator()->SetSeed(seed);
  TRandom3 randGen (seed);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("recoMCDataset_b%i_%i.root", q2Bin, years[iy]);
    if (!localFiles) filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/", years[iy]) + filename_data;

    // import data (or MC as data proxy)
    fin_data.push_back( TFile::Open( filename_data.c_str() ) );
    if ( !fin_data[iy] || !fin_data[iy]->IsOpen() ) {
      cout << "File not found: " << filename_data << endl;
      return;
    }
    wsp.push_back( (RooWorkspace*)fin_data[iy]->Get(Form("ws_b%ip0", q2Bin ) ) );
    if ( !wsp[iy] || wsp[iy]->IsZombie() ) {
      cout<<"Workspace not found in file: "<<filename_data<<endl;
      return;
    }
  

    // import KDE efficiency histograms and partial integral histograms
    string filename = Form("KDEeff_b%i_od_%i.root",q2Bin,years[iy]);
    if (!localFiles) filename = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/",years[iy]) + filename;
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



    // define angular PDF for signal, using the custom class
    // efficiency function and integral values are passed as arguments
    PDF_sig_ang_fullAngular.push_back( new PdfSigAng(("PDF_sig_ang_fullAngular_"+shortString+"_"+year).c_str(),
                                                     ("PDF_sig_ang_fullAngular_"+year).c_str(),
      		                                     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
      		                                     *effC[iy], *effW[iy], intCVec[iy],intWVec[iy]));
    // define PDF with penalty term
    PDF_sig_ang_fullAngular_penalty.push_back( new PdfSigAng(("PDF_sig_ang_fullAngular_penalty_"+shortString+"_"+year).c_str(),
							     ("PDF_sig_ang_fullAngular_penalty_"+year).c_str(),
							     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
							     *effC[iy], *effW[iy], intCVec[iy],intWVec[iy],*penTerm));

    // create roodataset (in case data-like option is selected, only import the correct % of data)
    std::vector<RooDataSet*> data_isample;

    int numCT = scale_to_data[years[iy]] * ((RooDataSet*)wsp[iy]->data(Form("data_ctRECO_ev_b%i",q2Bin)))->numEntries() ;
    int numWT = scale_to_data[years[iy]] * ((RooDataSet*)wsp[iy]->data(Form("data_wtRECO_ev_b%i",q2Bin)))->numEntries() ;
    if ( numCT==0 || numWT==0 ) {
      cout<<"Empty dataset found in file: "<<filename_data<<endl;
      return;
    }

    for (uint is = 0; is < nSample; is++) {
      
      RooDataSet* datatmp = PDF_sig_ang_fullAngular.back()->generate(vars,numCT+numWT,Name(("data_"+shortString + Form("_subs%i", is)).c_str()));

      datatmp->Print();
      data_isample.push_back (datatmp);
    }

    data.push_back(data_isample) ;

    // insert sample in the category map, to be imported in the combined dataset
    // and associate model with the data
    for (uint is = 0; is < nSample; is++) {
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
      simPdf->addPdf(*PDF_sig_ang_fullAngular[iy], ("data"+year+Form("_subs%d",is)).c_str());
      simPdf_penalty->addPdf(*PDF_sig_ang_fullAngular_penalty[iy], ("data"+year+Form("_subs%d",is)).c_str());
    }

  }

  string fout_name = "toyFitResults/simFitResult_toy2_fullAngular_" + all_years + Form("_b%i_",q2Bin);
  for (int iPar=0; iPar<pars.getSize(); ++iPar)
    fout_name = fout_name + Form((iPar>0?"-%.3f":"%.3f"),genPars[iPar]);
  fout_name = fout_name + Form("_s%i.root",seed);
  TFile* fout = new TFile(fout_name.c_str(),"UPDATE");

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
  RooRealVar* minTime = new RooRealVar("minTime","MINOS time",0,"s");
  RooRealVar* co1 = new RooRealVar("co1","Coefficient 1",0);
  RooRealVar* co4 = new RooRealVar("co4","Coefficient 4",0);
  RooRealVar* co5 = new RooRealVar("co5","Coefficient 5",0);
  RooRealVar* boundDist = new RooRealVar("boundDist","Distance from boundary",0);
  RooArgSet savePars (*co1,*co4,*co5,*fitTime,*minTime,*boundDist);
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

  // TTree with the MINOS output
  vector<double> vConfInterLow  (pars.getSize());
  vector<double> vConfInterHigh (pars.getSize());
  vector<double> vFitResult  (pars.getSize());
  vector<double> vFitErrLow  (pars.getSize());
  vector<double> vFitErrHigh (pars.getSize());
  fout->cd();
  TTree* MINOS_output = (TTree*)fout->Get("MINOS_output");
  if (!MINOS_output || MINOS_output->IsZombie()) {
    MINOS_output = new TTree("MINOS_output","MINOS_output");
    for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
      RooRealVar* par = (RooRealVar*)pars.at(iPar);
      MINOS_output->Branch(Form("%s_low",par->GetName()),&vConfInterLow[iPar]);
      MINOS_output->Branch(Form("%s_high",par->GetName()),&vConfInterHigh[iPar]);
      MINOS_output->Branch(Form("%s_best",par->GetName()),&vFitResult[iPar]);
    }
  } else {
    for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
      RooRealVar* par = (RooRealVar*)pars.at(iPar);
      MINOS_output->SetBranchAddress(Form("%s_low",par->GetName()),&vConfInterLow[iPar]);
      MINOS_output->SetBranchAddress(Form("%s_high",par->GetName()),&vConfInterHigh[iPar]);
      MINOS_output->SetBranchAddress(Form("%s_best",par->GetName()),&vFitResult[iPar]);
    }
  }

  // Timer for fitting time
  TStopwatch subTime;

  // counters to monitor results' status
  int cnt[9];
  for (int iCnt=0; iCnt<9; ++iCnt) cnt[iCnt] = 0;

  for (uint is = 0; is < nSample; is++) {

    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }

    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
    cout<<"Fitting toy "<<is<<" with "<<combData->numEntries()<<" entries"<<endl;

    // set penalty term power parameter
    int combEntries = combData->numEntries();
    penTerm->setPower(power/combEntries);

    double base1_corr = base1*sqrt(combEntries);
    double base4_corr = base4*sqrt(combEntries);
    if (base1_corr<min_base) base1_corr = min_base;
    if (base4_corr<min_base) base4_corr = min_base;

    // to start the fit, parameters are restored to the center of the parameter space
    Fl ->setVal(0.5);
    P1 ->setVal(0);
    P2 ->setVal(0);
    P3 ->setVal(0);
    P4p->setVal(0);
    P5p->setVal(0);
    P6p->setVal(0);
    P8p->setVal(0);

    // run the fit
    subTime.Start(true);
    RooFitResult* fitResult = fit(combData,simPdf,simPdf_penalty,nll,nll_penalty,boundary,penTerm,fac1,fac4,base1_corr,base4_corr,max1,max4);
    subTime.Stop();

    // include fit time in dataset with per-toy informations
    fitTime->setVal(subTime.CpuTime());
    // fitTime->setVal(subTime.RealTime());

    co1->setVal(0);
    co4->setVal(0);
    co5->setVal(0);
    boundDist->setVal(0);
    minTime->setVal(0);

    bool convCheck = false;
    bool boundCheck = false;

    if (fitResult) {
      
      convCheck = true;
      boundCheck = boundary->getValV() == 0;

      if (usedPenalty) {
	// include coefficient values in dataset with per-toy informations
	co1->setVal(coeff1);
	co4->setVal(coeff4);
	co5->setVal(coeff5);
      }

      // Compute distance from boundary, print it
      // and save it in dataset with per-toy informations
      double boundDistVal = bound_dist->getValV();
      boundDist->setVal(boundDistVal);

      // Improve global fit result
      if (usedPenalty) {
	vector<double> vTestPar(pars.getSize());
	vector<double> vImprovPar(pars.getSize());
	for (int iPar = 0; iPar < pars.getSize(); ++iPar)
	  vImprovPar[iPar] = ((RooRealVar*)pars.at(iPar))->getValV();
	double NLL_before = nll->getValV();
	double improvNLL = NLL_before;
	double testNLL = 0;
	int iImprove = 0;
	do {
	  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
	    RooRealVar* par = (RooRealVar*)pars.at(iPar);
	    do vTestPar[iPar] = randGen.Gaus(vImprovPar[iPar],TMath::Max(boundDistVal,0.002));
	    while (vTestPar[iPar]>par->getMax() || vTestPar[iPar]<par->getMin());
	    par->setVal(vTestPar[iPar]);
	  }
	  if (boundary->getValV()>0) continue;
	  testNLL = nll->getValV();
	  if (improvNLL>testNLL) {
	    improvNLL = testNLL;
	    for (int iPar = 0; iPar < pars.getSize(); ++iPar)
	      vImprovPar[iPar] = vTestPar[iPar];
	  }
	  ++iImprove;
	} while (iImprove<1e4);
	for (int iPar = 0; iPar < pars.getSize(); ++iPar)
	  ((RooRealVar*)pars.at(iPar))->setVal(vImprovPar[iPar]);

	// double improvDistVal = bound_dist->getValV();
	// cout<<"Improved fit result: deltaNLL = "<<NLL_before-improvNLL<<" bound dist: "<<boundDistVal<<" -> "<<improvDistVal<<endl;
      }

      // run MINOS error
      TStopwatch minosTime;
      minosTime.Start(true);

      // NLL of the best-fit result
      double NLL_min = nll->getValV();

      double probedNLL;

      // get best-fit results and errors from the fit
      for (int iPar = 0; iPar < pars.getSize(); ++iPar) {

	RooRealVar* par = (RooRealVar*)pars.at(iPar);
	vFitResult [iPar] = par->getValV();
	vFitErrLow [iPar] = par->getErrorLo();
	vFitErrHigh[iPar] = par->getErrorHi();

      }

      // Loop over the parameters
      for (int iPar = 0; iPar < pars.getSize(); ++iPar) {

	RooRealVar* par = (RooRealVar*)pars.at(iPar);

	// get and print the best-fit result
	double p_best = vFitResult[iPar];

	// firstly low, then high error
	for (int isErrHigh=0; isErrHigh<2; ++isErrHigh) {

	  vector<double> vLastHit(0);
	  for (int iPar1 = 0; iPar1 < pars.getSize(); ++iPar1)
	    vLastHit.push_back(vFitResult[iPar1]);

	  TH1D* parRandomPool = 0;
	  int nHistBins = 0;
	  if (isErrHigh>0) {
	    nHistBins = (int)((par->getMax()-p_best)/0.0005);
	    // nHistBins = (int)((par->getMax()-p_best+0.0005)/0.001);
	    if (nHistBins<1) nHistBins=1;
	    parRandomPool = new TH1D(Form("hRandPoolH%i%i",iPar,is),
				     Form("hRandPoolH%i%i",iPar,is),
				     nHistBins,par->getMax()-0.0005*nHistBins,par->getMax());
	    // nHistBins,par->getMax()+0.0005-0.001*nHistBins,par->getMax()+0.0005);
	  } else {
	    nHistBins = (int)((p_best-par->getMin())/0.0005);
	    // nHistBins = (int)((p_best-par->getMin()+0.0005)/0.001);
	    if (nHistBins<1) nHistBins=1;
	    parRandomPool = new TH1D(Form("hRandPoolL%i%i",iPar,is),
				     Form("hRandPoolL%i%i",iPar,is),
				     nHistBins,par->getMin(),par->getMin()+0.0005*nHistBins);
	    // nHistBins,par->getMin()-0.0005,par->getMin()-0.0005+0.001*nHistBins);
	  }
	  double sigma = fabs(isErrHigh>0?vFitErrHigh[iPar]:vFitErrLow[iPar]);
	  if (sigma<minParError) sigma = minParError;
	  for (int iBin=1; iBin<=nHistBins; ++iBin) {
	    double x = (parRandomPool->GetBinCenter(iBin)-p_best) / sigma;
	    parRandomPool->SetBinContent(iBin,fabs(x)*exp(-0.5*x*x));
	  }

	  double p_in = p_best;
	  double p_test = 0;

	  int iPnt=0;

	  do {

	    do p_test = parRandomPool->GetRandom();
	    while (p_test>par->getMax() || p_test<par->getMin());
	    par->setVal(p_test);
	    for (int iPar1 = 0; iPar1 < pars.getSize(); ++iPar1) {
	      if (iPar1==iPar) continue;
	      RooRealVar* par1 = (RooRealVar*)pars.at(iPar1);
	      double par1val = 0;
	      do par1val = randGen.Gaus(vLastHit[iPar1],widthScale*TMath::Max(0.5*(vFitErrHigh[iPar1]-vFitErrLow[iPar1]),minParError));
	      while (par1val>par1->getMax() || par1val<par1->getMin());
	      par1->setVal(par1val);
	    }
	    // check if the point is physical
	    if (boundary->getValV()>0) continue;
	    // get and test the local likelihood
	    probedNLL = nll->getValV();
	    if (probedNLL<=NLL_min+0.5) {
	      p_in = p_test;
	      if ( isErrHigh > 0 ) { if ( p_in > par->getMax()-parRandomPool->GetBinWidth(1) ) break; }
	      else if ( p_in < par->getMin()+parRandomPool->GetBinWidth(1) ) break;
	      for (int iPar1 = 0; iPar1 < pars.getSize(); ++iPar1) {
		RooRealVar* par1 = (RooRealVar*)pars.at(iPar1);
		vLastHit[iPar1] = par1->getValV();
	      }
	      if (isErrHigh>0)
		for (int iBin=1; iBin<=parRandomPool->FindBin(p_test); ++iBin)
		  parRandomPool->SetBinContent(iBin,0);
	      else 
		for (int iBin=parRandomPool->FindBin(p_test); iBin<=nHistBins; ++iBin)
		  parRandomPool->SetBinContent(iBin,0);
	    } else {
	      parRandomPool->Fill(p_test,0.02/(probedNLL-NLL_min-0.5));
	    }
	  
	    ++iPnt;
	    // apply conditions
	  } while ( iPnt < nGenMINOS );

	  if (isErrHigh>0) vConfInterHigh[iPar] = p_in;
	  else vConfInterLow[iPar] = p_in;

	}

      }

      minosTime.Stop();
      minTime->setVal(minosTime.CpuTime());

      // save MINOS errors
      MINOS_output->Fill();

      cout<<"CPU time: "<<subTime.CpuTime()<<"\t"<<minosTime.CpuTime()<<endl;
    
    }

    // save fit status and times
    if (!boundCheck)
      if (convCheck) subNegConv->add(savePars);
      else subNegNotc->add(savePars);
    else
      if (convCheck)
	if (usedPenalty) subPosConv->add(savePars);
	else subNoPen->add(savePars);
      else subPosNotc->add(savePars);

    // fill fit-status-dependent counters
    ++cnt[8];
    int iCnt = 0;
    if (!convCheck) iCnt += 4;
    if (!boundCheck) iCnt += 2;
    if (usedPenalty) iCnt += 1;
    ++cnt[iCnt];

  }

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

  if (save) {
    RooWorkspace* wksp = new RooWorkspace("wsToy","Workspace with set of toy fit results");

    wksp->import(*subResults);

    fout->cd();
    wksp->Write();
    MINOS_output->Write();

  }

  fout->Close();

}


int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = -1;

  if ( argc > 1 ) q2Bin   = atoi(argv[1]);

  vector<double> genPars(8);
  for (int iPar=0; iPar<8; ++iPar)
    if ( argc > 2+iPar ) genPars[iPar] = atof(argv[2+iPar]);
    else genPars[iPar] = 0;

  double fac1 = 1;
  double fac4 = 1;
  double base1 = 3;
  double base4 = 3;
  double max1 = 0;
  double max4 = 0;

  if ( argc > 10 ) fac1  = atof(argv[10]) / 1000.0;
  if ( argc > 11 ) fac4  = atof(argv[11]) / 1000.0;
  if ( argc > 12 ) base1 = atof(argv[12]) / 1000.0;
  if ( argc > 13 ) base4 = atof(argv[13]) / 1000.0;
  if ( argc > 14 ) max1  = atof(argv[14]);
  if ( argc > 15 ) max4  = atof(argv[15]);

  uint seed = 1;
  uint nSample = 1;
  if ( argc > 16 ) seed = atoi(argv[16]);
  if ( argc > 17 ) nSample = atoi(argv[17]);

  bool localFiles = false;
  if ( argc > 18 && atoi(argv[18]) > 0 ) localFiles = true;

  bool save = true;

  if ( argc > 19 && atoi(argv[19]) == 0 ) save = false;

  std::vector<int> years;
  if ( argc > 20 && atoi(argv[20]) != 0 ) years.push_back(atoi(argv[20]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 21 && atoi(argv[21]) != 0 ) years.push_back(atoi(argv[21]));
  if ( argc > 22 && atoi(argv[22]) != 0 ) years.push_back(atoi(argv[22]));

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;

  std::map<int,float> scale_to_data;
  // https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
  scale_to_data.insert(std::make_pair(2016, 0.006*2 /2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
  scale_to_data.insert(std::make_pair(2017, 0.005*2 /2.05 ));
  scale_to_data.insert(std::make_pair(2018, 0.007*2 /1.9  ));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_toy_fullAngularBin(q2Bin, genPars, seed, nSample, localFiles, save, years, scale_to_data, fac1, fac4, base1, base4, max1, max4);
  else
    simfit_toy_fullAngularBin(q2Bin, genPars, seed, nSample, localFiles, save, years, scale_to_data, fac1, fac4, base1, base4, max1, max4);

  return 0;

}

RooFitResult* fit (RooDataSet* combData,
		   RooAbsPdf* simPdf,
		   RooAbsPdf* simPdf_penalty,
		   RooAbsReal* & nll,
		   RooAbsReal* & nll_penalty,
		   BoundCheck* boundary,
		   Penalty* penTerm,
		   double fac1,
		   double fac4,
		   double base1,
		   double base4,
		   double max1,
		   double max4)
{

    coeff1 = 0;
    coeff4 = 0;
    coeff5 = 0;

    // set up free fit
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

    // free fit
    m.setStrategy(0);
    // m.setEvalErrorWall(false);
    m.migrad() ;
    m.hesse() ;
    // std::cout << std::endl;
    // std::cout << "######################### now strategy 2 #########################"<< std::endl;
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;
    m.minos() ;
    
    RooFitResult* fitResult = m.save("result") ;

    RooFitResult* fitResult_penalty = 0;
    usedPenalty = false;

    // if free fit is good return its result
    if ( fitResult->status()==0 && fitResult->covQual()==3 && boundary->getValV() == 0 ) return fitResult;

    usedPenalty = true;

    // optional: if a partial boundary is satisfied
    // do not apply the corresponding penalty term
    // bool inCTL4  = true;
    // bool inCTL15 = true;
    // if ( !boundary->isInCTL4()  ) inCTL4  = false;
    // if ( !boundary->isInCTL15() ) inCTL15 = false;

    for (int totCoeff=0; fac1*pow(base1,totCoeff)<=maxCoeff; ++totCoeff) {

      for (int iCoeff1=totCoeff; iCoeff1>=0; --iCoeff1) {

	// set penalty coefficients
	coeff1 = fac1 * pow(base1,iCoeff1);
	if (max1>0 && coeff1>max1) continue;

	coeff4 = fac4 * pow(base4,totCoeff-iCoeff1);
	if (max4>0 && coeff4>max4) continue;

	coeff5 = pow(coeff1,1.5) / 316.2;

	// optional: if a partial boundary is satisfied
	// do not apply the corresponding penalty term
	// if ( inCTL15 ) {
	//   if ( iCoeff1>0 ) continue;
	//   coeff1 = 0;
	//   coeff5 = 0;
	// }
	// if ( inCTL4 ) {
	//   if ( totCoeff-iCoeff1>0 ) continue;
	//   coeff4 = 0;
	// }

	penTerm->setCoefficient(1,coeff1);
	penTerm->setCoefficient(4,coeff4);
	penTerm->setCoefficient(5,coeff5);

	// set up the penalised fit
	nll_penalty = simPdf_penalty->createNLL(*combData,
						RooFit::Extended(kFALSE),
						RooFit::NumCPU(1)
						);

	RooMinimizer m_penalty (*nll_penalty) ;
	m_penalty.optimizeConst(kTRUE);
	m_penalty.setOffsetting(kTRUE);
	// m_penalty.setVerbose(kTRUE);
	m_penalty.setMinimizerType("Minuit2");
	// m_penalty.setProfile(kTRUE);
	m_penalty.setPrintLevel(-1);
	m_penalty.setPrintEvalErrors(-1);
	m_penalty.setStrategy(2);
    
	// penalised fit
	m_penalty.migrad() ;
	m_penalty.hesse() ;
	fitResult_penalty = m_penalty.save("result");
	    
	// cout<<penTerm->getCoefficient(1)<<"\t"<<penTerm->getCoefficient(5)<<"\t"<<P5p->getValV()<<endl;
	// fitResult_penalty->Print("v");

	// if a good fit is found return its result
	if ( fitResult_penalty->status()==0 && fitResult_penalty->covQual()==3 ) {
	  if ( boundary->getValV()==0 ) {
	    // cout<<"P "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
	    return fitResult_penalty;
	  } // else cout<<"O "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
	} // else cout<<"N "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;

      }

    }
    
    // if no good fit is found return a null pointer
    return 0;

}
