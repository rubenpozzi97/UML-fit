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
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>
#include <RooAddition.h>  /// check if needed
#include <RooRandom.h>

#include "ShapeSigAng.h"
#include "PdfSigAng.h"
#include "BoundCheck.h"
#include "BoundDist.h"
#include "Penalty.h"
#include "Fitter.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* cPen;
TCanvas* c [4*nBins];

double power = 1.0;

void simfit_recoMC_fullAngularBin(int q2Bin, int parity, bool multiSample, uint nSample, bool localFiles, bool plot, bool save, std::vector<int> years, std::map<int,float> scale_to_data)
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

  // Random generators
  RooRandom::randomGenerator()->SetSeed(1);

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
    wsp.push_back( (RooWorkspace*)fin_data[iy]->Get(Form("ws_b%ip%i", q2Bin, 1-parity ) ) );
    if ( !wsp[iy] || wsp[iy]->IsZombie() ) {
      cout<<"Workspace not found in file: "<<filename_data<<endl;
      return;
    }
  

    // import KDE efficiency histograms and partial integral histograms
    string filename = Form((parity==0 ? "KDEeff_b%i_ev_%i.root" : "KDEeff_b%i_od_%i.root"),q2Bin,years[iy]);
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
  RooArgSet *params      = (RooArgSet*) simPdf->getParameters(vars);
  RooArgSet* savedParams = (RooArgSet*) params->snapshot() ;

  // Results' containers
  RooRealVar* fitTime = new RooRealVar("fitTime","fit time",0,"s");
  RooRealVar* minTime = new RooRealVar("minTime","MINOS time",0,"s");
  RooRealVar* co1 = new RooRealVar("co1","Coefficient 1",0);
  RooRealVar* co4 = new RooRealVar("co4","Coefficient 4",0);
  RooRealVar* co5 = new RooRealVar("co5","Coefficient 5",0);
  RooRealVar* boundDist = new RooRealVar("boundDist","Distance from boundary",0);
  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
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
  vector<double> vFitResult  (pars.getSize());
  vector<double> vConfInterLow  (pars.getSize());
  vector<double> vConfInterHigh (pars.getSize());
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

  Fitter* fitter = 0;
  vector<Fitter*> vFitter (0);

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

    // set penalty term power parameter
    int combEntries = combData->numEntries();
    penTerm->setPower(power/combEntries);

//     for(auto it = map.cbegin(); it != map.cend(); ++it)
//       std::cout << "dataset: " << it->first << ", with n entries: " << it->second->sumEntries() << "\n";

    // to start the fit, parameters are restored to the center of the parameter space
//     Fl ->setVal(0.5);
//     P1 ->setVal(0);
//     P2 ->setVal(0);
//     P3 ->setVal(0);
//     P4p->setVal(0);
//     P5p->setVal(0);
//     P6p->setVal(0);
//     P8p->setVal(0);
    *params = *savedParams ;

    // run the fit
    fitter = new Fitter (Form("fitter%i",is),Form("fitter%i",is),pars,combData,simPdf,simPdf_penalty,boundary,bound_dist,penTerm);
    vFitter.push_back(fitter);
    subTime.Start(true);
    int status = fitter->fit();
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

    if (status==0) {
      
      convCheck = true;
      boundCheck = boundary->getValV() == 0;

      fitter->result()->SetName (Form("result_%s_subs%i",shortString.c_str(),is));
      fitter->result()->SetTitle(Form("result_%s_subs%i",shortString.c_str(),is));
      fitter->result()->Print("v");

      if (fitter->usedPenalty) {
	// include coefficient values in dataset with per-toy informations
	co1->setVal(fitter->coeff1);
	co4->setVal(fitter->coeff4);
	co5->setVal(fitter->coeff5);

	// Compute distance from boundary, print it
	// and save it in dataset with per-toy informations
	TStopwatch distTime;
	distTime.Start(true);
	double boundDistVal = bound_dist->getValV();
	distTime.Stop();
	cout<<"Distance from boundary: "<<boundDistVal<<" (computed in "<<distTime.CpuTime()<<" s)"<<endl;
	boundDist->setVal(boundDistVal);

	fitter->improveAng();
      }

      if (nSample>0) {
	// run MINOS error
	TStopwatch minosTime;
	minosTime.Start(true);

	fitter->MinosAng();

	minosTime.Stop();
	minTime->setVal(minosTime.CpuTime());
	cout<<"MINOS errors computed in "<<minosTime.CpuTime()<<" s"<<endl;

	// cout<<"Error difference [custMINOS - fit], lower and higher:"<<endl;
	// for (int iPar = 0; iPar < pars.getSize(); ++iPar)
	// 	cout<<vFitResult[iPar]-vConfInterLow[iPar]+vFitErrLow[iPar]<<"   \t"
	// 	    <<vConfInterHigh[iPar]-vFitResult[iPar]-vFitErrHigh[iPar]<<endl;

	// save MINOS errors
	for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
	  vFitResult[iPar] = fitter->vFitResult[iPar];
	  vConfInterLow[iPar] = fitter->vConfInterLow[iPar];
	  vConfInterHigh[iPar] = fitter->vConfInterHigh[iPar];
	}
	MINOS_output->Fill();
      }

    }

    // fill fit-status-dependent counters
    ++cnt[8];
    int iCnt = 0;
    if (!convCheck) iCnt += 4;
    if (!boundCheck) iCnt += 2;
    if (fitter->usedPenalty) iCnt += 1;
    ++cnt[iCnt];

    // print fit status and time
    if (!boundCheck) {
      if (convCheck) {
	subNegConv->add(savePars);
	cout<<"Converged in unphysical region";
      } else {
	subNegNotc->add(savePars);
	cout<<"Not converged";
      }
    } else {
      if (convCheck) {
	if (fitter->usedPenalty) {
	  subPosConv->add(savePars);
	  cout<<"Converged with penalty term with coeff: "<<fitter->coeff1<<" "<<fitter->coeff4<<" "<<fitter->coeff5;
	} else {
	  subNoPen->add(savePars);
	  cout<<"Converged without penalty";
	}
      } else {
	subPosNotc->add(savePars);
	cout<<"This should never be printed";
      }
    }
    cout<<" ("<<fitTime->getValV()<<"s)"<<endl;

    // Save fit results in file
    if (save && convCheck) {
      fout->cd();
      fitter->result()->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
    }

  }  

  if (multiSample) {
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

  if (save) {
    RooWorkspace* wksp = new RooWorkspace(((multiSample?"wsMulti_":"ws_")+shortString+Form("_s%i_pow%.1f",nSample,power)).c_str(),
					  (multiSample?"Workspace with set of RECO subsample fit results":
					   (nSample>0?"Workspace with RECO subsample fit result":
					    "Workspace with full RECO fit result")));

    if (multiSample) {
      wksp->import(*subResults);
    } else {
      wksp->import(*combData,Rename("data"));
      wksp->import(*simPdf,RenameVariable(simPdf->GetName(),"pdf"),Silence());
      if (fitter->usedPenalty) {
	wksp->import(*simPdf_penalty,RenameVariable(simPdf_penalty->GetName(),"pdfPen"),Silence(),RecycleConflictNodes());
	wksp->import(*penTerm,Silence(),RecycleConflictNodes());
      }

    }

    fout->cd();
    wksp->Write();
    if (nSample>0) MINOS_output->Write();

  }

  fout->Close();

  if (multiSample) {
    TCanvas* cDist = new TCanvas (("cDist_"+shortString).c_str(),("cDist_"+shortString).c_str(),1800,1800);
    RooPlot* fDist = boundDist->frame(Name("fDist"),Title("Distribution of results' distance fram boundary"),Range(0,0.1));
    subNoPen->plotOn(fDist,Binning(50,0,0.1),LineColor(kBlue),MarkerColor(kBlue),MarkerStyle(19),DrawOption("XL"));
    subPosConv->plotOn(fDist,Binning(50,0,0.1),LineColor(kRed),MarkerColor(kRed),MarkerStyle(19),DrawOption("XL"));
    cDist->cd();
    fDist->Draw();
    cDist->SaveAs( ("plotSimFit_d/recoBoundDist_" + shortString + "_" + all_years + ".pdf").c_str() );
  }

  if (!plot || multiSample) return;

  // For plotting the effective penalty term is used
  Penalty* penTerm_eff = new Penalty(*penTerm,"penTerm_eff");
  penTerm_eff->setPower(power);
  RooFormulaVar* penLog = new RooFormulaVar("penLog","penLog","-1.0 * log(penTerm_eff)",RooArgList(*penTerm_eff));

  double xZoom = 200.0;
  if (nSample>0) xZoom = 2.0;

  cnll  = new TCanvas (("cnll_"+shortString).c_str(),("cnll_"+shortString).c_str(),1800,1800);
  cZoom = new TCanvas (("cZoom_"+shortString).c_str(),("cZoom_"+shortString).c_str(),1800,1800);
  cPen = new TCanvas (("cPen_"+shortString).c_str(),("cPen_"+shortString).c_str(),1800,1800);
  cnll->Divide(3,3);
  cZoom->Divide(3,3);
  cPen->Divide(3,3);

  RooPlot* frame [8];
  RooPlot* fZoom [8];
  RooPlot* fPenTerm [8];

  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)pars.at(iPar);

    frame[iPar] = par->frame(Name(Form("f1%s",par->GetName())),Title(Form("-log(L) scan vs %s",par->GetTitle()))) ;
    fZoom[iPar] = par->frame(Name(Form("f2%s",par->GetName())),Title(Form("zoom on -log(L) scan vs %s",par->GetTitle())),
			     Range(TMath::Max(par->getMin(),par->getValV()+xZoom*par->getErrorLo()),
				   TMath::Min(par->getMax(),par->getValV()+xZoom*par->getErrorHi()) )) ;

    fitter->nll->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;
    fitter->nll->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;

    if (iPar>0) {

      double hMax = frame[iPar]->GetMaximum();

      boundary->plotOn(frame[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
      boundary->plotOn(fZoom[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));

      if (fitter->usedPenalty) {

	fitter->nll_penalty->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll_penalty->getVal()+10),LineColor(kBlue),LineWidth(2));
	fitter->nll_penalty->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(fitter->nll_penalty->getVal()+10),LineColor(kBlue),LineWidth(2));

	penLog->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));
	penLog->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));

      }

      frame[iPar]->SetMaximum(hMax);

      fPenTerm[iPar] = par->frame(Name(Form("f3%s",par->GetName())),Title(Form("Penalty term vs %s",par->GetTitle()))) ;
      penTerm_eff->plotOn(fPenTerm[iPar],LineColor(4),LineWidth(2)) ;
      double hMaxP = fPenTerm[iPar]->GetMaximum();
      boundary->plotOn(fPenTerm[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMaxP,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
      fPenTerm[iPar]->SetMaximum(hMaxP);
      cPen->cd(iPar+1);
      fPenTerm[iPar]->Draw();

    }

    fZoom[iPar]->SetMaximum(0.5*xZoom*xZoom);

    cnll->cd(iPar+1);
    frame[iPar]->Draw();

    cZoom->cd(iPar+1);
    fZoom[iPar]->Draw();


  }

  string plotString = shortString + "_" + all_years;
  if (nSample>0) plotString = plotString + Form("_s%i",nSample);

  cnll->Update();
  cnll->SaveAs( ("plotSimFit_d/recoNLL_scan_" + plotString + ".pdf").c_str() );

  cZoom->Update();
  cZoom->SaveAs( ("plotSimFit_d/recoNLL_scan_" + plotString + "_zoom.pdf").c_str() );

  cPen->Update();
  cPen->SaveAs( ("plotSimFit_d/recoPenTerm_" + plotString + ".pdf").c_str() );


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


void simfit_recoMC_fullAngularBin1(int q2Bin, int parity, bool multiSample, uint nSample, bool localFiles, bool plot, bool save, std::vector<int> years, std::map<int,float> scale_to_data)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullAngularBin(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years, scale_to_data);
  else
    simfit_recoMC_fullAngularBin(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years, scale_to_data);
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

  bool multiSample = false;
  uint nSample = 0;
  if ( argc > 3 && atoi(argv[3]) > 0 ) multiSample = true;
  if ( argc > 4 ) nSample = atoi(argv[4]);

  if (nSample==0) multiSample = false;

  bool localFiles = false;
  if ( argc > 5 && atoi(argv[5]) > 0 ) localFiles = true;

  bool plot = true;
  bool save = true;

  if ( argc > 6 && atoi(argv[6]) == 0 ) plot = false;
  if ( argc > 7 && atoi(argv[7]) == 0 ) save = false;

  std::vector<int> years;
  if ( argc > 8 && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));
  else {
    cout << "No specific years selected, using default: 2016" << endl;
    years.push_back(2016);
  }
  if ( argc > 9  && atoi(argv[9])  != 0 ) years.push_back(atoi(argv[9]));
  if ( argc > 10 && atoi(argv[10]) != 0 ) years.push_back(atoi(argv[10]));

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
      simfit_recoMC_fullAngularBin1(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years, scale_to_data);
  else
    simfit_recoMC_fullAngularBin1(q2Bin, parity, multiSample, nSample, localFiles, plot, save, years, scale_to_data);

  return 0;

}
