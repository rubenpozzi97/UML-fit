#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH1C.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooMinimizer.h>
#include <RooCategory.h>

#include "DecayRate.h"
#include "DecayRate_Pen.h"
#include "BoundCheck.h"
#include "Penalty.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cTime;

double genStat  [nBins] = {1700246,3190316,2416720,4413921,1,5184529,1,3094894,1};
double dataStat [nBins] = {994,2288,1979,3403,0,5381,0,3506,0};
double d16Stat  [nBins] = {246,605,519,873,0,1301,0,908,0};

double fac1 = 10;
double fac4 = 10;
double fac5 = 10;
double base1 = 3;
double base4 = 3;
double base5 = 3;
double maxCoeff = 1e8;

void fit_genMCBin(int q2Bin, int parity, bool plot, bool save, int nSamples, double power, double coeff1, double coeff4, double coeff5)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  int stat = nSamples / 10000;
  nSamples = nSamples % 10000;
  double statFactor = 0;
  switch(stat) {
  case 1:  statFactor = d16Stat [q2Bin]/genStat[q2Bin]; break;
  default: statFactor = dataStat[q2Bin]/genStat[q2Bin]; break;
  }
  
  if ( nSamples>1 ) {
    if ( stat > -1 && stat < 2 )
      cout<<"Fitting "<<nSamples<<Form(" subsamples with about 1/%.0f events wrt full GEN sample",1.0/statFactor)<<endl;
    else {
      cout<<"Error: nSamples / 10000 = "<<stat<<" not valid."<<endl;
      return;
    }
  } else {
    cout<<"Error: nSamples % 10000 = "<<nSamples<<" not valid."<<endl;
    return;
  }

  // Load variables and dataset
  // import the "other parity" dataset, to stay coherent with fit_recoMC notation
  // gen without cuts is the same for the 3 years
  string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/2016/lmnr/effDataset_GENOnly_b%i.root",q2Bin);
  TFile* fin_data = TFile::Open( filename_data.c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: "<<filename_data<<endl;
    return;
  }
  // import the "other parity" dataset, to stay coherent with fit_recoMC notation
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(Form("ws_b%ip%i",q2Bin,1-parity));
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace not found in file: "<<filename_data<<endl;
    return;
  }
  RooRealVar* ctK = wsp->var("ctK");
  RooRealVar* ctL = wsp->var("ctL");
  RooRealVar* phi = wsp->var("phi");
  RooRealVar* rand = wsp->var("rand");
  if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
    cout<<"Variables not found in file: "<<filename_data<<endl;
    return;
  }
  RooArgList vars (* ctK,* ctL,* phi);
  RooArgSet varsPlusOne (*ctK, *ctL, *phi, *rand);
  string datasetString = Form((parity==1?"data_genDen_ev_b%i":"data_genDen_od_b%i"),q2Bin);
  RooDataSet* fullData = (RooDataSet*)wsp->data(datasetString.c_str());
  if ( !fullData || fullData->IsZombie() ) {
    cout<<"Dataset "<<datasetString<<" not found in file: "<<filename_data<<endl;
    return;
  }
  RooDataSet* data = 0;

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl = new RooRealVar("Fl","F_{L}",0.5,0.0,1.0);
  RooRealVar* P1 = new RooRealVar("P1","P_{1}",0.0,-1.0,1.0);
  RooRealVar* P2 = new RooRealVar("P2","P_{2}",0.0,-0.5,0.5);
  RooRealVar* P3 = new RooRealVar("P3","P_{3}",0.0,-0.5,0.5);
  RooRealVar* P4p = new RooRealVar("P4p","P'_{4}",0.0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p = new RooRealVar("P5p","P'_{5}",0.0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p = new RooRealVar("P6p","P'_{6}",0.0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p = new RooRealVar("P8p","P'_{8}",0.0,-1*sqrt(2),sqrt(2));

  // Physical region representation
  BoundCheck* boundary = new BoundCheck("boundary","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // Define penalty term
  Penalty* penTerm = new Penalty("penTerm","Penalty term",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,power,coeff1,coeff4,coeff5);

  // define angular PDF for signal and generation-level distributions, using the custom class
  RooAbsPdf* PDF_sig_ang_decayRate = new DecayRate(("PDF_sig_ang_decayRate_"+shortString).c_str(),"PDF_sig_ang_decayRate",
						   *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  // define PDF including penalty term
  RooAbsPdf* PDF_sig_ang_decayRate_pen = new DecayRate_Pen(("PDF_sig_ang_decayRate_pen_"+shortString).c_str(),"PDF_sig_ang_decayRate_pen",
							   *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*penTerm);

  // Results' containers
  RooRealVar* fitTime = new RooRealVar("fitTime","fit time",0,"s");
  RooRealVar* co1 = new RooRealVar("co1","Coefficient 1",0);
  RooRealVar* co4 = new RooRealVar("co4","Coefficient 4",0);
  RooRealVar* co5 = new RooRealVar("co5","Coefficient 5",0);
  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  RooArgSet savePars (*co1,*co4,*co5,*fitTime);
  savePars.add(pars);
  RooDataSet* subPosConv = new RooDataSet("subPosConv","subPosConv",savePars);
  RooDataSet* subPosNotc = new RooDataSet("subPosNotc","subPosNotc",savePars);
  RooDataSet* subNegConv = new RooDataSet("subNegConv","subNegConv",savePars);
  RooDataSet* subNegNotc = new RooDataSet("subNegNotc","subNegNotc",savePars);

  // Measure time of the full fit sequence
  TStopwatch subTime;

  // Use adaptive coefficient tuning in case no coefficient is given
  bool isAdaptive = true;
  if ( coeff1>0 || coeff4>0 || coeff5>0 ) isAdaptive = false;


  // Produce and fit subsamples
  for (int iSam=0; iSam<nSamples; ++iSam) {

    // create the subsample
    string selection = Form("rand > %1.6f && rand < %1.6f", iSam*statFactor, (iSam+1)*statFactor);
    data = (RooDataSet*)fullData->reduce( Name(Form("data_%s_%i",shortString.c_str(),iSam)),SelectVars(vars),Cut(selection.c_str()) );
    cout<<"Fitting subsample "<<iSam<<" with "<<data->numEntries()<<" events"<<endl;

    // set penalty power with sample-specific correction
    penTerm->setPower(power/data->numEntries());

    // reset parameters to initial values
    Fl->setVal(0.5);
    P1->setVal(0);
    P2->setVal(0);
    P3->setVal(0);
    P4p->setVal(0);
    P5p->setVal(0);
    P6p->setVal(0);
    P8p->setVal(0);

    // define free NLL
    RooAbsReal* nll = PDF_sig_ang_decayRate->createNLL(*data,
						       RooFit::Extended(kFALSE),
						       RooFit::NumCPU(1)
						       );
    // Penalised NLL will be initialised in each tuning step
    RooAbsReal* nll_pen = 0;
    
    // Start timer
    subTime.Start(true);

    // Prepare free minimizer
    RooMinimizer m (*nll) ;
    m.optimizeConst(kTRUE);
    m.setOffsetting(kTRUE);
    m.setMinimizerType("Minuit2");
    m.setPrintLevel(-1);

    // Fit step-0: quick free fit (strategy=0)
    m.setStrategy(0);
    m.migrad() ;
    m.hesse() ;   

    // Fit step-1: accurate free fit (strategy=2)
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;   

    // Save free-fit result
    RooFitResult* fitResult = m.save(Form("fitResult_%s_%i",shortString.c_str(),iSam),
				     Form("Free fit result to GEN subsample %s %i",shortString.c_str(),iSam));

    RooFitResult* penFitResult = 0;
    bool usedPenalty = false;
    bool isPhysical = false;

    // use penalised fit in case the free fit is not good
    if ( boundary->getValV() > 0 ||
	 fitResult->status() != 0 ||
	 fitResult->covQual() != 3 ) {
    
      usedPenalty = true;

      // Fit step-2: accurate penalalised fit

      // if penalty coefficients are given, a single fit with stable penalty is performed
      // otherwise, the adaptive approach with automatic coefficient tuning is run
      if ( !isAdaptive ) {

	nll_pen = PDF_sig_ang_decayRate_pen->createNLL(*data,
						       RooFit::Extended(kFALSE),
						       RooFit::NumCPU(1)
						       );
	// Prepare penalised minimizer
	RooMinimizer mPen (*nll_pen) ;
	mPen.optimizeConst(kTRUE);
	mPen.setOffsetting(kTRUE);
	mPen.setMinimizerType("Minuit2");
	mPen.setPrintLevel(-1);
	mPen.setStrategy(2);

	// Run fit
	mPen.migrad() ;
	mPen.hesse() ;   

	// Save result
	penFitResult = mPen.save(Form("penFitResult_%s_%i",shortString.c_str(),iSam),
				 Form("Penalised fit result to GEN subsample %s %i",shortString.c_str(),iSam));

	if ( boundary->getValV() <= 0 &&
	     penFitResult->status() == 0 &&
	     penFitResult->covQual() == 3 ) isPhysical = true;

      } else {

	int totCoeff, totCoeffPro, phase1, phase2;
	double sinPh1;

	// Each set of three coefficients is identified by spherical coordinatesof the 3D space

	// Radius
	for (totCoeff=0; coeff1<maxCoeff; ++totCoeff) {

	  // Polar angle = pi/2 * phase1 / totCoeff
	  for (phase1=0; phase1<=totCoeff; ++phase1) {

	    // coeff1 is represented by the z axis
	    coeff1 = fac1 * pow(base1,totCoeff*cos(0.5*TMath::Pi()*phase1/(totCoeff>0?totCoeff:1)));
	    penTerm->setCoefficient(1,coeff1);

	    // Azimuthal angle = pi/2 * phase2 / ( totCoeff * sin(theta) )
	    sinPh1 = sin(0.5*TMath::Pi()*phase1/(totCoeff>0?totCoeff:1));
	    totCoeffPro = (int)(totCoeff*sinPh1);
	    for (phase2=0; phase2<=totCoeffPro; ++phase2) {
	    
	      // coeff4 is represented by the y axis
	      coeff4 = fac4 * pow(base4,totCoeff*sinPh1*sin(0.5*TMath::Pi()*phase2/(totCoeffPro>0?totCoeffPro:1)));
	      penTerm->setCoefficient(4,coeff4);

	      // coeff5 is represented by the x axis
	      coeff5 = fac5 * pow(base5,totCoeff*sinPh1*cos(0.5*TMath::Pi()*phase2/(totCoeffPro>0?totCoeffPro:1)));
	      penTerm->setCoefficient(5,coeff5);

	      nll_pen = PDF_sig_ang_decayRate_pen->createNLL(*data,
							     RooFit::Extended(kFALSE),
							     RooFit::NumCPU(1)
							     );
	      // Prepare penalised minimizer
	      RooMinimizer mPen (*nll_pen) ;
	      mPen.optimizeConst(kTRUE);
	      mPen.setOffsetting(kTRUE);
	      mPen.setMinimizerType("Minuit2");
	      mPen.setPrintLevel(-1);
	      mPen.setStrategy(2);
	    
	      // Run fit
	      mPen.migrad() ;
	      mPen.hesse() ;   
	    
	      // Save result
	      penFitResult = mPen.save(Form("penFitResult_%s_%i",shortString.c_str(),iSam),
				       Form("Penalised fit result to GEN subsample %s %i",shortString.c_str(),iSam));

	      // If fit is good stop tuning
	      if ( boundary->getValV() <= 0 &&
		   penFitResult->status() == 0 &&
		   penFitResult->covQual() == 3 ) {
		isPhysical = true;
		break;
	      }

	    }
	    if (isPhysical) break;

	  }
	  if (isPhysical) break;

	}

      }

    }

    // Stop timer and save output
    subTime.Stop();
    fitTime->setVal(subTime.CpuTime());
    cout<<"Total fit time: "<<subTime.RealTime()<<" s (real), "<<subTime.CpuTime()<<" s (CPU)"<<endl;
  
    // Print results

    co1->setVal(0);
    co4->setVal(0);
    co5->setVal(0);

    double boundCheck = boundary->getValV();

    bool convCheck = false;
    if (usedPenalty) {

      if (penFitResult->status()==0 && penFitResult->covQual()==3) convCheck = true;

      if (isPhysical) {
	cout<<"Physical result with coeff: "<<coeff1<<" "<<coeff4<<" "<<coeff5<<endl<<endl;
	co1->setVal(coeff1);
	co4->setVal(coeff4);
	co5->setVal(coeff5);
      } else cout<<"No physical result"<<endl<<endl;

    } else {
      if (fitResult->status()==0 && fitResult->covQual()==3) convCheck = true;
      cout<<"Physical result with free fit"<<endl<<endl;
    }

    // Fill results' dataset
    if (boundCheck>0) {
      if (convCheck) subNegConv->add(savePars);
      else subNegNotc->add(savePars);
    } else {
      if (convCheck) subPosConv->add(savePars);
      else subPosNotc->add(savePars);
    }

  }

  RooCategory resStatus ("resStatus","Status of the fit result");
  resStatus.defineType("convergent-positive",0);
  resStatus.defineType("convergent-negative",1);
  resStatus.defineType("notconvergent-positive",2);
  resStatus.defineType("notconvergent-negative",3);
  RooDataSet* subResults = new RooDataSet("subResults","Results of GEN sub-sample fitting",
					  savePars,Index(resStatus),
					  Import("convergent-positive",*subPosConv),
					  Import("convergent-negative",*subNegConv),
					  Import("notconvergent-positive",*subPosNotc),
					  Import("notconvergent-negative",*subNegNotc));

  if (save) {

    string foutName = "simFitResults/fitResult_genMC_penalty.root";
    if (!isAdaptive) foutName = "simFitResults/fitResult_genMC_staticPenalty.root";
    TFile* fout = new TFile(foutName.c_str(),"UPDATE");

    RooWorkspace* wksp = new RooWorkspace(Form("wsMulti_%s_s%i_n%i_pow%.1f",shortString.c_str(),stat,nSamples,power),"Workspace with collection of fit result on GEN subsamples");
    wksp->import(*subResults);

    fout->cd();
    wksp->Write();
    fout->Close();
    
  }

  double time90quant = 0;
  double quant = 0;
  double totEntries = subResults->sumEntries();
  // for (time90quant = 0; quant<0.9; time90quant += 0.1)
  while (quant<0.9) {
    time90quant += 0.1;
    quant = subResults->sumEntries(Form("fitTime<%.2f",time90quant))/totEntries;
  }
  // for (timeMax = time90quant; quant<1; timeMax += 0.1)
  double timeMax = time90quant;
  while (quant<1) {
    timeMax += 0.1;
    quant = subResults->sumEntries(Form("fitTime<%.2f",timeMax))/totEntries;
  }
  cout<<"Average fit time: "<<subResults->mean(*fitTime)<<" sec (90% quantile: "<<time90quant<<" sec, max: "<<timeMax<<" sec)"<<endl;


  cout<<"Fitted subsamples: "<<subResults->numEntries()<<" of which good: "<<subPosConv->numEntries()<<endl;
  cout<<"Bad fits: "<<subNegConv->numEntries()<<" converging outside physical region, "
      <<subPosNotc->numEntries()+subNegNotc->numEntries()<<" not converged ("<<subPosNotc->numEntries()<<" in ph region)"<<endl;

  if (!plot) return;

  // Plotting the fit-time distribution
  
  cTime = new TCanvas (("cTime_"+shortString).c_str(),("cTime_"+shortString).c_str(),1800,1800);
  // cTime->Divide(2,1);
 
  // RooBinning timeBins (0.001);

  // RooPlot* fP5p  = P5p    ->frame(Name("fP5p") ,Title("P5p distribution"     ));
  RooPlot* fTime = fitTime->frame(Name("fTime"),Title("Fit-time distribution"),Range(0,1.05*timeMax));
  // RooPlot* fTime = fitTime->frame(Name("fTime"),Title("Fit-time distribution"),AutoRange(*subResults,0.1));

  // subResults->plotOn(fP5p ,Binning(50));
  subResults->plotOn(fTime,Binning(50,0,1.05*timeMax));

  cTime->cd();
  fTime->Draw();

  // cTime->cd(2);
  // fP5p->Draw();

  string penString = "";
  if (isAdaptive) penString = Form("_adPen-%.1f",power);
  else penString = Form("_stPen-%.1f-%.0f-%.0f-%.0f",power,coeff1,coeff4,coeff5);
  
  cTime->Update();
  cTime->SaveAs( ("plotSimFit_d/gen_fitTime_" + shortString + penString + Form("_s%i_n%i.pdf",stat,nSamples)).c_str() );

}

void fit_genMCBin1(int q2Bin, int parity, bool plot, bool save, int nSamples, double power, double coeff1, double coeff4, double coeff5)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      fit_genMCBin(q2Bin, parity, plot, save, nSamples, power, coeff1, coeff4, coeff5);
  else
    fit_genMCBin(q2Bin, parity, plot, save, nSamples, power, coeff1, coeff4, coeff5);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin    = -1;
  int parity   = -1; 
  int nSamples = 100;

  if ( argc >= 2 ) q2Bin    = atoi(argv[1]);
  if ( argc >= 3 ) parity   = atoi(argv[2]);
  if ( argc >= 4 ) nSamples = atoi(argv[3]);

  double power  = 1.0;
  double coeff1 = 0.0;
  double coeff4 = 0.0;
  double coeff5 = 0.0;

  if ( argc >= 5 ) power  = atof(argv[4]);
  if ( argc >= 6 ) coeff1 = atof(argv[5]);
  if ( argc >= 7 ) coeff4 = atof(argv[6]);
  if ( argc >= 8 ) coeff5 = atof(argv[7]);

  bool plot = true;
  bool save = true;

  if ( argc >= 9  && atoi(argv[8]) == 0 ) plot = false;
  if ( argc >= 10 && atoi(argv[9]) == 0 ) save = false;

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout<<"Running all the q2 bins"<<endl;
  if ( parity==-1 )  cout<<"Running both the parity datasets"<<endl;

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      fit_genMCBin1(q2Bin, parity, plot, save, nSamples, power, coeff1, coeff4, coeff5);
  else
    fit_genMCBin1(q2Bin, parity, plot, save, nSamples, power, coeff1, coeff4, coeff5);

  return 0;

}
