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

#include "DecayRate.h"
#include "DecayRate_Pen.h"
#include "BoundCheck.h"
#include "Penalty.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* cPen;
TCanvas* c [2*nBins];

double genStat  [nBins] = {1700246,3190316,2416720,4413921,1,5184529,1,3094894,1};
double dataStat [nBins] = {994,2288,1979,3403,0,5381,0,3506,0};
double d16Stat  [nBins] = {246,605,519,873,0,1301,0,908,0};

double zoomMax = 20000.0;
double zoomMaxSub = 2.0;

double fac1 = 1.0;
double fac4 = 1.0;
double fac5 = 1.0;
double base1 = 1.2;
double base4 = 1.2;
double base5 = 1.2;
double maxCoeff = 1e5;

void fit_genMCBin(int q2Bin, int parity, bool plot, bool save, int nSample, double power, double coeff1, double coeff4, double coeff5)
{

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

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
  if (nSample==0) data = new RooDataSet(*fullData);
  else if (nSample<=10000) data = (RooDataSet*)fullData->reduce( RooArgSet(vars),
							       Form("rand > %1.6f && rand < %1.6f",
								    (nSample-1)*dataStat[q2Bin]/genStat[q2Bin],
								    nSample    *dataStat[q2Bin]/genStat[q2Bin]));
  else data = (RooDataSet*)fullData->reduce( RooArgSet(vars),
					     Form("rand > %.6f && rand < %.6f",
						  (nSample%10000-1)*d16Stat[q2Bin]/genStat[q2Bin],
						  (nSample%10000)  *d16Stat[q2Bin]/genStat[q2Bin]));

  cout<<"Fitting "<<data->numEntries()<<" events"<<endl;

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
  Penalty* penTerm = new Penalty("penTerm","Penalty term",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,power/data->numEntries(),coeff1,coeff4,coeff5);
  RooFormulaVar* penLog = new RooFormulaVar("penLog","penLog","-1.0 * log(penTerm)",RooArgList(*penTerm));

  // define angular PDF for signal and generation-level distributions, using the custom class
  RooAbsPdf* PDF_sig_ang_decayRate = new DecayRate(("PDF_sig_ang_decayRate_"+shortString).c_str(),"PDF_sig_ang_decayRate",
						   *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  // define PDF including penalty term
  RooAbsPdf* PDF_sig_ang_decayRate_pen = new DecayRate_Pen(("PDF_sig_ang_decayRate_pen_"+shortString).c_str(),"PDF_sig_ang_decayRate_pen",
							   *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*penTerm);

  // Measure time of the full fit sequence
  TStopwatch subTime;
  subTime.Start(true);

  // define free NLL
  RooAbsReal* nll = PDF_sig_ang_decayRate->createNLL(*data,
						     RooFit::Extended(kFALSE),
						     RooFit::NumCPU(1)
						     );
  // Penalised NLL will be initialised in each tuning step
  RooAbsReal* nll_pen = 0;
    
  // Prepare free minimizer
  RooMinimizer m (*nll) ;
  m.optimizeConst(kTRUE);
  m.setOffsetting(kTRUE);
  // m.setVerbose(kTRUE);
  m.setMinimizerType("Minuit2");
  m.setProfile(kTRUE);

  // Fit step-0: quick free fit (strategy=0)
  m.setStrategy(0);
  m.migrad() ;
  m.hesse() ;   

  // Fit step-1: accurate free fit (strategy=2)
  m.setStrategy(2);
  m.migrad() ;
  m.hesse() ;   

  // Save free-fit result
  RooFitResult* fitResult = m.save("fitResult",Form("Free fit result to GEN subsample %s %i",shortString.c_str(),nSample));

  RooFitResult* penFitResult = 0;
  bool usedPenalty = false;
  bool isPhysical = false;
  bool isAdaptive = false;

  // use penalised fit in case the free fit is not good
  if ( boundary->getValV() > 0 ||
       fitResult->status() != 0 ||
       fitResult->covQual() != 3 ) {
    
    usedPenalty = true;

    // Fit step-2: accurate penalalised fit

    // if penalty coefficients are parsed, a single fit with stable penalty is performed
    // otherwise, the adaptive approach with automatic coefficient tuning is run
    if ( coeff1>0 || coeff4>0 || coeff5>0 ) {
      
      nll_pen = PDF_sig_ang_decayRate_pen->createNLL(*data,
							 RooFit::Extended(kFALSE),
							 RooFit::NumCPU(1)
							 );
      // Prepare penalised minimizer
      RooMinimizer mPen (*nll_pen) ;
      mPen.optimizeConst(kTRUE);
      mPen.setOffsetting(kTRUE);
      mPen.setMinimizerType("Minuit2");
      mPen.setProfile(kTRUE);
      mPen.setStrategy(2);

      // Run fit
      mPen.migrad() ;
      mPen.hesse() ;   

      // Save result
      penFitResult = mPen.save("penFitResult",Form("Penalised fit result to GEN subsample %s %i with coefficients %.1f %.0f %.0f %.0f",shortString.c_str(),nSample,power,coeff1,coeff4,coeff5));

      if ( boundary->getValV() <= 0 &&
	   penFitResult->status() == 0 &&
	   penFitResult->covQual() == 3 ) isPhysical = true;

    } else {

      isAdaptive = true;

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
	    penFitResult = mPen.save("penFitResult",Form("Penalised fit result to GEN subsample %s %i with coefficients %.1f %.0f %.0f %.0f",shortString.c_str(),nSample,power,coeff1,coeff4,coeff5));

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

  // Stop timer and print output
  subTime.Stop();
  cout<<"Total fit time: "<<subTime.RealTime()<<" s (real), "<<subTime.CpuTime()<<" s (CPU)"<<endl;
  
  // Print results
  cout<<endl<<"========== FREE FIT ==========="<<endl;
  fitResult->Print("v");
  cout<<endl;
  
  if (usedPenalty) {

    cout<<"========== PENALISED FIT ==========="<<endl;
    penFitResult->Print("v");
    cout<<endl;

    if (isPhysical) cout<<"Physical result with coeff: "<<coeff1<<" "<<coeff4<<" "<<coeff5<<endl<<endl;
    else cout<<"No physical result"<<endl<<endl;

  } else cout<<"Penalty term not needed"<<endl<<endl;

  if (boundary->getValV()>0) cout<<"Result OUTSIDE physical region"<<endl<<endl;
  else cout<<"Result INSIDE physical region"<<endl<<endl;

  if (save) {

    string foutName = "simFitResults/fitResult_genMC_penalty.root";
    if (usedPenalty && !isAdaptive) foutName = "simFitResults/fitResult_genMC_staticPenalty.root";
    TFile* fout = new TFile(foutName.c_str(),"UPDATE");

    RooWorkspace* wksp = new RooWorkspace(Form("ws_%s_s%i_pow%.1f",shortString.c_str(),nSample,power),"Workspace with GEN fit result");
    wksp->import(*nll);
    wksp->import(*fitResult);
    if (usedPenalty) {
      wksp->import(*penLog);
      wksp->import(*nll_pen);
      wksp->import(*penFitResult);
    }

    fout->cd();
    wksp->Write();
    fout->Close();
    
  }

  if (!plot) return;

  // For plotting the effective penalty term is used
  penTerm->setPower(power);
  
  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  if (nSample>0) zoomMax = zoomMaxSub;
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

    nll->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;
    nll->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed),LineWidth(2)) ;

    if (iPar>0) {

      double hMax = frame[iPar]->GetMaximum();
      
      boundary->plotOn(frame[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
      boundary->plotOn(fZoom[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMax,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));

      if (usedPenalty) {

	nll_pen->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll_pen->getVal()+10),LineColor(kBlue),LineWidth(2));
	nll_pen->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll_pen->getVal()+10),LineColor(kBlue),LineWidth(2));

	penLog->plotOn(frame[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));
	penLog->plotOn(fZoom[iPar],PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(penLog->getVal()+10),LineColor(8),LineWidth(2));

      }
      
      frame[iPar]->SetMaximum(hMax);

      fPenTerm[iPar] = par->frame(Name(Form("f3%s",par->GetName())),Title(Form("Penalty term vs %s",par->GetTitle()))) ;
      penTerm->plotOn(fPenTerm[iPar],LineColor(4),LineWidth(2)) ;
      double hMaxP = fPenTerm[iPar]->GetMaximum();
      boundary->plotOn(fPenTerm[iPar],LineColor(13),FillColor(13),FillStyle(3545),Normalization(1.1*hMaxP,RooAbsReal::Raw),DrawOption("LF"),VLines(),LineWidth(2));
      fPenTerm[iPar]->SetMaximum(hMaxP);
      cPen->cd(iPar+1);
      fPenTerm[iPar]->Draw();

    }

    fZoom[iPar]->SetMaximum(zoomMax);

    cnll->cd(iPar+1);
    frame[iPar]->Draw();

    cZoom->cd(iPar+1);
    fZoom[iPar]->Draw();


  }

  string penString = "";
  if (usedPenalty) {
    if (isAdaptive) penString = Form("_adPen-%.1f",power);
    else penString = Form("_stPen-%.1f-%.0f-%.0f-%.0f",power,coeff1,coeff4,coeff5);
  }
  
  cnll->Update();
  cnll->SaveAs( ("plotSimFit_d/genNLL_scan_" + shortString + penString + (nSample>0?Form("_s%i.pdf",nSample):".pdf")).c_str() );

  cZoom->Update();
  cZoom->SaveAs( ("plotSimFit_d/genNLL_scan_" + shortString + penString + (nSample>0?Form("_s%i_zoom.pdf",nSample):"_zoom.pdf")).c_str() );

  cPen->Update();
  cPen->SaveAs( ("plotSimFit_d/genPenTerm_" + shortString + penString + (nSample>0?Form("_s%i.pdf",nSample):".pdf")).c_str() );


  int confIndex = nBins*parity + q2Bin;
  string longString  = "Fit to generation-level distributions";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to GEN-level MC - "+longString).c_str(),2000,700);
  TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);
  RooPlot* xframe = ctK->frame(Title(longString.c_str()));
  RooPlot* yframe = ctL->frame(Title(longString.c_str()));
  RooPlot* zframe = phi->frame(Title(longString.c_str()));
  data->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(100),Name("plData"));
  data->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(100));
  data->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(100));
  PDF_sig_ang_decayRate->plotOn(xframe,LineWidth(1),Name("plPDF"));
  PDF_sig_ang_decayRate->plotOn(yframe,LineWidth(1));
  PDF_sig_ang_decayRate->plotOn(zframe,LineWidth(1));
  xframe->GetYaxis()->SetTitleOffset(1.8);
  yframe->GetYaxis()->SetTitleOffset(1.8);
  zframe->GetYaxis()->SetTitleOffset(1.8);
  xframe->SetMaximum(xframe->GetMaximum()*1.15);
  yframe->SetMaximum(yframe->GetMaximum()*1.15);
  zframe->SetMaximum(zframe->GetMaximum()*1.15);
  xframe->SetMinimum(0);
  yframe->SetMinimum(0);
  zframe->SetMinimum(0);
  leg->SetTextSize(0.03);
  leg->AddEntry(xframe->findObject("plData"),"Generation-level distribution" ,"lep");
  leg->AddEntry(xframe->findObject("plPDF" ),"Angular decay rate" ,"l");

  c[confIndex]->Divide(3,1);
  c[confIndex]->cd(1);
  gPad->SetLeftMargin(0.15);
  xframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(2);
  gPad->SetLeftMargin(0.15);
  yframe->Draw();
  leg->Draw("same");
  c[confIndex]->cd(3);
  gPad->SetLeftMargin(0.15);
  zframe->Draw();
  leg->Draw("same");

  c[confIndex]->SaveAs( ("plotSimFit_d/fitProj_genMC_" + shortString + penString + (nSample>0?Form("_s%i.pdf",nSample):".pdf")).c_str() );

}

void fit_genMCBin1(int q2Bin, int parity, bool plot, bool save, int nSample, double power, double coeff1, double coeff4, double coeff5)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      fit_genMCBin(q2Bin, parity, plot, save, nSample, power, coeff1, coeff4, coeff5);
  else
    fit_genMCBin(q2Bin, parity, plot, save, nSample, power, coeff1, coeff4, coeff5);
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
  int nSample = 0;

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);
  if ( argc >= 4 ) nSample = atoi(argv[3]);

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
      fit_genMCBin1(q2Bin, parity, plot, save, nSample, power, coeff1, coeff4, coeff5);
  else
    fit_genMCBin1(q2Bin, parity, plot, save, nSample, power, coeff1, coeff4, coeff5);

  return 0;

}
