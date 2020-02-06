#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <TStopwatch.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>

#include "ParBound.h"

using namespace RooFit;
using namespace std;

void integrate_boundary(double power, double shift1, double shift5, double coeff)
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
 
  RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* P1    = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2    = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3    = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* R45   = new RooRealVar("R45","R45",0.01,0,1);
  RooRealVar* R68   = new RooRealVar("R68","R68",0.01,0,1);
  RooRealVar* phi45 = new RooRealVar("phi45","phi45",TMath::Pi(),0,2*TMath::Pi());
  RooRealVar* phi68 = new RooRealVar("phi68","phi68",TMath::Pi(),0,2*TMath::Pi());

  // Define the PDF defining the physical region, used as exexternal contraint in the fit
  RooAbsPdf* PDF_phys_bound = new ParBound("PDF_phys_bound","PDF_phys_bound",*P1,*P2,*P3,*R45,*R68,*phi45,*phi68,power,shift1,shift5,coeff);

  PDF_phys_bound->getIntegratorConfig()->methodND().setLabel("RooMCIntegrator");
  PDF_phys_bound->getIntegratorConfig()->setEpsRel(1e-7);
  PDF_phys_bound->getIntegratorConfig()->setEpsAbs(1e-7);
  int nInt = 2e7;
  int nRefine = 5;
  PDF_phys_bound->getIntegratorConfig()->getConfigSection("RooMCIntegrator").setRealValue("nIntPerDim",nInt);
  PDF_phys_bound->getIntegratorConfig()->getConfigSection("RooMCIntegrator").setRealValue("nRefinePerDim",nInt/nRefine);
  PDF_phys_bound->getIntegratorConfig()->getConfigSection("RooMCIntegrator").setRealValue("nRefineIter",nRefine);
  cout<<"Number of points per dimension: "<<nInt<<endl;
  cout<<"Number of grid-refinement steps: "<<nRefine<<endl;
  TStopwatch t;
  t.Start();
  double intValue=PDF_phys_bound->getNorm(RooArgSet(*P1,*P2,*P3,*R45,*R68,*phi45,*phi68));
  t.Stop();
  t.Print();
  cout<<"source simfit_recoMC_fullAngular.sh "<<power<<" "<<shift1<<" "<<shift5<<" "<<coeff<<" "<<intValue<<endl;
  return;

}

int main(int argc, char** argv)
{

  double power  = 2;
  double shift1 = 1000;
  double shift5 = 100;
  double coeff  = 50;

  if ( argc >= 2 ) power  = atoi(argv[1]);
  if ( argc >= 3 ) shift1 = atof(argv[2]);
  if ( argc >= 4 ) shift5 = atof(argv[3]);
  if ( argc >= 5 ) coeff  = atof(argv[4]);

  integrate_boundary(power,shift1,shift5,coeff);

  return 0;

}
