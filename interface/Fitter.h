#ifndef FITTER
#define FITTER

#include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TH1D.h>
#include <TString.h>
#include <TStopwatch.h>

#include <RooFitResult.h>
#include <RooDataSet.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooMinimizer.h>
#include <RooArgList.h>
#include <RooRealVar.h>

#include "BoundCheck.h"
#include "BoundDist.h"
#include "Penalty.h"

class Fitter {

 protected:

  TString name;
  TString title;

  RooArgList angPars;

  RooDataSet* combData;
  RooAbsPdf* simPdf;
  RooAbsPdf* simPdf_penalty;

  BoundCheck* boundary;
  BoundDist* bound_dist;
  Penalty* penTerm;

  RooArgSet* constrVars;

  Double_t fac1;
  Double_t fac4;
  Double_t base1;
  Double_t base4;
  Double_t max1;
  Double_t max4;
  Double_t base1_corr;
  Double_t base4_corr;
  Double_t min_base;

  void fillResultContainers(bool fromImprov = false) ;

  Double_t computeBoundaryDistance() ;

 public:

  Double_t maxCoeff;
  Double_t coeff1;
  Double_t coeff4;
  Double_t coeff5;
  Bool_t usedPenalty;

  RooFitResult* result_free;
  RooFitResult* result_penalty;

  RooAbsReal* nll;
  RooAbsReal* nll_penalty;

  Double_t boundDist;
  Double_t boundDistTime;

  // lower threshold to parameters' uncertainties
  // to build the randomisation models (too small leads to many useless points)
  double minParError;

  // MINOS parameters
  double widthScale;

  std::vector<Double_t> vFitResult;
  std::vector<Double_t> vFitErrLow;
  std::vector<Double_t> vFitErrHigh;
  std::vector<Double_t> vImprovResult;

  std::vector<Double_t> vResult; /* Containing always the latest result, either from a fit or fit-improv */
  std::vector<Double_t> vConfInterLow; /* Custom MINOS error */
  std::vector<Double_t> vConfInterHigh; /* Custom MINOS error */

  Fitter() {} ; 
  Fitter(const char *_name, const char *_title,
	 RooArgList _angPars,
	 RooDataSet* _combData,
	 RooAbsPdf* _simPdf,
	 RooAbsPdf* _simPdf_penalty,
	 BoundCheck* _boundary,
	 BoundDist* _bound_dist,
	 Penalty* _penTerm,
	 RooArgSet* _constrVars = nullptr
	 );
  Fitter(const char *_name, const char *_title,
	 RooArgList _angPars,
	 RooDataSet* _combData,
	 RooAbsPdf* _simPdf,
	 RooAbsPdf* _simPdf_penalty,
	 BoundCheck* _boundary,
	 BoundDist* _bound_dist,
	 RooArgSet* _constrVars = nullptr
	 );
  Fitter(const Fitter& other, const char* name=0) ;
  virtual Fitter* clone(const char* newname) const { return new Fitter(*this,newname); }
  inline virtual ~Fitter() { }

  void SetDefConf() ;

  Int_t fit() ;

  Int_t improveAng(int seed = 1, int nGen = 10000) ;

  Int_t MinosAng(int seed = 1, int nGenMINOS = 20000) ;

  RooFitResult* result () { if (usedPenalty) return result_penalty; return result_free; };

  ClassDef(Fitter,1) // Code to run the fit and statistical uncertainty
};
 
#endif
