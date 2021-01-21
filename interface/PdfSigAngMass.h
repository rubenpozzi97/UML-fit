/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef PDFSIGANGMASS
#define PDFSIGANGMASS

#include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooCategoryProxy.h"
#include "RooSetProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooFit.h"
#include "Riostream.h"
#include "RooObjCacheManager.h"
 
class PdfSigAngMass : public RooAbsPdf {
 protected:

  RooRealProxy ctK ;
  RooRealProxy ctL ;
  RooRealProxy phi ;
  RooRealProxy m ;
  RooRealProxy Fl ;
  RooRealProxy P1 ;
  RooRealProxy P2 ;
  RooRealProxy P3 ;
  RooRealProxy P4p ;
  RooRealProxy P5p ;
  RooRealProxy P6p ;
  RooRealProxy P8p ;

  RooRealProxy mean_rt     ;
  RooRealProxy sigma_rt1   ;
  RooRealProxy sigma_rt2   ;
  RooRealProxy alpha_rt1   ;
  RooRealProxy alpha_rt2   ;
  RooRealProxy n_rt1       ;
  RooRealProxy n_rt2       ;
  RooRealProxy f1rt        ;

  RooRealProxy mean_wt     ;
  RooRealProxy sigma_wt1   ;
  RooRealProxy alpha_wt1   ;
  RooRealProxy alpha_wt2   ;
  RooRealProxy n_wt1       ;
  RooRealProxy n_wt2       ;
  
  RooRealProxy mFrac ;

  RooRealProxy PenTerm;
  RooRealProxy rtAngTerm;
  RooRealProxy wtAngTerm;
  RooRealProxy rtMassTerm;
  RooRealProxy wtMassTerm;
  
  bool isPenalised;
  
  const RooAbsReal* penTermVal() const {
    // Return pointer to penalty term function
    return (RooAbsReal*) PenTerm.absArg() ;
  }

  Double_t evaluate() const ;

 public:
  PdfSigAngMass() {} ; 
  // 1 sigma RT, penalty
  PdfSigAngMass(const char *name, const char *title,
	    RooAbsReal& _ctK,
	    RooAbsReal& _ctL,
	    RooAbsReal& _phi,
	    RooAbsReal& _m,
	    RooAbsReal& _Fl,
	    RooAbsReal& _P1,
	    RooAbsReal& _P2,
	    RooAbsReal& _P3,
	    RooAbsReal& _P4p,
	    RooAbsReal& _P5p,
	    RooAbsReal& _P6p,
	    RooAbsReal& _P8p,
	    RooAbsReal& _mean_rt  ,
	    RooAbsReal& _sigma_rt1,
	    RooAbsReal& _alpha_rt1,
	    RooAbsReal& _alpha_rt2,
	    RooAbsReal& _n_rt1    ,
	    RooAbsReal& _n_rt2    ,
	    RooAbsReal& _mean_wt  ,
	    RooAbsReal& _sigma_wt1,
	    RooAbsReal& _alpha_wt1,
	    RooAbsReal& _alpha_wt2,
	    RooAbsReal& _n_wt1    ,
	    RooAbsReal& _n_wt2    ,
	    RooAbsReal& _mFrac,
	    RooAbsReal& _PenTerm,
	    RooAbsReal& _rtAngTerm,
	    RooAbsReal& _wtAngTerm,
	    RooAbsReal& _rtMassTerm,
	    RooAbsReal& _wtMassTerm
	    );
  // 1 sigma RT, no penalty
  PdfSigAngMass(const char *name, const char *title,
	    RooAbsReal& _ctK,
	    RooAbsReal& _ctL,
	    RooAbsReal& _phi,
	    RooAbsReal& _m,
	    RooAbsReal& _Fl,
	    RooAbsReal& _P1,
	    RooAbsReal& _P2,
	    RooAbsReal& _P3,
	    RooAbsReal& _P4p,
	    RooAbsReal& _P5p,
	    RooAbsReal& _P6p,
	    RooAbsReal& _P8p,
	    RooAbsReal& _mean_rt  ,
	    RooAbsReal& _sigma_rt1,
	    RooAbsReal& _alpha_rt1,
	    RooAbsReal& _alpha_rt2,
	    RooAbsReal& _n_rt1    ,
	    RooAbsReal& _n_rt2    ,
	    RooAbsReal& _mean_wt  ,
	    RooAbsReal& _sigma_wt1,
	    RooAbsReal& _alpha_wt1,
	    RooAbsReal& _alpha_wt2,
	    RooAbsReal& _n_wt1    ,
	    RooAbsReal& _n_wt2    ,
	    RooAbsReal& _mFrac,
	    RooAbsReal& _rtAngTerm,
	    RooAbsReal& _wtAngTerm,
	    RooAbsReal& _rtMassTerm,
	    RooAbsReal& _wtMassTerm
	    );
  // 2 sigma RT,  penalty
  PdfSigAngMass(const char *name, const char *title,
	    RooAbsReal& _ctK,
	    RooAbsReal& _ctL,
	    RooAbsReal& _phi,
	    RooAbsReal& _m,
	    RooAbsReal& _Fl,
	    RooAbsReal& _P1,
	    RooAbsReal& _P2,
	    RooAbsReal& _P3,
	    RooAbsReal& _P4p,
	    RooAbsReal& _P5p,
	    RooAbsReal& _P6p,
	    RooAbsReal& _P8p,
	    RooAbsReal& _mean_rt  ,
	    RooAbsReal& _sigma_rt1,
	    RooAbsReal& _sigma_rt2,
	    RooAbsReal& _alpha_rt1,
	    RooAbsReal& _alpha_rt2,
	    RooAbsReal& _n_rt1    ,
	    RooAbsReal& _n_rt2    ,
            RooAbsReal& _f1rt     ,
	    RooAbsReal& _mean_wt  ,
	    RooAbsReal& _sigma_wt1,
	    RooAbsReal& _alpha_wt1,
	    RooAbsReal& _alpha_wt2,
	    RooAbsReal& _n_wt1    ,
	    RooAbsReal& _n_wt2    ,
	    RooAbsReal& _mFrac,
	    RooAbsReal& _PenTerm,
	    RooAbsReal& _rtAngTerm,
	    RooAbsReal& _wtAngTerm,
	    RooAbsReal& _rtMassTerm,
	    RooAbsReal& _wtMassTerm
	    );

  // 2 sigma RT, no penalty
  PdfSigAngMass(const char *name, const char *title,
	    RooAbsReal& _ctK,
	    RooAbsReal& _ctL,
	    RooAbsReal& _phi,
	    RooAbsReal& _m,
	    RooAbsReal& _Fl,
	    RooAbsReal& _P1,
	    RooAbsReal& _P2,
	    RooAbsReal& _P3,
	    RooAbsReal& _P4p,
	    RooAbsReal& _P5p,
	    RooAbsReal& _P6p,
	    RooAbsReal& _P8p,
	    RooAbsReal& _mean_rt  ,
	    RooAbsReal& _sigma_rt1,
	    RooAbsReal& _sigma_rt2,
	    RooAbsReal& _alpha_rt1,
	    RooAbsReal& _alpha_rt2,
	    RooAbsReal& _n_rt1    ,
	    RooAbsReal& _n_rt2    ,
	    RooAbsReal& _f1rt     ,  
	    RooAbsReal& _mean_wt  ,
	    RooAbsReal& _sigma_wt1,
	    RooAbsReal& _alpha_wt1,
	    RooAbsReal& _alpha_wt2,
	    RooAbsReal& _n_wt1    ,
	    RooAbsReal& _n_wt2    ,
	    RooAbsReal& _mFrac,
	    RooAbsReal& _rtAngTerm,
	    RooAbsReal& _wtAngTerm,
	    RooAbsReal& _rtMassTerm, 
	    RooAbsReal& _wtMassTerm
	    );

  PdfSigAngMass(const PdfSigAngMass& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new PdfSigAngMass(*this,newname); }
  inline virtual ~PdfSigAngMass() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  ClassDef(PdfSigAngMass,1) // PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events
    };
 
#endif
