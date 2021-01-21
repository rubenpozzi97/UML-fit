/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef SHAPESIGANG
#define SHAPESIGANG

#include <math.h>
#include "Math/SpecFunc.h"
#include "TMath.h"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooFit.h"
#include "Riostream.h"
#include "RooObjCacheManager.h"
 
class ShapeSigAng : public RooAbsPdf {
 protected:

  RooRealProxy ctK ;
  RooRealProxy ctL ;
  RooRealProxy phi ;
  RooRealProxy Fl ;
  RooRealProxy P1 ;
  RooRealProxy P2 ;
  RooRealProxy P3 ;
  RooRealProxy P4p ;
  RooRealProxy P5p ;
  RooRealProxy P6p ;
  RooRealProxy P8p ;

  RooRealProxy EffC ;
  RooRealProxy EffW ;

  std::vector<double> intCPart;
  std::vector<double> intWPart;

//   RooRealProxy isC;

  bool isC;
  
  const RooAbsReal* effCVal() const { 
    // Return pointer to efficiency function in product
    return (RooAbsReal*) EffC.absArg() ; 
  }

  const RooAbsReal* effWVal() const { 
    // Return pointer to efficiency function in product
    return (RooAbsReal*) EffW.absArg() ; 
  }

  Double_t evaluate() const ;

 public:
  ShapeSigAng() {} ; 
  ShapeSigAng(const char *name, const char *title,
	    RooAbsReal& _ctK,
	    RooAbsReal& _ctL,
	    RooAbsReal& _phi,
	    RooAbsReal& _Fl,
	    RooAbsReal& _P1,
	    RooAbsReal& _P2,
	    RooAbsReal& _P3,
	    RooAbsReal& _P4p,
	    RooAbsReal& _P5p,
	    RooAbsReal& _P6p,
	    RooAbsReal& _P8p,
	    RooAbsReal& _EffC,
	    RooAbsReal& _EffW,
	    std::vector<double> _intCPart,
	    std::vector<double> _intWPart,
	    bool _isC);
  ShapeSigAng(const ShapeSigAng& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new ShapeSigAng(*this,newname); }
  inline virtual ~ShapeSigAng() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  ClassDef(ShapeSigAng,1) // PDF for (angular decay rate x efficiency) of both correctly-tagged and wrongly-tagged events
    };
 
#endif
