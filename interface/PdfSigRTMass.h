#include <RooAddPdf.h>
#include <RooRealSumPdf.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooWorkspace.h>


RooDoubleCBFast* createRTMassShape(RooRealVar* x,
                                   RooRealVar* mean_rt,
                                   RooRealVar* sigma_rt,
                                   RooRealVar* alpha_rt1,
                                   RooRealVar* alpha_rt2,
                                   RooRealVar* n_rt1,
                                   RooRealVar* n_rt2,
                                   RooWorkspace *w,
                                   int year
                                   ){

    RooDoubleCBFast* dcb_rt = new RooDoubleCBFast ( Form("dcb_rt_%i", year)  , 
                                                   "dcb_rt"      , 
                                                   *x, 
                                                   *mean_rt, *sigma_rt, *alpha_rt1, *n_rt1, *alpha_rt2, *n_rt2
                                                   );
    return dcb_rt;                                                   
}




RooRealSumPdf* createRTMassShape( RooRealVar* x,
                                  RooRealVar* mean_rt,
                                  RooRealVar* sigma_rt,
                                  RooRealVar* sigma_rt2,
                                  RooRealVar* alpha_rt1,
                                  RooRealVar* alpha_rt2,
                                  RooRealVar* n_rt1,
                                  RooRealVar* n_rt2,
                                  RooRealVar* f1rt,
                                  RooWorkspace *w,
                                  int year
                                  ){

    RooCBShape* cbshape_rt1 = new RooCBShape (Form("cbshape_rt1_%i", year) , 
                                              Form("cbshape_rt1_%i", year) ,  
                                              *x, 
                                              *mean_rt, *sigma_rt , *alpha_rt1, *n_rt1
                                              );
    RooCBShape* cbshape_rt2 = new RooCBShape (Form("cbshape_rt2_%i", year) , 
                                              Form("cbshape_rt2_%i", year) ,  
                                              *x, 
                                              *mean_rt, *sigma_rt2, *alpha_rt2, *n_rt2
                                              );
    RooRealSumPdf* dcb_rt = new RooRealSumPdf (Form("dcb_rt_%i", year) , 
                                       Form("dcb_rt_%i", year) ,  
                                       RooArgList(*cbshape_rt1,*cbshape_rt2), 
                                       RooArgList(*f1rt));

    return dcb_rt;                                                   
}