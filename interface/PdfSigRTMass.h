#include <RooAddPdf.h>
#include <RooRealSumPdf.h>
#include <RooCBShape.h>
#include "RooDoubleCBFast.h"
#include <RooRealVar.h>
#include <RooWorkspace.h>

RooDoubleCBFast* createRTMassShape(int q2Bin,
                                   RooRealVar* x,
                                   RooRealVar* mean_rt,
                                   RooRealVar* sigma_rt,
                                   RooRealVar* alpha_rt1,
                                   RooRealVar* alpha_rt2,
                                   RooRealVar* n_rt1,
                                   RooRealVar* n_rt2,
                                   RooWorkspace *w,
                                   int year,
                                   bool constrainVars, 
                                   RooArgSet &c_vars,
                                   RooArgSet &c_pdfs
                                   ){

    RooDoubleCBFast* dcb_rt = new RooDoubleCBFast ( Form("dcb_rt_%i", year)  , 
                                                   "dcb_rt"      , 
                                                   *x, 
                                                   *mean_rt, *sigma_rt, *alpha_rt1, *n_rt1, *alpha_rt2, *n_rt2
                                                   );

    if (constrainVars){
        constrainVar2(mean_rt  , Form("mean_{RT}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(sigma_rt , Form("#sigma_{RT1}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(alpha_rt1, Form("#alpha_{RT1}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(alpha_rt2, Form("#alpha_{RT2}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(n_rt1    , Form("n_{RT1}^{%i}",q2Bin)      , w, year, true, c_vars, c_pdfs);
        constrainVar2(n_rt2    , Form("n_{RT2}^{%i}",q2Bin)      , w, year, true, c_vars, c_pdfs);
    }

    return dcb_rt;                                                   
}




RooRealSumPdf* createRTMassShape( int q2Bin,
                                  RooRealVar* x,
                                  RooRealVar* mean_rt,
                                  RooRealVar* sigma_rt,
                                  RooRealVar* sigma_rt2,
                                  RooRealVar* alpha_rt1,
                                  RooRealVar* alpha_rt2,
                                  RooRealVar* n_rt1,
                                  RooRealVar* n_rt2,
                                  RooRealVar* f1rt,
                                  RooWorkspace *w,
                                  int year,
                                  bool constrainVars, 
                                  RooArgSet &c_vars,
                                  RooArgSet &c_pdfs
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

    if (constrainVars){
        constrainVar2(mean_rt ,  Form("mean_{RT}^{%i}",q2Bin)   , w, year, true, c_vars, c_pdfs);
        constrainVar2(sigma_rt , Form("#sigma_{RT1}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(alpha_rt1, Form("#alpha_{RT1}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(alpha_rt2, Form("#alpha_{RT2}^{%i}",q2Bin) , w, year, true, c_vars, c_pdfs);
        constrainVar2(n_rt1    , Form("n_{RT1}^{%i}",q2Bin)      , w, year, true, c_vars, c_pdfs);
        constrainVar2(n_rt2    , Form("n_{RT2}^{%i}",q2Bin)      , w, year, true, c_vars, c_pdfs);
        constrainVar2(sigma_rt2 , Form("#sigma_{RT2}^{%i}",q2Bin), w, year, true, c_vars, c_pdfs);
        constrainVar2(f1rt      , Form("f^{RT%i}"         ,q2Bin), w, year, true, c_vars, c_pdfs);
    }

    return dcb_rt;                                                   
}