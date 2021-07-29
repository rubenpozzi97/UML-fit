#include <RooAddPdf.h>
#include <RooCBShape.h>
#include <RooGaussian.h>
#include "RooDoubleCBFast.h"
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

RooDoubleCBFast* createRTMassShape(int q2Bin,
                                   RooAbsReal* x,
                                   RooAbsReal* mean_rt,
                                   RooAbsReal* sigma_rt,
                                   RooAbsReal* alpha_rt1,
                                   RooAbsReal* alpha_rt2,
                                   RooAbsReal* n_rt1,
                                   RooAbsReal* n_rt2,
                                   TString input_file,
                                   int year,
                                   bool constrainVars,
                                   RooArgSet &c_vars,
                                   RooArgSet &c_pdfs
                                   ){

    RooDoubleCBFast* dcb_rt = new RooDoubleCBFast (Form("dcb_rt_%i", year),
                                                   "dcb_rt",
                                                   *x,
                                                   *mean_rt, *sigma_rt, *alpha_rt1, *n_rt1, *alpha_rt2, *n_rt2
                                                   );

    if (constrainVars){
        constrainVar3(input_file, Form("#sigma_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{RT2}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{RT1}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{RT2}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs);
    }

    return dcb_rt;
}

RooAddPdf* createRTMassShape2(int q2Bin,
                              RooAbsReal* x,
                              RooAbsReal* mean_rt,
                              RooAbsReal* sigma_rt,
                              RooAbsReal* sigma_rt2,
                              RooAbsReal* alpha_rt1,
                              RooAbsReal* alpha_rt2,
                              RooAbsReal* n_rt1,
                              RooAbsReal* n_rt2,
                              RooAbsReal* f1rt,
                              TString input_file,
                              int year,
                              bool constrainVars,
                              RooArgSet &c_vars,
                              RooArgSet &c_pdfs
                              ){

    RooCBShape* cbshape_rt1 = new RooCBShape (Form("cbshape_rt1_%i", year),
                                              Form("cbshape_rt1_%i", year),
                                              *x,
                                              *mean_rt, *sigma_rt , *alpha_rt1, *n_rt1
                                              );
    RooCBShape* cbshape_rt2 = new RooCBShape (Form("cbshape_rt2_%i", year),
                                              Form("cbshape_rt2_%i", year),
                                              *x,
                                              *mean_rt, *sigma_rt2, *alpha_rt2, *n_rt2
                                              );

    RooAddPdf* dcb_rt = new RooAddPdf (Form("dcb_rt_%i", year) ,
                                       Form("dcb_rt_%i", year) ,
                                       RooArgList(*cbshape_rt1,*cbshape_rt2),
                                       RooArgList(*f1rt));

    if (constrainVars){
        constrainVar3(input_file, Form("#sigma_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{RT2}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{RT1}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{RT2}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#sigma_{RT2}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("f^{RT%i}"         ,year) , year, q2Bin, true, c_vars, c_pdfs);
    }

    return dcb_rt;
}

RooAddPdf* createRTMassShape3(int q2Bin,
                              RooAbsReal* x,
                              RooAbsReal* mean_rt,
                              RooAbsReal* sigma_rt,
                              RooAbsReal* sigma_rt2,
                              RooAbsReal* alpha_rt1,
                              RooAbsReal* n_rt1,
                              RooAbsReal* f1rt,
                              TString input_file,
                              int year,
                              bool constrainVars,
                              RooArgSet &c_vars,
                              RooArgSet &c_pdfs
                              ){

    RooCBShape* cbshape_rt = new RooCBShape (Form("cbshape_rt_%i", year),
                                              Form("cbshape_rt_%i", year),
                                              *x,
                                              *mean_rt, *sigma_rt , *alpha_rt1, *n_rt1
                                              );
    RooGaussian* gaussian_rt = new RooGaussian (Form("gaussian_rt_%i", year),
                                                Form("gaussian_rt_%i", year),
                                                *x,
                                                *mean_rt, *sigma_rt2
                                                );

    RooAddPdf* dcb_rt = new RooAddPdf (Form("dcb_rt_%i", year) ,
                                       Form("dcb_rt_%i", year) ,
                                       RooArgList(*cbshape_rt,*gaussian_rt),
                                       RooArgList(*f1rt));

    if (constrainVars){
        constrainVar3(input_file, Form("#sigma_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{RT1}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#sigma_{RT2}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("f^{RT%i}"         ,year) , year, q2Bin, true, c_vars, c_pdfs);
    }

    return dcb_rt;
}

RooAddPdf* createRTMassShape4(int q2Bin,
                              RooAbsReal* x,
                              RooAbsReal* mean_rt,
                              RooAbsReal* sigma_rt,
                              RooAbsReal* sigma_rt2,
                              RooAbsReal* f1rt,
                              TString input_file,
                              int year,
                              bool constrainVars,
                              RooArgSet &c_vars,
                              RooArgSet &c_pdfs
                              ){

    RooGaussian* gaussian_rt1 = new RooGaussian (Form("gaussian_rt1_%i", year),
                                                 Form("gaussian_rt1_%i", year),
                                                 *x,
                                                 *mean_rt, *sigma_rt
                                                 );
    RooGaussian* gaussian_rt2 = new RooGaussian (Form("gaussian_rt2_%i", year),
                                                 Form("gaussian_rt2_%i", year),
                                                 *x,
                                                 *mean_rt, *sigma_rt2
                                                 );

    RooAddPdf* dcb_rt = new RooAddPdf (Form("dcb_rt_%i", year) ,
                                       Form("dcb_rt_%i", year) ,
                                       RooArgList(*gaussian_rt1,*gaussian_rt2),
                                       RooArgList(*f1rt));

    if (constrainVars){
        constrainVar3(input_file, Form("#sigma_{RT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#sigma_{RT2}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("f^{RT%i}"         ,year) , year, q2Bin, true, c_vars, c_pdfs);
    }

    return dcb_rt;
}
