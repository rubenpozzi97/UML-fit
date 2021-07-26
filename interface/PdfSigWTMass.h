#include "RooDoubleCBFast.h"
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

RooDoubleCBFast* createWTMassShape(int q2Bin,
                                   RooAbsReal* x,
                                   RooAbsReal* mean_wt,
                                   RooAbsReal* sigma_wt,
                                   RooAbsReal* alpha_wt1,
                                   RooAbsReal* alpha_wt2,
                                   RooAbsReal* n_wt1,
                                   RooAbsReal* n_wt2,
                                   TString input_file,
                                   int year,
                                   bool constrainVars, 
                                   RooArgSet &c_vars,
                                   RooArgSet &c_pdfs
                                   ){

    RooDoubleCBFast* dcb_wt = new RooDoubleCBFast (Form("dcb_wt_%i", year), 
                                                   "dcb_wt", 
                                                   *x, 
                                                   *mean_wt, *sigma_wt, *alpha_wt1, *n_wt1, *alpha_wt2, *n_wt2
                                                   );
   
    if (constrainVars){
        constrainVar3(input_file, Form("#sigma_{WT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{WT1}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("#alpha_{WT2}^{%i}",year) , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{WT1}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs);
        constrainVar3(input_file, Form("n_{WT2}^{%i}",year)      , year, q2Bin, true, c_vars, c_pdfs); 
    }
    
    return dcb_wt;                                                   
}




