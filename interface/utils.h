#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

using namespace std;
using namespace RooFit;

std::map<int,std::vector<float>> frt_sigmas = {
  {2016, {0.0040, 0.0046, 0.0043, 0.0040, 0.0038, 0.0037, 0.0037, 0.0040}},
  {2017, {0.0039, 0.0039, 0.0042, 0.0035, 0.0035, 0.0038, 0.0042, 0.0042}},
  {2018, {0.0025, 0.0028, 0.0028, 0.0025, 0.0025, 0.0026, 0.0026, 0.0026}},
};


RooPlot* prepareFrame(RooPlot* frame){
    frame->GetYaxis()->SetTitleOffset(1.8);
    frame->SetMaximum(frame->GetMaximum()*1.15);
    frame->SetMinimum(0);
    return frame;
}


RooGaussian* constrainVar(RooRealVar* var, 
                          string inVarName,  
                          RooWorkspace *w, 
                          int year,
                          bool addToList,
                          RooArgSet &c_vars,
                          RooArgSet &c_pdfs
                          ){
    RooGaussian* gauss_constr = new RooGaussian(  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  Form("c_%s_%i", inVarName.c_str(), year) , 
                                                  *var,  
                                                  RooConst( w->var(inVarName.c_str())->getVal()  ), 
                                                  RooConst( w->var(inVarName.c_str())->getError())
                                                 ); 
//     gauss_constr ->Print();
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
    return gauss_constr;                 
}


void retrieveWorkspace(string filename, std::vector<RooWorkspace*> &ws, std::string ws_name){

    TFile* f =  TFile::Open( filename.c_str() ) ;
    if ( !f || !f->IsOpen() ) {
      cout << "File not found: " << filename << endl;
      return;
    }
    RooWorkspace* open_w = (RooWorkspace*)f->Get(ws_name.c_str());
    if ( !open_w || open_w->IsZombie() ) {
      cout<<"Workspace "<< ws_name <<  "not found in file: " << filename << endl;
      return;
    }
    ws.push_back( open_w );
    f->Close();
}



