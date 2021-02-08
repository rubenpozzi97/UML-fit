#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <map>

using namespace std;
using namespace RooFit;

extern std::map<int, float> scale_to_data;

// redo jpsi and psi2s
std::map<int,std::vector<float>> frt_sigmas = {
  {2016, {0.021, 0.015, 0.016, 0.011, 0.009, 0.009, 0.008, 0.011}},
  {2017, {0.018, 0.013, 0.014, 0.010, 0.008, 0.008, 0.007, 0.011}},
  {2018, {0.013, 0.010, 0.011, 0.007, 0.006, 0.006, 0.006, 0.007}},
};

std::map<int,std::vector<float>> fM_sigmas = {
  {2016, {0.024, 0.017, 0.018, 0.013, 0.005, 0.010, 0.006, 0.010}},
  {2017, {0.020, 0.015, 0.016, 0.011, 0.004, 0.009, 0.005, 0.013}},
  {2018, {0.015, 0.011, 0.013, 0.008, 0.002, 0.006, 0.003, 0.008}},
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

void constrainVar2(RooRealVar* var, 
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
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
}


bool retrieveWorkspace(string filename, std::vector<RooWorkspace*> &ws, std::string ws_name){

    TFile* f =  TFile::Open( filename.c_str() ) ;
    if ( !f || !f->IsOpen() ) {
      cout << "File not found: " << filename << endl;
      return false;
    }
    RooWorkspace* open_w = (RooWorkspace*)f->Get(ws_name.c_str());
    if ( !open_w || open_w->IsZombie() ) {
      cout<<"Workspace "<< ws_name <<  "not found in file: " << filename << endl;
      return false;
    }
    ws.push_back( open_w );
    f->Close();
    return true;
}


std::vector<RooDataSet*> createDataset(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, 
                                       int q2Bin, int parity, int year, //std::map<int,float> scale_to_data,
                                       RooArgSet reco_vars, std::string shortString  ){

    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

	RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i", is)).c_str(), 
					     ("data_"+shortString + Form("_subs%i", is)).c_str(), 
					     RooArgSet(reco_vars));

        dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
        dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;

        isample->append(*dataCT);
        isample->append(*dataWT);
        datasample.push_back (isample);
      }
    }
    else{
      RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
					   ("data_"+shortString + "_subs0").c_str(), 
					   RooArgSet(reco_vars));
      dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
      isample->append(*dataCT);
      isample->append(*dataWT);
      datasample.push_back (isample);
    }
    return datasample;
}

