#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <map>
#include <TStyle.h>
#include <RooMCStudy.h>
#include <TF1.h>
#include <RooAbsPdf.h>
#include <RooCategory.h> 
#include <RooArgSet.h>

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
    RooGaussian* gauss_constr = new RooGaussian(Form("c_%s_%i", inVarName.c_str(), year), 
                                                Form("c_%s_%i", inVarName.c_str(), year), 
                                                *var,  
                                                RooConst(w->var(inVarName.c_str())->getVal()), 
                                                RooConst(w->var(inVarName.c_str())->getError())
                                                ); 
    //gauss_constr ->Print();
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
    RooGaussian* gauss_constr = new RooGaussian(Form("c_%s_%i", inVarName.c_str(), year) , 
        					Form("c_%s_%i", inVarName.c_str(), year),
                                                *var,
                                                RooConst(w->var(inVarName.c_str())->getVal()), 
                                                RooConst(w->var(inVarName.c_str())->getError())
                                                ); 
    if (addToList){
      c_vars.add(*var);    
      c_pdfs.add(*gauss_constr);
    }
}


bool retrieveWorkspace(string filename, std::vector<RooWorkspace*> &ws, std::string ws_name, std::vector<RooWorkspace*> &ws1, std::string ws_name_odd, int parity){

    TFile* f =  TFile::Open( filename.c_str() ) ;
    if ( !f || !f->IsOpen() ) {
      cout << "File not found: " << filename << endl;
      return false;
    }
    RooWorkspace* open_w = (RooWorkspace*)f->Get(ws_name.c_str());
    RooWorkspace* open_w_odd = (RooWorkspace*)f->Get(ws_name_odd.c_str()); 

   if ( !open_w || open_w->IsZombie() ) {
      cout<<"Workspace "<< ws_name <<  "not found in file: " << filename << endl;
      return false;
    }
    if(parity < 2){ws.push_back(open_w);}
    else{
      ws.push_back(open_w); // even
      ws1.push_back(open_w_odd); // odd
    }
    f->Close();
    return true;
}


std::vector<RooDataSet*> createDataset(int nSample, uint firstSample, uint lastSample, RooWorkspace *ws, RooWorkspace *ws1, 
                                       int q2Bin, int parity, int year, //std::map<int,float> scale_to_data,
                                       RooArgSet reco_vars, std::string shortString, int comp){

    RooDataSet* dataCT, *dataWT;
    RooDataSet* dataCT_even, *dataWT_even;
    RooDataSet* dataCT_odd, *dataWT_odd;
    std::vector<RooDataSet*> datasample;

    if (nSample>0){  
      for (uint is = firstSample; is <= lastSample; is++) {

	RooDataSet* isample = new RooDataSet(("data_"+shortString + Form("_subs%i", is)).c_str(), 
					     ("data_"+shortString + Form("_subs%i", is)).c_str(), 
					     RooArgSet(reco_vars));

        if(parity < 2){
          dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
            ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
          dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
            ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
        }
        else{
          dataCT_even = (RooDataSet*)ws->data(Form("data_ctRECO_ev_b%i",q2Bin))
            ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
          dataCT_odd = (RooDataSet*)ws1->data(Form("data_ctRECO_od_b%i",q2Bin))
            ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
          dataWT_even = (RooDataSet*)ws->data(Form("data_wtRECO_ev_b%i",q2Bin))
            ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
          dataWT_odd = (RooDataSet*)ws1->data(Form("data_wtRECO_od_b%i",q2Bin))
            ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[year], (is+1)*scale_to_data[year] )) ;
        }

        if(comp == 0){
          if(parity < 2){ 
            isample->append(*dataCT);
          }
          else{
            isample->append(*dataCT_even);
            isample->append(*dataCT_odd);
          }
        }
        else if(comp == 1){
          if(parity < 2){
            isample->append(*dataWT);
          }
          else{
            isample->append(*dataWT_even);
            isample->append(*dataWT_odd);
          }
        } 
        else{
          if(parity < 2){
            isample->append(*dataCT);
            isample->append(*dataWT);
          }
          else{
            isample->append(*dataCT_even);
            isample->append(*dataCT_odd);
            isample->append(*dataWT_even);
            isample->append(*dataWT_odd);
          }
        }
        datasample.push_back (isample);
      }
    }
    else{
      RooDataSet* isample = new RooDataSet(("data_"+shortString + "_subs0").c_str(), 
					   ("data_"+shortString + "_subs0").c_str(), 
					   RooArgSet(reco_vars));

      if(parity < 2){
        dataCT = (RooDataSet*)ws->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
        dataWT = (RooDataSet*)ws->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
      }
      else{
        dataCT_even = (RooDataSet*)ws->data(Form(("data_ctRECO_ev_b%i"),q2Bin));
        dataCT_odd = (RooDataSet*)ws1->data(Form(("data_ctRECO_od_b%i"),q2Bin));
        dataWT_even = (RooDataSet*)ws->data(Form(("data_wtRECO_ev_b%i"),q2Bin));
        dataWT_odd = (RooDataSet*)ws1->data(Form(("data_wtRECO_od_b%i"),q2Bin));
      }

      if(comp == 0){
        if(parity < 2){isample->append(*dataCT);}
        else{
          isample->append(*dataCT_even);
          isample->append(*dataCT_odd);
        }
      }
      else if(comp == 1){
        if(parity < 2){isample->append(*dataWT);}
        else{
          isample->append(*dataWT_even);
          isample->append(*dataWT_odd);
        }
      }
      else{
        if(parity < 2){
          isample->append(*dataCT);
          isample->append(*dataWT);
        }
        else{
          isample->append(*dataCT_even);
          isample->append(*dataCT_odd);
          isample->append(*dataWT_even);
          isample->append(*dataWT_odd);
        }
      }
      datasample.push_back (isample); 
    }
 
    return datasample;
}

void validate_fit(RooWorkspace* w, RooCategory sample, int year, int q2Bin, int parity){

  RooRealVar mass = *(w->var("mass"));

  RooAbsPdf* model = w->pdf("simPdf");

  vector<RooRealVar> params;
  params.push_back(*(w->var( Form("sig_yield^{%i}",year) )));
  string var_name = (w->var( Form("sig_yield^{%i}",year)) )->GetName();

  int params_size = params.size();
  cout << "params size = " << params_size << endl;
  cout << "variable value = " << params[0].getVal() << endl;

  cout << "Starting toy MC study" << endl;
  RooMCStudy* mcstudy = new RooMCStudy(*model, RooArgSet(mass,sample), Binned(kTRUE), Silence(kTRUE), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));
  mcstudy->generateAndFit(5000);
  cout << "Ending toy MC study" << endl;

  vector<RooPlot*> framesPull, framesParam;

  for(int i = 0; i < params_size; ++i){
    framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(200)));//FrameRange(-5,5)
    framesPull[i]->SetTitle("");
    framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
    framesParam[i]->SetTitle("");
  }

  vector<TGraph*> h;

  for(int i = 0; i < params_size; ++i){
    h.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
  }

  gStyle->SetOptFit(0111);

  TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

  gPad->SetLeftMargin(0.15);

  for(int i = 0; i < params_size; ++i){
    h[i]->SetTitle("");
    h[i]->Draw();
    c_pull->Update();
    h[i]->Fit("gaus","","",-20,20);
    h[i]->GetFunction("gaus")->SetLineColor(4);
    h[i]->GetFunction("gaus")->SetLineWidth(5);
    h[i]->SetTitle( ("Pull - " + var_name + Form(" b%ip%i year %i",q2Bin, parity, year) ).c_str() );
    h[i]->GetXaxis()->SetTitle("Pull");
    h[i]->GetYaxis()->SetTitle("Toy MCs");
    h[i]->GetXaxis()->SetRangeUser(-10,10);
    h[i]->Draw("same");
  }

  TCanvas* c_params = new TCanvas("params", "params", 900, 800);

  for(int i = 0; i < params_size; ++i){
    c_params->cd();
    framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
    framesParam.at(i)->Draw();
  }

  c_pull->SaveAs( ("plotSimMassFit_pulls/pulls_poisson_" + var_name + Form("b%ip%i_%i",q2Bin,parity,year) + ".gif").c_str() );
  c_params->SaveAs( ("plotSimMassFit_pulls/pulls_params_poisson_" + var_name + Form("b%ip%i_%i",q2Bin,parity,year) + ".gif").c_str() );

}


