#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooAbsReal.h>
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

void constrainVar2(RooAbsReal* var, 
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

void constrainVar3(TString input_file,
		   int year,
		   int q2Bin,
                   bool addToList,
                   RooArgSet &c_vars,
		   RooArgSet &c_pdfs
                   ){

  TFile* f = TFile::Open(input_file);
  RooFitResult* fitresult = (RooFitResult*)f->Get(Form("simFitResult_b%ip2c1subs0",q2Bin));

  RooRealVar* mean_difference = (RooRealVar*)fitresult->floatParsFinal().find(Form("mean_difference^{%i}",year));

  RooGaussian* gauss_constr = new RooGaussian(Form("c_%s_%i", mean_difference->GetName(), year) ,
                                              Form("c_%s_%i", mean_difference->GetName(), year),
                                              *mean_difference,
                                              RooConst(mean_difference->getVal()),
                                              RooConst(mean_difference->getError())
                                              );

    if (addToList){
      c_vars.add(*mean_difference);
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
                                       RooArgSet reco_vars, std::string shortString, int comp, int dat){

    RooDataSet* data;
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

      if(dat == 0){
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
      }

      else if(dat == 1){ 
        data = (RooDataSet*)ws->data(Form("Dataset/data_b%i",q2Bin));
        isample->append(*data);
      }

      datasample.push_back (isample); 
    }
 
    return datasample;
}

void validate_fit(RooWorkspace* w, RooCategory sample, RooArgSet c_vars, int year, int q2Bin, int parity){

  RooRealVar mass = *(w->var("mass"));
  RooAbsPdf* model = w->pdf("simPdf");

  vector<RooRealVar> params;
  params.push_back(*(w->var( Form("sig_yield^{%i}",year) )));
  params.push_back(*(w->var( Form("mFrac^{%i}",year) )));
  params.push_back(*(w->var( Form("mean_{RT}^{%i}",year) )));
  params.push_back(*(w->var( Form("#sigma_{RT1}^{%i}",year) )));
  params.push_back(*(w->var( Form("#alpha_{RT1}^{%i}",year) )));
  params.push_back(*(w->var( Form("#alpha_{RT2}^{%i}",year) )));
  params.push_back(*(w->var( Form("n_{RT1}^{%i}",year) )));
  params.push_back(*(w->var( Form("n_{RT2}^{%i}",year) )));
  if(q2Bin >= 5){
    params.push_back(*(w->var( Form("#sigma_{RT2}^{%i}",year) )));
    params.push_back(*(w->var( Form("f^{RT%i}",year) )));
  }
  params.push_back(*(w->var( Form("mean_difference^{%i}",year) )));
  params.push_back(*(w->var( Form("#sigma_{WT1}^{%i}",year) )));
  params.push_back(*(w->var( Form("#alpha_{WT1}^{%i}",year) )));
  params.push_back(*(w->var( Form("#alpha_{WT2}^{%i}",year) )));
  params.push_back(*(w->var( Form("n_{WT1}^{%i}",year) )));
  params.push_back(*(w->var( Form("n_{WT2}^{%i}",year) )));

  string var_name;
  int params_size = params.size();

  cout << "Starting toy MC study" << endl;
  RooMCStudy* mcstudy = new RooMCStudy(*model, RooArgSet(mass,sample), Constrain(c_vars), Binned(kFALSE), Silence(kTRUE), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));
  mcstudy->generateAndFit(1000);
  cout << "Ending toy MC study" << endl;

  vector<RooPlot*> framesPull, framesParam;

  for(int i = 0; i < params_size; ++i){
    framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(100)));//200
    framesPull[i]->SetTitle("");
    framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
    framesParam[i]->SetTitle("");
  }
  
  vector<TGraph*> h;
  for(int i = 0; i < params_size; ++i){
    h.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
  }

  gStyle->SetOptFit(0111);

  TCanvas* c_pull;
  gPad->SetLeftMargin(0.15);

  ofstream output(Form("simFitPullResults/out_b%i_year_%i.dat", q2Bin, year));
  output << '|' << setw(15) <<  "var_name" << '|' << setw(15) <<  "  pull mean  " << '|' << setw(15) <<  "  pull width  " << endl;

  for(int i = 0; i < params_size; ++i){
    c_pull = new TCanvas( "pulls", "pulls", 900, 800);
    c_pull->cd();
    var_name = params[i].GetName();

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

    TF1* fit = h[i]->GetFunction("gaus");
    
    output << '|' << setw(15) <<  var_name << '|' << setw(15) <<  fit->GetParameter(1) << '|' << setw(15) <<  fit->GetParameter(2) << endl;
    
    c_pull->SaveAs(( "plotSimMassFit_pulls/pulls_poisson_" + var_name + Form("b%ip%i_%i",q2Bin,parity,year) + ".gif").c_str() );
    delete c_pull;
  }
  output.close();

  TCanvas* c_params = new TCanvas("params", "params", 900, 800);

  for(int i = 0; i < params_size; ++i){
    var_name = params[i].GetName();

    c_params->cd();
    framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
    framesParam.at(i)->Draw();

    c_params->SaveAs(( "plotSimMassFit_pulls/pulls_params_poisson_" + var_name + Form("b%ip%i_%i",q2Bin,parity,year) + ".gif").c_str() );
  }

}


