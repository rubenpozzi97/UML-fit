#include <TH1F.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <TAttMarker.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>

double read_weights(TH1F* histo_variable, double var_value);
double getWeight(double var_value, TH1F* h_weight);

void efficiency(int year){

  double q2_bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  const int n_q2_bins = 8;

  if(year < 2016 || year > 2018){return;}

  // RECO MC (PUweight)
  TString input_file_mc_cuts_jpsi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_JPSI.root",year,year);
  TString input_file_mc_cuts_psi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_PSI.root",year,year);
  TString input_file_mc_cuts_lmnr = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_LMNR.root",year,year);
  TFile* f_mc_cuts_jpsi = new TFile(input_file_mc_cuts_jpsi);
  TFile* f_mc_cuts_psi = new TFile(input_file_mc_cuts_psi);
  TFile* f_mc_cuts_lmnr = new TFile(input_file_mc_cuts_lmnr);

  TString input_file_mc_nocuts_jpsi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_JPSI.root",year,year);
  TString input_file_mc_nocuts_psi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_PSI.root",year,year);
  TString input_file_mc_nocuts_lmnr = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_LMNR.root",year,year);
  TFile* f_mc_nocuts_jpsi = new TFile(input_file_mc_nocuts_jpsi);
  TFile* f_mc_nocuts_psi = new TFile(input_file_mc_nocuts_psi);
  TFile* f_mc_nocuts_lmnr = new TFile(input_file_mc_nocuts_lmnr);
 
  TTree* t_cuts_jpsi = (TTree*)f_mc_cuts_jpsi->Get("ntuple");
  TTree* t_cuts_psi = (TTree*)f_mc_cuts_psi->Get("ntuple");
  TTree* t_cuts_lmnr = (TTree*)f_mc_cuts_lmnr->Get("ntuple");

  TTree* t_nocuts_jpsi = (TTree*)f_mc_nocuts_jpsi->Get("ntuple");
  TTree* t_nocuts_psi = (TTree*)f_mc_nocuts_psi->Get("ntuple");
  TTree* t_nocuts_lmnr = (TTree*)f_mc_nocuts_lmnr->Get("ntuple");

  // GEN-FILTER MC (acceptance)
  TString input_file_mc_gen_jpsi = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0JpsiKstar.root";
  TString input_file_mc_gen_psi = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0PsiKstar.root";
  TString input_file_mc_gen_lmnr = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0MuMuKstar_p*.root/ntuple";
  TFile* f_mc_gen_jpsi = new TFile(input_file_mc_gen_jpsi);
  TFile* f_mc_gen_psi = new TFile(input_file_mc_gen_psi);
  TChain* t_gen_lmnr = new TChain();

  TTree* t_gen_jpsi = (TTree*)f_mc_gen_jpsi->Get("ntuple");
  TTree* t_gen_psi = (TTree*)f_mc_gen_psi->Get("ntuple");
  t_gen_lmnr->Add(input_file_mc_gen_lmnr);

  // data(splot)/MC weights for systematics
  TFile* weight_b0 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b0.root",year));
  TFile* weight_b1 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b1.root",year));
  TFile* weight_b2 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b2.root",year));
  TFile* weight_b3 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b3.root",year));
  TFile* weight_b4 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b4.root",year));
  TFile* weight_b5 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b5.root",year));
  TFile* weight_b6 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b6.root",year));
  TFile* weight_b7 = new TFile(Form("~/public/UML-fit/results/mc_validation_plots/weights/weights_%i_b7.root",year));

  TH1F* histo_wei_b0 = (TH1F*)weight_b0->Get("weights_bdt_prob");
  TH1F* histo_wei_b1 = (TH1F*)weight_b1->Get("weights_bdt_prob");
  TH1F* histo_wei_b2 = (TH1F*)weight_b2->Get("weights_bdt_prob");
  TH1F* histo_wei_b3 = (TH1F*)weight_b3->Get("weights_bdt_prob");
  TH1F* histo_wei_b4 = (TH1F*)weight_b4->Get("weights_bdt_prob");
  TH1F* histo_wei_b5 = (TH1F*)weight_b5->Get("weights_bdt_prob");
  TH1F* histo_wei_b6 = (TH1F*)weight_b6->Get("weights_bdt_prob");
  TH1F* histo_wei_b7 = (TH1F*)weight_b7->Get("weights_bdt_prob");

  double mumuMass_jpsi;
  double mumuMass_psi;
  double mumuMass_lmnr;

  double gen_mumuMass_jpsi;
  double gen_mumuMass_psi;
  double gen_mumuMass_lmnr;

  double bfilter_mumuMass_jpsi;
  double bfilter_mumuMass_psi;
  double bfilter_mumuMass_lmnr;

  //variables for applying GEN-filter
  double genmupEta_jpsi;
  double genmumEta_jpsi;
  double genkstTrkpEta_jpsi;
  double genkstTrkmEta_jpsi;
  double genmupPt_jpsi;
  double genmumPt_jpsi;
  double genkstTrkpPt_jpsi;
  double genkstTrkmPt_jpsi;

  double genmupEta_psi;
  double genmumEta_psi;
  double genkstTrkpEta_psi;
  double genkstTrkmEta_psi;
  double genmupPt_psi;
  double genmumPt_psi;
  double genkstTrkpPt_psi;
  double genkstTrkmPt_psi;

  double genmupEta_lmnr;
  double genmumEta_lmnr;
  double genkstTrkpEta_lmnr;
  double genkstTrkmEta_lmnr;
  double genmupPt_lmnr;
  double genmumPt_lmnr;
  double genkstTrkpPt_lmnr;
  double genkstTrkmPt_lmnr;

  double genbEta_lmnr;
  double genbEta_jpsi;
  double genbEta_psi;

  double genmupEta_jpsi1;
  double genmumEta_jpsi1;
  double genkstTrkpEta_jpsi1;
  double genkstTrkmEta_jpsi1;
  double genmupPt_jpsi1;
  double genmumPt_jpsi1;
  double genkstTrkpPt_jpsi1;
  double genkstTrkmPt_jpsi1;

  double genmupEta_psi1;
  double genmumEta_psi1;
  double genkstTrkpEta_psi1;
  double genkstTrkmEta_psi1;
  double genmupPt_psi1;
  double genmumPt_psi1;
  double genkstTrkpPt_psi1;
  double genkstTrkmPt_psi1;

  double genmupEta_lmnr1;
  double genmumEta_lmnr1;
  double genkstTrkpEta_lmnr1;
  double genkstTrkmEta_lmnr1;
  double genmupPt_lmnr1;
  double genmumPt_lmnr1;
  double genkstTrkpPt_lmnr1;
  double genkstTrkmPt_lmnr1;

  double genbEta_lmnr1;
  double genbEta_jpsi1;
  double genbEta_psi1;

  t_cuts_jpsi->SetBranchAddress("mumuMass",&mumuMass_jpsi);
  t_cuts_psi->SetBranchAddress("mumuMass",&mumuMass_psi);
  t_cuts_lmnr->SetBranchAddress("mumuMass",&mumuMass_lmnr);

  t_nocuts_jpsi->SetBranchAddress("genQ",&gen_mumuMass_jpsi);
  t_nocuts_psi->SetBranchAddress("genQ",&gen_mumuMass_psi);
  t_nocuts_lmnr->SetBranchAddress("genQ",&gen_mumuMass_lmnr);

  t_gen_jpsi->SetBranchAddress("genQ",&bfilter_mumuMass_jpsi);
  t_gen_psi->SetBranchAddress("genQ",&bfilter_mumuMass_psi);
  t_gen_lmnr->SetBranchAddress("genQ",&bfilter_mumuMass_lmnr);

  t_gen_jpsi->SetBranchAddress("genmupEta",&genmupEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genmumEta",&genmumEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_jpsi);
  t_gen_jpsi->SetBranchAddress("genmupPt",&genmupPt_jpsi);
  t_gen_jpsi->SetBranchAddress("genmumPt",&genmumPt_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_jpsi);
  t_gen_jpsi->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_jpsi);

  t_gen_psi->SetBranchAddress("genmupEta",&genmupEta_psi);
  t_gen_psi->SetBranchAddress("genmumEta",&genmumEta_psi);
  t_gen_psi->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_psi);
  t_gen_psi->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_psi);
  t_gen_psi->SetBranchAddress("genmupPt",&genmupPt_psi);
  t_gen_psi->SetBranchAddress("genmumPt",&genmumPt_psi);
  t_gen_psi->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_psi);
  t_gen_psi->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_psi);

  t_gen_lmnr->SetBranchAddress("genmupEta",&genmupEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genmumEta",&genmumEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_lmnr);
  t_gen_lmnr->SetBranchAddress("genmupPt",&genmupPt_lmnr);
  t_gen_lmnr->SetBranchAddress("genmumPt",&genmumPt_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_lmnr);
  t_gen_lmnr->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_lmnr);

  t_gen_lmnr->SetBranchAddress("genbEta",&genbEta_lmnr);
  t_gen_jpsi->SetBranchAddress("genbEta",&genbEta_jpsi);
  t_gen_psi->SetBranchAddress("genbEta",&genbEta_psi);

  t_nocuts_jpsi->SetBranchAddress("genmupEta",&genmupEta_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genmumEta",&genmumEta_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genmupPt",&genmupPt_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genmumPt",&genmumPt_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_jpsi1);
  t_nocuts_jpsi->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_jpsi1);

  t_nocuts_psi->SetBranchAddress("genmupEta",&genmupEta_psi1);
  t_nocuts_psi->SetBranchAddress("genmumEta",&genmumEta_psi1);
  t_nocuts_psi->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_psi1);
  t_nocuts_psi->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_psi1);
  t_nocuts_psi->SetBranchAddress("genmupPt",&genmupPt_psi1);
  t_nocuts_psi->SetBranchAddress("genmumPt",&genmumPt_psi1);
  t_nocuts_psi->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_psi1);
  t_nocuts_psi->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_psi1);

  t_nocuts_lmnr->SetBranchAddress("genmupEta",&genmupEta_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genmumEta",&genmumEta_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genkstTrkpEta",&genkstTrkpEta_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genkstTrkmEta",&genkstTrkmEta_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genmupPt",&genmupPt_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genmumPt",&genmumPt_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genkstTrkpPt",&genkstTrkpPt_lmnr1);
  t_nocuts_lmnr->SetBranchAddress("genkstTrkmPt",&genkstTrkmPt_lmnr1);

  t_nocuts_lmnr->SetBranchAddress("genbEta",&genbEta_lmnr1);
  t_nocuts_jpsi->SetBranchAddress("genbEta",&genbEta_jpsi1);
  t_nocuts_psi->SetBranchAddress("genbEta",&genbEta_psi1);

  double reco_entries[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double gen_entries[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double den_entries[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double num_entries[] = {0., 0., 0., 0., 0., 0., 0., 0.};

  // systematics 
  double reco_entries_weight[] = {0., 0., 0., 0., 0., 0., 0., 0.};

  float reco_bdt_lmnr;
  float reco_bdt_jpsi;
  float reco_bdt_psi;

  t_cuts_lmnr->SetBranchAddress("bdt_prob", &reco_bdt_lmnr);
  t_cuts_jpsi->SetBranchAddress("bdt_prob", &reco_bdt_jpsi);
  t_cuts_psi->SetBranchAddress("bdt_prob", &reco_bdt_psi);

  cout << "RECO JPSI"  << endl; 
  for(int evt = 0; evt < t_cuts_jpsi->GetEntries(); evt++){
    t_cuts_jpsi->GetEntry(evt);
     
    if( (pow(mumuMass_jpsi,2) > q2_bins[4]) && (pow(mumuMass_jpsi,2) < q2_bins[5]) ){
      reco_entries[4] += 1;
      reco_entries_weight[4] += read_weights(histo_wei_b4,reco_bdt_jpsi);
    }
  }

  cout << "RECO PSI"  << endl;
  for(int evt = 0; evt < t_cuts_psi->GetEntries(); evt++){
    t_cuts_psi->GetEntry(evt);

    if( (pow(mumuMass_psi,2) > q2_bins[6]) && (pow(mumuMass_psi,2) < q2_bins[7]) ){
      reco_entries[6] += 1;
      reco_entries_weight[6] += read_weights(histo_wei_b6,reco_bdt_psi);
    }
  }

  cout << "RECO LMNR"  << endl;
  for(int evt = 0; evt < t_cuts_lmnr->GetEntries(); evt++){
    t_cuts_lmnr->GetEntry(evt);

    if( (pow(mumuMass_lmnr,2) > q2_bins[0]) && (pow(mumuMass_lmnr,2) < q2_bins[1]) ){
      reco_entries[0] += 1;
      reco_entries_weight[0] += read_weights(histo_wei_b0,reco_bdt_lmnr);
    }
    else if( (pow(mumuMass_lmnr,2) > q2_bins[1]) && (pow(mumuMass_lmnr,2) < q2_bins[2]) ){
      reco_entries[1] += 1;
      reco_entries_weight[1] += read_weights(histo_wei_b1,reco_bdt_lmnr);
    }
    else if( (pow(mumuMass_lmnr,2) > q2_bins[2]) && (pow(mumuMass_lmnr,2) < q2_bins[3]) ){
      reco_entries[2] += 1;
      reco_entries_weight[2] += read_weights(histo_wei_b2,reco_bdt_lmnr);
    }
    else if( (pow(mumuMass_lmnr,2) > q2_bins[3]) && (pow(mumuMass_lmnr,2) < q2_bins[4]) ){
      reco_entries[3] += 1;
      reco_entries_weight[3] += read_weights(histo_wei_b3,reco_bdt_lmnr);
    }
    else if( (pow(mumuMass_lmnr,2) > q2_bins[5]) && (pow(mumuMass_lmnr,2) < q2_bins[6]) ){
      reco_entries[5] += 1;
      reco_entries_weight[5] += read_weights(histo_wei_b5,reco_bdt_lmnr);
    }
    else if( (pow(mumuMass_lmnr,2) > q2_bins[7]) && (pow(mumuMass_lmnr,2) < q2_bins[8]) ){
      reco_entries[7] += 1;
      reco_entries_weight[7] += read_weights(histo_wei_b7,reco_bdt_lmnr);
    }
  }

  cout << "GEN JPSI" << endl;
  for(int evt = 0; evt < t_nocuts_jpsi->GetEntries(); evt++){
      t_nocuts_jpsi->GetEntry(evt);

      if( (pow(gen_mumuMass_jpsi,2) > q2_bins[4]) && (pow(gen_mumuMass_jpsi,2) < q2_bins[5]) ){
        if(genbEta_jpsi1 < 3){
          if(fabs(genmupEta_jpsi1)<2.5 && fabs(genmumEta_jpsi1)<2.5 &&
            fabs(genkstTrkpEta_jpsi1)<2.5 && fabs(genkstTrkmEta_jpsi1)<2.5 &&
            genmupPt_jpsi1>2.5 && genmumPt_jpsi1>2.5 &&
            genkstTrkpPt_jpsi1>0.4 && genkstTrkmPt_jpsi1>0.4){
              gen_entries[4] += 1;
          }
        }
      }
  }

  cout << "GEN PSI" << endl;
  for(int evt = 0; evt < t_nocuts_psi->GetEntries(); evt++){
      t_nocuts_psi->GetEntry(evt);

      if( (pow(gen_mumuMass_psi,2) > q2_bins[6]) && (pow(gen_mumuMass_psi,2) < q2_bins[7]) ){
        if(genbEta_psi1 < 3){
          if(fabs(genmupEta_psi1)<2.5 && fabs(genmumEta_psi1)<2.5 &&
            fabs(genkstTrkpEta_psi1)<2.5 && fabs(genkstTrkmEta_psi1)<2.5 &&
            genmupPt_psi1>2.5 && genmumPt_psi1>2.5 &&
            genkstTrkpPt_psi1>0.4 && genkstTrkmPt_psi1>0.4){
              gen_entries[6] += 1;
          }
        }
      }
  }

  cout << "GEN LMNR" << endl;
  for(int evt = 0; evt < t_nocuts_lmnr->GetEntries(); evt++){
      t_nocuts_lmnr->GetEntry(evt);

    if( (pow(gen_mumuMass_lmnr,2) > q2_bins[0]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[1]) ){
      if(genbEta_lmnr1 < 3){
        if(fabs(genmupEta_lmnr1)<2.5 && fabs(genmumEta_lmnr1)<2.5 &&
          fabs(genkstTrkpEta_lmnr1)<2.5 && fabs(genkstTrkmEta_lmnr1)<2.5 &&
          genmupPt_lmnr1>2.5 && genmumPt_lmnr1>2.5 &&
          genkstTrkpPt_lmnr1>0.4 && genkstTrkmPt_lmnr1>0.4){
            gen_entries[0] += 1;
        }
      }
    }

    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[1]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[2]) ){
      if(genbEta_lmnr1 < 3){
        if(fabs(genmupEta_lmnr1)<2.5 && fabs(genmumEta_lmnr1)<2.5 &&
          fabs(genkstTrkpEta_lmnr1)<2.5 && fabs(genkstTrkmEta_lmnr1)<2.5 &&
          genmupPt_lmnr1>2.5 && genmumPt_lmnr1>2.5 &&
          genkstTrkpPt_lmnr1>0.4 && genkstTrkmPt_lmnr1>0.4){
            gen_entries[1] += 1;
        }
      }
    }
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[2]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[3]) ){
      if(genbEta_lmnr1 < 3){
        if(fabs(genmupEta_lmnr1)<2.5 && fabs(genmumEta_lmnr1)<2.5 &&
          fabs(genkstTrkpEta_lmnr1)<2.5 && fabs(genkstTrkmEta_lmnr1)<2.5 &&
          genmupPt_lmnr1>2.5 && genmumPt_lmnr1>2.5 &&
          genkstTrkpPt_lmnr1>0.4 && genkstTrkmPt_lmnr1>0.4){
            gen_entries[2] += 1;
        }
      }
    }
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[3]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[4]) ){
      if(genbEta_lmnr1 < 3){
        if(fabs(genmupEta_lmnr1)<2.5 && fabs(genmumEta_lmnr1)<2.5 &&
          fabs(genkstTrkpEta_lmnr1)<2.5 && fabs(genkstTrkmEta_lmnr1)<2.5 &&
          genmupPt_lmnr1>2.5 && genmumPt_lmnr1>2.5 &&
          genkstTrkpPt_lmnr1>0.4 && genkstTrkmPt_lmnr1>0.4){
            gen_entries[3] += 1;
        }
      }
    }
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[5]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[6]) ){
      if(genbEta_lmnr1 < 3){
        if(fabs(genmupEta_lmnr1)<2.5 && fabs(genmumEta_lmnr1)<2.5 &&
          fabs(genkstTrkpEta_lmnr1)<2.5 && fabs(genkstTrkmEta_lmnr1)<2.5 &&
          genmupPt_lmnr1>2.5 && genmumPt_lmnr1>2.5 &&
          genkstTrkpPt_lmnr1>0.4 && genkstTrkmPt_lmnr1>0.4){
            gen_entries[5] += 1;
        }
      }
    }
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[7]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[8]) ){
      if(genbEta_lmnr1 < 3){
        if(fabs(genmupEta_lmnr1)<2.5 && fabs(genmumEta_lmnr1)<2.5 &&
          fabs(genkstTrkpEta_lmnr1)<2.5 && fabs(genkstTrkmEta_lmnr1)<2.5 &&
          genmupPt_lmnr1>2.5 && genmumPt_lmnr1>2.5 &&
          genkstTrkpPt_lmnr1>0.4 && genkstTrkmPt_lmnr1>0.4){
            gen_entries[7] += 1;
        }
      } 
    }
  }

  cout << "GEN B-FILTER JPSI" << endl;
  for(int evt = 0; evt < t_gen_jpsi->GetEntries(); evt++){
      t_gen_jpsi->GetEntry(evt);

      if( (pow(bfilter_mumuMass_jpsi,2) > q2_bins[4]) && (pow(bfilter_mumuMass_jpsi,2) < q2_bins[5]) ){
        if(genbEta_jpsi < 3){
          den_entries[4] += 1;
          if(fabs(genmupEta_jpsi)<2.5 && fabs(genmumEta_jpsi)<2.5 &&
             fabs(genkstTrkpEta_jpsi)<2.5 && fabs(genkstTrkmEta_jpsi)<2.5 &&
             genmupPt_jpsi>2.5 && genmumPt_jpsi>2.5 &&
             genkstTrkpPt_jpsi>0.4 && genkstTrkmPt_jpsi>0.4){
               num_entries[4] += 1;
          }
        }
      }
  }

  cout << "GEN B-FILTER PSI" << endl;
  for(int evt = 0; evt < t_gen_psi->GetEntries(); evt++){
      t_gen_psi->GetEntry(evt);

      if( (pow(bfilter_mumuMass_psi,2) > q2_bins[6]) && (pow(bfilter_mumuMass_psi,2) < q2_bins[7]) ){
        if(genbEta_psi < 3){
          den_entries[6] += 1;
          if(fabs(genmupEta_psi)<2.5 && fabs(genmumEta_psi)<2.5 &&
             fabs(genkstTrkpEta_psi)<2.5 && fabs(genkstTrkmEta_psi)<2.5 &&
             genmupPt_psi>2.5 && genmumPt_psi>2.5 &&
             genkstTrkpPt_psi>0.4 && genkstTrkmPt_psi>0.4){
               num_entries[6] += 1;
          }
        }
      }
  }

  cout << "GEN B-FILTER LMNR" << endl;
  for(int evt = 0; evt < t_gen_lmnr->GetEntries(); evt++){
      t_gen_lmnr->GetEntry(evt);

      if( (pow(bfilter_mumuMass_lmnr,2) > q2_bins[0]) && (pow(bfilter_mumuMass_lmnr,2) < q2_bins[1]) ){
        if(genbEta_lmnr < 3){
          den_entries[0] += 1;
          if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
             fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
             genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
             genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
               num_entries[0] += 1;
          }
        }
      }
      else if( (pow(bfilter_mumuMass_lmnr,2) > q2_bins[1]) && (pow(bfilter_mumuMass_lmnr,2) < q2_bins[2]) ){
        if(genbEta_lmnr < 3){
          den_entries[1] += 1;
          if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
             fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
             genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
             genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
               num_entries[1] += 1;
          }
        }
      }
      else if( (pow(bfilter_mumuMass_lmnr,2) > q2_bins[2]) && (pow(bfilter_mumuMass_lmnr,2) < q2_bins[3]) ){
        if(genbEta_lmnr < 3){
          den_entries[2] += 1;
          if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
             fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
             genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
             genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
               num_entries[2] += 1;
          }
        }
      }
      else if( (pow(bfilter_mumuMass_lmnr,2) > q2_bins[3]) && (pow(bfilter_mumuMass_lmnr,2) < q2_bins[4]) ){
        if(genbEta_lmnr < 3){
          den_entries[3] += 1;
          if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
             fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
             genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
             genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
               num_entries[3] += 1;
          }
        }
      }
      else if( (pow(bfilter_mumuMass_lmnr,2) > q2_bins[5]) && (pow(bfilter_mumuMass_lmnr,2) < q2_bins[6]) ){
        if(genbEta_lmnr < 3){
          den_entries[5] += 1;
          if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
             fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
             genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
             genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
               num_entries[5] += 1;
          }
        }
      }
      else if( (pow(bfilter_mumuMass_lmnr,2) > q2_bins[7]) && (pow(bfilter_mumuMass_lmnr,2) < q2_bins[8]) ){
        if(genbEta_lmnr < 3){
          den_entries[7] += 1;
          if(fabs(genmupEta_lmnr)<2.5 && fabs(genmumEta_lmnr)<2.5 &&
             fabs(genkstTrkpEta_lmnr)<2.5 && fabs(genkstTrkmEta_lmnr)<2.5 &&
             genmupPt_lmnr>2.5 && genmumPt_lmnr>2.5 &&
             genkstTrkpPt_lmnr>0.4 && genkstTrkmPt_lmnr>0.4){
               num_entries[7] += 1;
          }
        }
      }
  }

  cout << '|' << setw(15) << "q2Bin" << '|' << setw(15) << "Efficiency" << '|' << setw(15) << "Acceptance" << '|' << setw(15) << "Efficiency x acceptance" << '|' << setw(15) << "Weighted eff" << '|' << setw(15) << "syst" << '|' << endl;

  double eff[n_q2_bins];
  double acc[n_q2_bins];
  double wei_eff[n_q2_bins];

  double eff_x_acc[n_q2_bins];
  double syst[n_q2_bins]; 
  double absolute_syst[n_q2_bins];

  for(int i = 0; i < n_q2_bins; i++){
    eff[i] = reco_entries[i]/gen_entries[i];
    acc[i] = num_entries[i]/den_entries[i];
    eff_x_acc[i] = eff[i]*acc[i];

    wei_eff[i] = reco_entries_weight[i]/gen_entries[i];

    syst[i] = abs(wei_eff[i]-eff[i])/eff[i];
    absolute_syst[i] = syst[i]*eff_x_acc[i];

  cout << '|' << setw(15) << q2_bins[i] << " - " << q2_bins[i+1] << '|' << setw(15) << eff[i] << '|' << setw(15) << acc[i] << '|' << setw(15) << eff_x_acc[i] << '|' << setw(15) << wei_eff[i] << '|' << setw(15) << syst[i] << '|' << endl;
  }

  double q2Bins_half[n_q2_bins];
  double q2Bins_err[n_q2_bins];

  for(int i = 0; i < n_q2_bins; i++){
    q2Bins_half[i] = q2_bins[i] + 0.5* (q2_bins[i+1]-q2_bins[i]);
    q2Bins_err[i] = 0.5* (q2_bins[i+1]-q2_bins[i]);
  }

  TGraph* geff = new TGraph(n_q2_bins,q2Bins_half,eff);
  TGraph* gacc = new TGraph(n_q2_bins,q2Bins_half,acc);
  TGraph* gwei = new TGraph(n_q2_bins,q2Bins_half,wei_eff);
  TGraphErrors* g_eff = new TGraphErrors(n_q2_bins,q2Bins_half,eff_x_acc,q2Bins_err,absolute_syst);

  TCanvas c;
  c.cd();
  g_eff->SetTitle(Form("Efficiency x Acceptance - %i",year));
  g_eff->GetXaxis()->SetTitle("q^{2} [GeV]");
  g_eff->GetYaxis()->SetTitle("#epsilon");
  g_eff->SetMarkerStyle(8);
  g_eff->SetMarkerColor(1);
  g_eff->Draw("AP");
  c.SaveAs(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.gif",year));

  TFile* f;
  f = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",year),"UPDATE");
  f->cd();
  geff->Write();
  gacc->Write();
  gwei->Write();
  g_eff->Write();
  f->Close();

}

double read_weights(TH1F* histo_variable, double var_value){

  double weight;
  double variable_min;
  double variable_max;

  variable_min = histo_variable->GetXaxis()->GetXmin();
  variable_max = histo_variable->GetXaxis()->GetXmax();

  //if the event is not in the range its weight is 1.
  if(var_value>=variable_min && var_value<=variable_max){weight = getWeight(var_value,histo_variable);}
  else{weight = 1;}

  return weight;
}

double getWeight(double var_value, TH1F* h_weight){
  int bin = h_weight->FindBin(var_value);
  return h_weight->GetBinContent(bin);
}

