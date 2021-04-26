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

void efficiency(int year, int q2Bin = -1){

  double q2_bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  const int n_q2_bins = 8;

  if(q2Bin < -1 || q2Bin >= n_q2_bins){return;}
  if(year < 2016 || year > 2018){return;}


  // RECO MC
  TString input_file_mc_cuts_jpsi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_JPSI.root",year,year);
  TString input_file_mc_cuts_psi = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_PSI.root",year,year);
  TString input_file_mc_cuts_lmnr = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_LMNR.root",year,year);
  TFile* f_mc_cuts_jpsi = new TFile(input_file_mc_cuts_jpsi);
  TFile* f_mc_cuts_psi = new TFile(input_file_mc_cuts_psi);
  TFile* f_mc_cuts_lmnr = new TFile(input_file_mc_cuts_lmnr);

  // GEN-LEVEL MC
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

  double mumuMass_jpsi;
  double mumuMass_psi;
  double mumuMass_lmnr;

  double gen_mumuMass_jpsi;
  double gen_mumuMass_psi;
  double gen_mumuMass_lmnr;

  t_cuts_jpsi->SetBranchAddress("mumuMass",&mumuMass_jpsi);
  t_cuts_psi->SetBranchAddress("mumuMass",&mumuMass_psi);
  t_cuts_lmnr->SetBranchAddress("mumuMass",&mumuMass_lmnr);

  t_nocuts_jpsi->SetBranchAddress("genQ",&gen_mumuMass_jpsi);
  t_nocuts_psi->SetBranchAddress("genQ",&gen_mumuMass_psi);
  t_nocuts_lmnr->SetBranchAddress("genQ",&gen_mumuMass_lmnr);

  double reco_entries[] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double gen_entries[] = {0., 0., 0., 0., 0., 0., 0., 0.};

  cout << "RECO JPSI"  << endl; 
  for(int evt = 0; evt < t_cuts_jpsi->GetEntries(); evt++){
    t_cuts_jpsi->GetEntry(evt);
     
    if( (pow(mumuMass_jpsi,2) > q2_bins[4]) && (pow(mumuMass_jpsi,2) < q2_bins[5]) ){reco_entries[4] += 1;}
  }

  cout << "RECO PSI"  << endl;
  for(int evt = 0; evt < t_cuts_psi->GetEntries(); evt++){
    t_cuts_psi->GetEntry(evt);

    if( (pow(mumuMass_psi,2) > q2_bins[6]) && (pow(mumuMass_psi,2) < q2_bins[7]) ){reco_entries[6] += 1;}
  }

  cout << "RECO LMNR"  << endl;
  for(int evt = 0; evt < t_cuts_lmnr->GetEntries(); evt++){
    t_cuts_lmnr->GetEntry(evt);

    if( (pow(mumuMass_lmnr,2) > q2_bins[0]) && (pow(mumuMass_lmnr,2) < q2_bins[1]) ){reco_entries[0] += 1;}
    else if( (pow(mumuMass_lmnr,2) > q2_bins[1]) && (pow(mumuMass_lmnr,2) < q2_bins[2]) ){reco_entries[1] += 1;}
    else if( (pow(mumuMass_lmnr,2) > q2_bins[2]) && (pow(mumuMass_lmnr,2) < q2_bins[3]) ){reco_entries[2] += 1;}
    else if( (pow(mumuMass_lmnr,2) > q2_bins[3]) && (pow(mumuMass_lmnr,2) < q2_bins[4]) ){reco_entries[3] += 1;}
    else if( (pow(mumuMass_lmnr,2) > q2_bins[5]) && (pow(mumuMass_lmnr,2) < q2_bins[6]) ){reco_entries[5] += 1;}
    else if( (pow(mumuMass_lmnr,2) > q2_bins[7]) && (pow(mumuMass_lmnr,2) < q2_bins[8]) ){reco_entries[7] += 1;}
  }

  cout << "GEN JPSI" << endl;
  for(int evt = 0; evt < t_nocuts_jpsi->GetEntries(); evt++){
      t_nocuts_jpsi->GetEntry(evt);

      if( (pow(gen_mumuMass_jpsi,2) > q2_bins[4]) && (pow(gen_mumuMass_jpsi,2) < q2_bins[5]) ){gen_entries[4] += 1;}
  }

  cout << "GEN PSI" << endl;
  for(int evt = 0; evt < t_nocuts_psi->GetEntries(); evt++){
      t_nocuts_psi->GetEntry(evt);

      if( (pow(gen_mumuMass_psi,2) > q2_bins[6]) && (pow(gen_mumuMass_psi,2) < q2_bins[7]) ){gen_entries[6] += 1;}
  }

  cout << "GEN LMNR" << endl;
  for(int evt = 0; evt < t_nocuts_lmnr->GetEntries(); evt++){
      t_nocuts_lmnr->GetEntry(evt);

    if( (pow(gen_mumuMass_lmnr,2) > q2_bins[0]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[1]) ){gen_entries[0] += 1;}
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[1]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[2]) ){gen_entries[1] += 1;}
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[2]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[3]) ){gen_entries[2] += 1;}
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[3]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[4]) ){gen_entries[3] += 1;}
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[5]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[6]) ){gen_entries[5] += 1;}
    else if( (pow(gen_mumuMass_lmnr,2) > q2_bins[7]) && (pow(gen_mumuMass_lmnr,2) < q2_bins[8]) ){gen_entries[7] += 1;}
  }


  cout << '|' << setw(15) << "q2Bin" << '|' << setw(15) << "Passed" << '|' << setw(15) << "Total" << '|' << endl;

  TH1F* eff_x_acc = new TH1F("eff_x_acc","eff_x_acc",n_q2_bins,q2_bins);

  for(int i = 0; i < n_q2_bins; i++){
    cout << '|' << setw(15) << i << '|' << setw(15) << reco_entries[i] << '|' << setw(15) << gen_entries[i] << '|' << endl;
    eff_x_acc->SetBinContent(i,reco_entries[i]/gen_entries[i]);
  }

  TCanvas c;
  c.cd();
  eff_x_acc->SetTitle(Form("Efficiency x Acceptance - %i",year));
  eff_x_acc->GetXaxis()->SetTitle("q^2 [GeV]");
  eff_x_acc->GetYaxis()->SetTitle("#epsilon");
  eff_x_acc->SetMarkerStyle(8);
  eff_x_acc->SetMarkerColor(1);
  eff_x_acc->Draw("P");
  c.SaveAs(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.gif",year));

  TFile* f;
  f = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",year),"UPDATE");
  f->cd();
  eff_x_acc->Write();
  f->Close();
}




