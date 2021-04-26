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

void efficiency(int year){

  double q2_bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  const int n_q2_bins = 8;

  TH1F* eff_x_acc = new TH1F("eff_x_acc","eff_x_acc",n_q2_bins,q2_bins);

  double bin0[] = {1,2};
  double bin1[] = {2,4.3};
  double bin2[] = {4.3,6};
  double bin3[] = {6,8.68};
  double bin4[] = {8.68,10.09};
  double bin5[] = {10.09,12.86};
  double bin6[] = {12.86,14.18};
  double bin7[] = {14.18,16};

  TH1F* pass0 = new TH1F("pass0","pass0",1,bin0);
  TH1F* pass1 = new TH1F("pass1","pass1",1,bin1);
  TH1F* pass2 = new TH1F("pass2","pass2",1,bin2);
  TH1F* pass3 = new TH1F("pass3","pass3",1,bin3);
  TH1F* pass4 = new TH1F("pass4","pass4",1,bin4);
  TH1F* pass5 = new TH1F("pass5","pass5",1,bin5);
  TH1F* pass6 = new TH1F("pass6","pass6",1,bin6);
  TH1F* pass7 = new TH1F("pass7","pass7",1,bin7);

  TH1F* gen0 = new TH1F("gen0","gen0",1,bin0);
  TH1F* gen1 = new TH1F("gen1","gen1",1,bin1);
  TH1F* gen2 = new TH1F("gen2","gen2",1,bin2);
  TH1F* gen3 = new TH1F("gen3","gen3",1,bin3);
  TH1F* gen4 = new TH1F("gen4","gen4",1,bin4);
  TH1F* gen5 = new TH1F("gen5","gen5",1,bin5);
  TH1F* gen6 = new TH1F("gen6","gen6",1,bin6);
  TH1F* gen7 = new TH1F("gen7","gen7",1,bin7);

  TEfficiency* eff0;
  TEfficiency* eff1;
  TEfficiency* eff2;
  TEfficiency* eff3;
  TEfficiency* eff4;
  TEfficiency* eff5;
  TEfficiency* eff6;
  TEfficiency* eff7;

  for(int q2Bin = 0; q2Bin < n_q2_bins; q2Bin++){

    // RECO MC
    TString input_file_mc_cuts;
    if(q2Bin == 4){input_file_mc_cuts = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_JPSI.root",year,year);}
    else if(q2Bin == 6){input_file_mc_cuts = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_PSI.root",year,year);}
    else{input_file_mc_cuts = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iMC_LMNR.root",year,year);}
    TFile* f_mc_cuts = new TFile(input_file_mc_cuts);

    // GEN-LEVEL MC
    TString input_file_mc_nocuts;
    if(q2Bin == 4){input_file_mc_nocuts = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_JPSI.root",year,year);}
    else if(q2Bin == 6){input_file_mc_nocuts = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_PSI.root",year,year);}
    else{input_file_mc_nocuts = Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iGEN_MC_LMNR.root",year,year);}
    TFile* f_mc_nocuts = new TFile(input_file_mc_nocuts);
 
    TTree* t_cuts = (TTree*)f_mc_cuts->Get("ntuple");
    TTree* t_nocuts = (TTree*)f_mc_nocuts->Get("ntuple");

    double mumuMass;
    double gen_mumuMass;

    t_cuts->SetBranchAddress("mumuMass",&mumuMass);
    t_nocuts->SetBranchAddress("genQ",&gen_mumuMass);

    cout << "RECO MC - " << "q2Bin = " << q2Bin << " - year = " << year << endl; 
    for(int evt = 0; evt < t_cuts->GetEntries(); evt++){
      t_cuts->GetEntry(evt);
     
      if( (q2Bin == 0) && (pow(mumuMass,2) > q2_bins[0]) && (pow(mumuMass,2) < q2_bins[1]) ){pass0->Fill(mumuMass);}
      else if( (q2Bin == 1) && (pow(mumuMass,2) > q2_bins[1]) && (pow(mumuMass,2) < q2_bins[2]) ){pass1->Fill(mumuMass);}
      else if( (q2Bin == 2) && (pow(mumuMass,2) > q2_bins[2]) && (pow(mumuMass,2) < q2_bins[3]) ){pass2->Fill(mumuMass);}
      else if( (q2Bin == 3) && (pow(mumuMass,2) > q2_bins[3]) && (pow(mumuMass,2) < q2_bins[4]) ){pass3->Fill(mumuMass);}
      else if( (q2Bin == 4) && (pow(mumuMass,2) > q2_bins[4]) && (pow(mumuMass,2) < q2_bins[5]) ){pass4->Fill(mumuMass);}
      else if( (q2Bin == 5) && (pow(mumuMass,2) > q2_bins[5]) && (pow(mumuMass,2) < q2_bins[6]) ){pass5->Fill(mumuMass);}
      else if( (q2Bin == 6) && (pow(mumuMass,2) > q2_bins[6]) && (pow(mumuMass,2) < q2_bins[7]) ){pass6->Fill(mumuMass);}
      else if( (q2Bin == 7) && (pow(mumuMass,2) > q2_bins[7]) && (pow(mumuMass,2) < q2_bins[8]) ){pass7->Fill(mumuMass);}
    }

    cout << "GEN MC - " << "q2Bin = " << q2Bin << " - year = " << year << endl;
    for(int evt = 0; evt < t_nocuts->GetEntries(); evt++){
      t_nocuts->GetEntry(evt);

      if( (q2Bin == 0) && (pow(gen_mumuMass,2) > q2_bins[0]) && (pow(gen_mumuMass,2) < q2_bins[1]) ){gen0->Fill(gen_mumuMass);}
      else if( (q2Bin == 1) && (pow(gen_mumuMass,2) > q2_bins[1]) && (pow(gen_mumuMass,2) < q2_bins[2]) ){gen1->Fill(gen_mumuMass);}
      else if( (q2Bin == 2) && (pow(gen_mumuMass,2) > q2_bins[2]) && (pow(gen_mumuMass,2) < q2_bins[3]) ){gen2->Fill(gen_mumuMass);}
      else if( (q2Bin == 3) && (pow(gen_mumuMass,2) > q2_bins[3]) && (pow(gen_mumuMass,2) < q2_bins[4]) ){gen3->Fill(gen_mumuMass);}
      else if( (q2Bin == 4) && (pow(gen_mumuMass,2) > q2_bins[4]) && (pow(gen_mumuMass,2) < q2_bins[5]) ){gen4->Fill(gen_mumuMass);}
      else if( (q2Bin == 5) && (pow(gen_mumuMass,2) > q2_bins[5]) && (pow(gen_mumuMass,2) < q2_bins[6]) ){gen5->Fill(gen_mumuMass);}
      else if( (q2Bin == 6) && (pow(gen_mumuMass,2) > q2_bins[6]) && (pow(gen_mumuMass,2) < q2_bins[7]) ){gen6->Fill(gen_mumuMass);}
      else if( (q2Bin == 7) && (pow(gen_mumuMass,2) > q2_bins[7]) && (pow(gen_mumuMass,2) < q2_bins[8]) ){gen7->Fill(gen_mumuMass);}
    }

    if(q2Bin == 0){eff0 = new TEfficiency(*pass0,*gen0);}
    if(q2Bin == 1){eff1 = new TEfficiency(*pass1,*gen1);}
    if(q2Bin == 2){eff2 = new TEfficiency(*pass2,*gen2);}
    if(q2Bin == 3){eff3 = new TEfficiency(*pass3,*gen3);}
    if(q2Bin == 4){eff4 = new TEfficiency(*pass4,*gen4);}
    if(q2Bin == 5){eff5 = new TEfficiency(*pass5,*gen5);}
    if(q2Bin == 6){eff6 = new TEfficiency(*pass6,*gen6);}
    if(q2Bin == 7){eff7 = new TEfficiency(*pass7,*gen7);}
  }

  eff_x_acc->SetBinContent(0,eff0->GetEfficiency(0));
  cout << "eff0 = " << eff0->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(1,eff1->GetEfficiency(0));
  cout << "eff1 = " << eff1->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(2,eff2->GetEfficiency(0));
  cout << "eff2 = " << eff2->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(3,eff3->GetEfficiency(0));
  cout << "eff3 = " << eff3->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(4,eff4->GetEfficiency(0));
  cout << "eff4 = " << eff4->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(5,eff5->GetEfficiency(0));
  cout << "eff5 = " << eff5->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(6,eff6->GetEfficiency(0));
  cout << "eff6 = " << eff6->GetEfficiency(0) << endl;
  eff_x_acc->SetBinContent(7,eff7->GetEfficiency(0));
  cout << "eff7 = " << eff7->GetEfficiency(0) << endl;

  TCanvas c;
  c.cd();
  eff_x_acc->SetTitle(Form("Efficiency x Acceptance - %i",year));
  eff_x_acc->GetXaxis()->SetTitle("q^2 [GeV]");
  eff_x_acc->GetYaxis()->SetTitle("#epsilon");
  eff_x_acc->SetMarkerStyle(7);
  eff_x_acc->SetMarkerColor(4);
  eff_x_acc->Draw("P");
  c.SaveAs(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.gif",year));

  TFile* f;
  f = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",year),"UPDATE");
  f->cd();
  eff_x_acc->Write();
  f->Close();
}




