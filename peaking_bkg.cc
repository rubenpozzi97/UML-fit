#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooGenericPdf.h>
#include <RooArgList.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include "RooFit.h"
#include <RooPlot.h>
#include <RooHist.h>
#include <TPad.h>
#include <TLine.h>
#include <RooProduct.h>
#include <RooCBShape.h>
#include <RooExponential.h>
#include <RooGaussian.h>
#include <RooArgSet.h>
#include <TMath.h>
#include <RooDataHist.h>
#include <TAxis.h>

using namespace RooFit;

int genBroPdgId [4];
double genBroPt  [4];
double genBroEta [4];
double genBroPhi [4];

int matchTrack (double trkPt, double trkEta, double trkPhi, bool isPos);

void peaking_bkg(){

  double kMass = 0.493677;
  double piMass = 0.13957039;
  double muMass = 0.10566;

  auto fin = TFile::Open("/eos/cms/store/user/fiorendi/p5prime/sidebands/2018MC_HBJPSIX_newphi_punzi_noTkMu_B0PsiFlag_addVars.root");
  auto tin = (TTree*)fin->Get("ntuple");

  int genJpsiAncestorPdgId;

  double genJpsiMupPt, genJpsiMupEta, genJpsiMupPhi;
  double genJpsiMumPt, genJpsiMumEta, genJpsiMumPhi;

  double kstTrkpPt, kstTrkpEta, kstTrkpPhi;
  double kstTrkmPt, kstTrkmEta, kstTrkmPhi;
  
  double tagged_mass;

  double genWeight;

  for (int i=0; i<4; ++i) {
    tin->SetBranchAddress(Form("genBro%iPdgId",i+1),&genBroPdgId[i]);
    tin->SetBranchAddress(Form("genBro%iPt",i+1),&genBroPt[i]);
    tin->SetBranchAddress(Form("genBro%iEta",i+1),&genBroEta[i]);
    tin->SetBranchAddress(Form("genBro%iPhi",i+1),&genBroPhi[i]);
  }

  tin->SetBranchAddress("genJpsiAncestorPdgId",&genJpsiAncestorPdgId);

  tin->SetBranchAddress("tagged_mass",&tagged_mass);

  tin->SetBranchAddress("genWeight",&genWeight);

  tin->SetBranchAddress("kstTrkpPt",&kstTrkpPt);
  tin->SetBranchAddress("kstTrkpEta",&kstTrkpEta);
  tin->SetBranchAddress("kstTrkpPhi",&kstTrkpPhi);
  tin->SetBranchAddress("kstTrkmPt",&kstTrkmPt);
  tin->SetBranchAddress("kstTrkmEta",&kstTrkmEta);
  tin->SetBranchAddress("kstTrkmPhi",&kstTrkmPhi);

  tin->SetBranchAddress("genJpsiMupPt",&genJpsiMupPt);
  tin->SetBranchAddress("genJpsiMupEta",&genJpsiMupEta);
  tin->SetBranchAddress("genJpsiMupPhi",&genJpsiMupPhi);
  tin->SetBranchAddress("genJpsiMumPt",&genJpsiMumPt);
  tin->SetBranchAddress("genJpsiMumEta",&genJpsiMumEta);
  tin->SetBranchAddress("genJpsiMumPhi",&genJpsiMumPhi);

  auto hBMass = new TH1D("hBMass","hBMass",100,4.9,5.7);
  auto hBMassPeak = new TH1D("hBMassPeak","hBMassPeak",100,5.0,5.6);
  auto hBMassPart = new TH1D("hBMassPart","hBMassPart",100,5.0,5.6);
  auto hBMassComb = new TH1D("hBMassComb","hBMassComb",100,4.9,5.7);
  auto hBMassDoubComb = new TH1D("hBMassDoubComb","hBMassDoubComb",100,4.9,5.7);
  auto hGenB0Mass = new TH1D("hGenB0Mass","hGenB0Mass",100,4.9,5.7);
  auto hGenBpMass = new TH1D("hGenBpMass","hGenBpMass",100,4.9,5.7);
  auto hGenBsMass = new TH1D("hGenBsMass","hGenBsMass",100,4.9,5.7);
  auto hGenBcMass = new TH1D("hGenBcMass","hGenBcMass",100,4.9,5.7);

  double sig_yield;
  double bkg_yield;

  for (int i=0; i<tin->GetEntries(); ++i) {
    tin->GetEntry(i);

    int matchP = matchTrack(kstTrkpPt,kstTrkpEta,kstTrkpPhi,true);
    int matchM = matchTrack(kstTrkmPt,kstTrkmEta,kstTrkmPhi,false);

    if ( matchP<0 && matchM<0 ) {
      hBMassDoubComb->Fill(tagged_mass,genWeight);
      continue;
    }
    
    if ( matchP<0 || matchM<0 ) {
      hBMassComb->Fill(tagged_mass,genWeight);
      continue;
    }

    bool isPartReco = false;
    for (int j=0; j<4; ++j) {
      
      if ( genBroPdgId[j] == 0 ) continue;
      if ( fabs(genBroPdgId[j]) == 11 ) isPartReco = true;

      if ( j == matchP || j == matchM ) continue;

      isPartReco = true;

    }

    if (!isPartReco) {

      double trkpMass = piMass;
      double trkmMass = piMass;
      if (fabs(genBroPdgId[matchP])==321) trkpMass = kMass;
      if (fabs(genBroPdgId[matchM])==321) trkmMass = kMass;
    
/*    ROOT::Math::PtEtaPhiMVector genTrkpVec (genBroPt[matchP],genBroEta[matchP],genBroPhi[matchP],trkpMass);
      ROOT::Math::PtEtaPhiMVector genTrkmVec (genBroPt[matchM],genBroEta[matchM],genBroPhi[matchM],trkmMass);

      ROOT::Math::PtEtaPhiMVector genMupVec (genJpsiMupPt,genJpsiMupEta,genJpsiMupPhi,muMass);
      ROOT::Math::PtEtaPhiMVector genMumVec (genJpsiMumPt,genJpsiMumEta,genJpsiMumPhi,muMass);

      if ( fabs(genJpsiAncestorPdgId) == 511 || fabs(genJpsiAncestorPdgId) == 513 ) {
	if ( (genTrkpVec+genTrkmVec+genMupVec+genMumVec).mass() < 5.16 ) isPartReco = true;
	hGenB0Mass->Fill( (genTrkpVec+genTrkmVec+genMupVec+genMumVec).mass(), genWeight );
      // } else if ( fabs(genJpsiAncestorPdgId) == 521 || fabs(genJpsiAncestorPdgId) == 523 ) {
      //   hGenBpMass->Fill( (genTrkpVec+genTrkmVec+genMupVec+genMumVec).mass(), genWeight );
      } else if ( fabs(genJpsiAncestorPdgId) == 531 || fabs(genJpsiAncestorPdgId) == 533 ) {
	if ( (genTrkpVec+genTrkmVec+genMupVec+genMumVec).mass() < 5.34 ) isPartReco = true;
	hGenBsMass->Fill( (genTrkpVec+genTrkmVec+genMupVec+genMumVec).mass(), genWeight );
      // } else if ( fabs(genJpsiAncestorPdgId) == 541 || fabs(genJpsiAncestorPdgId) == 543 ) {
      //   hGenBcMass->Fill( (genTrkpVec+genTrkmVec+genMupVec+genMumVec).mass(), genWeight );
      } // else cout<<genJpsiAncestorPdgId<<endl;
      
      // hGenJpsiMass->Fill( (genMupVec+genMumVec).mass(), genWeight );
      // hGenKstMass->Fill( (genBro1Vec+genBro2Vec).mass(), genWeight );
*/
    }

    if (isPartReco) {
      hBMassPart->Fill(tagged_mass,genWeight);
      if((tagged_mass > 5.0) && (tagged_mass < 5.6)){bkg_yield += 1;}
      continue;
    }

    if ( ( fabs(genJpsiAncestorPdgId) == 511 || fabs(genJpsiAncestorPdgId) == 513 ) &&
	 ( ( genBroPdgId[matchP] == 321 && genBroPdgId[matchM] == -211 ) ||
	   ( genBroPdgId[matchP] == 211 && genBroPdgId[matchM] == -321 ) ) ){
      hBMass->Fill(tagged_mass,genWeight);
      if((tagged_mass > 5.0) && (tagged_mass < 5.6)){sig_yield += 1;}
    }
    else{ 
      hBMassPeak->Fill(tagged_mass,genWeight);
    }
    
  }

  hBMassPeak->SetLineColor(8);
  hBMassPart->SetLineColor(2);
  hBMassComb->SetLineColor(7);
  hBMassDoubComb->SetLineColor(38);
  // hGenB0Mass->SetLineColor(3);
  // hGenBpMass->SetLineColor(3);
  // hGenBsMass->SetLineColor(6);
  // hGenBcMass->SetLineColor(6);

  auto canv = new TCanvas();
  canv->cd();
  hBMass->Draw();
  // hGenB0Mass->Draw();
  // hGenBpMass->Draw("same");
  // hGenBsMass->Draw("same");
  // hGenBcMass->Draw("same");
  hBMassPeak->Draw("same");
  hBMassPart->Draw("same");
  hBMassComb->Draw("same");
  hBMassDoubComb->Draw("same");
  canv->SaveAs("~/public/UML-fit/Peaking_bkg/bkg_components.gif");

  RooRealVar* mass = new RooRealVar("mass", "mass", 5.0, 5.6);

  RooDataHist* data_hist = new RooDataHist("data_hist", "data_hist", *mass, hBMassPart);


  cout << "sig yield = " << sig_yield << endl;
  cout << "bkg yield = " << bkg_yield << endl;

  double fraction = bkg_yield/sig_yield;
  cout << "fraction = " << fraction << endl;

  RooRealVar* frac = new RooRealVar("frac","frac",fraction,0.,1.);
  frac->setConstant();
  cout << "frac = " << frac->getVal() << endl;

  // PDF model
  RooRealVar* scale = new RooRealVar("scale", "scale", 12., 0., 30.);
  RooRealVar* shift = new RooRealVar("shift", "shift", 5.13, 4.0, 5.6);
  RooRealVar* yield = new RooRealVar("yield", "yield", data_hist->numEntries(), 0, data_hist->numEntries()*2);

  RooGenericPdf* erf_pdf = new RooGenericPdf("erf_pdf","erf_pdf","TMath::Erfc((mass-shift)*scale)",RooArgList(*mass,*shift,*scale));

  RooRealVar* lambda = new RooRealVar("lambda", "lambda", -1., -10., 5.);

  RooExponential* exp_pdf = new RooExponential("exp_pdf","exp_pdf",*mass,*lambda);

  RooRealVar* f_erf = new RooRealVar("f_erf","f_erf",0.5,0.,1.);

  RooAddPdf* final_pdf = new RooAddPdf("final_pdf", "final_pdf", RooArgList(*erf_pdf,*exp_pdf), *f_erf);

  TFile* fout = new TFile("~/public/UML-fit/Peaking_bkg/peaking_bkg.root", "RECREATE");
  fout->cd();
  RooFitResult* fitresult =  final_pdf->fitTo(*data_hist,Save());
  fitresult->Print("v");
  fitresult->Write();
  frac->Write();
  fout->Close();

  RooPlot* massframe = mass->frame();
  data_hist->plotOn(massframe,Name("Data"), MarkerColor(kBlack));
  final_pdf->plotOn(massframe, Name("Fit"), LineColor(kRed));
  final_pdf->paramOn(massframe,Layout(0.6,0.9,0.9));

  TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
  leg->AddEntry(massframe->findObject("Data"),"Data","lep");
  leg->AddEntry(massframe->findObject("Fit"),"Fit","lep");

  TCanvas c;
  c.cd();

  TPad *p1 = new TPad("p1","p1",0.0,0.27,0.99,0.99);
  p1->SetTitle("Part. reco. background");
  p1->SetBorderMode(1);
  p1->SetFrameBorderMode(0);
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.10);
  p1->Draw();

  TPad *p2 = new TPad("p2","p2",0.0,0.065,0.99,0.24);
  p2->SetTitle("");
  p2->SetTopMargin(0.);
  p2->SetBottomMargin(0.2);
  p2->SetBorderMode(1);
  p2->Draw();

  p1->cd();
  massframe->GetYaxis()->SetTitle("events / GeV");
  massframe->GetXaxis()->SetTitle("mass [GeV]");
  massframe->Draw();
  //leg->Draw("same");

  RooHist* pull_hist = massframe->pullHist("Data","Fit");
  RooPlot *pull_plot = mass->frame();

  pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"P");
  pull_plot->SetTitle("");

  pull_plot->GetXaxis()->SetTitle("");
  pull_plot->GetXaxis()->SetTitleFont(42);
  pull_plot->GetXaxis()->SetTitleSize(0.17);
  pull_plot->GetXaxis()->SetTitleOffset(1.09);
  pull_plot->GetXaxis()->SetLabelFont(42);
  pull_plot->GetXaxis()->SetLabelSize(0.15);
  pull_plot->GetXaxis()->SetLabelOffset(0.01);
  pull_plot->GetXaxis()->SetTickLength(0.13);

  pull_hist->GetYaxis()->SetTitle("Pull hist");
  pull_hist->GetYaxis()->SetTitleFont(42);
  pull_plot->GetYaxis()->SetTitleSize(0.10);
  pull_plot->GetYaxis()->SetTitleOffset(1.09);
  pull_plot->GetYaxis()->SetLabelFont(42);
  pull_plot->GetYaxis()->SetLabelSize(0.13);
  pull_plot->GetYaxis()->SetLabelOffset(0.005);
  pull_plot->GetYaxis()->SetNdivisions(305);

  gPad->Update();
  TLine *line = new TLine(gPad->GetUxmin(), 0, gPad->GetUxmax(), 0);
  line->SetLineStyle(2);
  line->SetLineColor(kBlue);

  p2->cd();
  pull_plot->Draw();
  line->Draw("same");

  c.SaveAs("~/public/UML-fit/Peaking_bkg/peaking_bkg_fit.gif");
  c.SaveAs("~/public/UML-fit/Peaking_bkg/peaking_bkg_fit.pdf");
  
}

int matchTrack (double trkPt, double trkEta, double trkPhi, bool isPos){

  for (int i=0; i<4; ++i) {

    if ( genBroPdgId[i]==0 ) continue;
    if ( isPos ) {
      if ( genBroPdgId[i]!=321 && genBroPdgId[i]!=211 && genBroPdgId[i]!=-11 ) continue;
    } else {
      if ( genBroPdgId[i]!=-321 && genBroPdgId[i]!=-211 && genBroPdgId[i]!=11 ) continue;
    }

    double deltaPhi = fabs( trkPhi - genBroPhi[i] );
    if (deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;

    if ( deltaPhi*deltaPhi + (trkEta-genBroEta[i])*(trkEta-genBroEta[i]) > 0.0001 ) continue;

    return i;

  }

  return -1;
  
}

