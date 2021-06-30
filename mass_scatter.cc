#include <TFile.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <TH1.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <RooRealVar.h>
#include <RooAbsRealLValue.h>
#include <RooGlobalFunc.h>

#include <THistPainter.h>
#include "TStyle.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "TH2.h"
#include "TPad.h"
#include <TTree.h>
#include <TH2F.h>
#include <TLine.h>

using namespace RooFit;

void mass_scatter(int year){

  TFile* f = new TFile(Form("/eos/cms/store/user/fiorendi/p5prime/%i/skims/newphi/%iData_All_finalSelection.root",year,year));
  TTree* t = (TTree*)f->Get("ntuple");

  double PDGB0Mass = 5.27958;
  double PDGJpsiMass = 3.096916;
  double PDGPsiPrimeMass = 3.686109;
  double PDGKstMass = 0.896;

  // mass variables
  double mass, mumuMass;
  t->SetBranchAddress("bMass",&mass);
  t->SetBranchAddress("mumuMass",&mumuMass);

  // cut variables
  bool passB0Psi_lmnr, passB0Psi_jpsi, passB0Psi_psip;
  t->SetBranchAddress( "passB0Psi_lmnr", &passB0Psi_lmnr );
  t->SetBranchAddress( "passB0Psi_jpsi", &passB0Psi_jpsi );
  t->SetBranchAddress( "passB0Psi_psip", &passB0Psi_psip );

  // cut to remove B+->Psi(2S)K->Jpsi pi pi K
  double wt_mass, wt_kstarmass, kaonPt, pionPt, mmpiMass, mmkMass;
  t->SetBranchAddress( "wt_mass",      &wt_mass      );
  t->SetBranchAddress( "wt_kstarmass", &wt_kstarmass );
  t->SetBranchAddress( "kaonPt",       &kaonPt       );
  t->SetBranchAddress( "pionPt",       &pionPt       );
  t->SetBranchAddress( "mmpiMass",     &mmpiMass     );
  t->SetBranchAddress( "mmkMass",      &mmkMass      );

  double x0Cut=-0.4;
  double y0Cut= 0.3;
  double x1Cut= 0.6;
  double y1Cut=-0.1;

  double x_0Cut=3;
  double y_0Cut=3.8;
  double x_1Cut=3.6;
  double y_1Cut=4.8;

  double CutX1=3.2;
  double CutX2=3.6;
  double CutY1=4.7;
  double CutY2=4.9;

  TH2F* histo = new TH2F("histo","histo",200,5.0,5.6,200,1.0,4.5);

  TH2F* histo_before = new TH2F("histo_before","histo_before",200,5.0,5.6,200,1.0,16);
  TH2F* histo_after = new TH2F("histo_after","histo_after",200,5.0,5.6,200,1.0,16);
  TH2F* histo_radiation = new TH2F("histo_radiation","histo_radiation",200,5.0,5.6,200,1.0,16);

  for(int i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    histo_before->Fill(mass,pow(mumuMass,2));

    if((pow(mumuMass,2) > 8.68) && (pow(mumuMass,2) < 10.09) && (passB0Psi_jpsi == 0)){continue;}
    else if((pow(mumuMass,2) > 12.86) && (pow(mumuMass,2) < 14.18)  && (passB0Psi_psip == 0)){continue;}
    else if(passB0Psi_lmnr == 0){continue;}
    histo_radiation->Fill(mass,pow(mumuMass,2));

    bool XCut= (( (PDGB0Mass - wt_mass) - y0Cut ) / (y1Cut-y0Cut)) < (((wt_kstarmass-PDGKstMass)-x0Cut) / (x1Cut-x0Cut)) && \
                  kaonPt > pionPt && \
                  (wt_kstarmass-PDGKstMass)>0 && \
                  (mmpiMass > CutX1 && mmpiMass < CutX2) && \
                  (mmkMass >  CutY1 && mmkMass  < CutY2) && \
                  ((mmkMass - y_0Cut) / (y_1Cut - y_0Cut)) > ((mmpiMass-x_0Cut)/(x_1Cut-x_0Cut));

    if (XCut && ((pow(mumuMass,2) > 8.68) && (pow(mumuMass,2) < 10.09)) ) continue;

    histo->Fill(mass,mumuMass);
    histo_after->Fill(mass,pow(mumuMass,2));
  }

  TCanvas c;
  c.cd();

  gStyle->SetOptStat(0);
  histo->SetTitle("");
  histo->GetXaxis()->SetTitle("m(K^{+}#pi^{-}#mu^{+}#mu^{-}) [GeV/c^{2}]");
  histo->GetYaxis()->SetTitle("m(#mu^{+}#mu^{-}) [GeV/c^{2}]");
  histo->Draw("COLZ");

  c.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_%i.gif",year));

  TLine* line0 = new TLine(5.0, 1.0, 5.6, 1.0);
  TLine* line1 = new TLine(5.0, 2.0, 5.6, 2.0);
  TLine* line2 = new TLine(5.0, 4.3, 5.6, 4.3);
  TLine* line3 = new TLine(5.0, 6.0, 5.6, 6.0);
  TLine* line4 = new TLine(5.0, 8.68, 5.6, 8.68);
  TLine* line5 = new TLine(5.0, 10.09, 5.6, 10.09);
  TLine* line6 = new TLine(5.0, 12.86, 5.6, 12.86);
  TLine* line7 = new TLine(5.0, 14.18, 5.6, 14.18);
  TLine* line8 = new TLine(5.0, 16.0, 5.6, 16.0);

  line0->SetLineStyle(1);
  line1->SetLineStyle(1);
  line2->SetLineStyle(1);
  line3->SetLineStyle(1);
  line4->SetLineStyle(1);
  line5->SetLineStyle(1);
  line6->SetLineStyle(1);
  line7->SetLineStyle(1);
  line8->SetLineStyle(1);

  line0->SetLineWidth(3);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line3->SetLineWidth(3);
  line4->SetLineWidth(2);
  line5->SetLineWidth(2);
  line6->SetLineWidth(2);
  line7->SetLineWidth(2);
  line8->SetLineWidth(3);

  line0->SetLineColor(kBlack);
  line1->SetLineColor(kBlack);
  line2->SetLineColor(kBlack);
  line3->SetLineColor(kBlack);
  line4->SetLineColor(kRed);
  line5->SetLineColor(kRed);
  line6->SetLineColor(kGreen);
  line7->SetLineColor(kGreen);
  line8->SetLineColor(kBlack);

  TCanvas c1;
  c1.cd();

  gStyle->SetOptStat(0);
  histo_before->SetTitle("");
  histo_before->GetXaxis()->SetTitle("m(K^{+}#pi^{-}#mu^{+}#mu^{-}) [GeV/c^{2}]");
  histo_before->GetYaxis()->SetTitle("m(#mu^{+}#mu^{-})^{2} [GeV^{2}/c^{4}]");
  histo_before->Draw("COLZ");  
  line0->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");

  c1.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_q2bins_before_%i.gif",year));
  c1.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_q2bins_before_%i.pdf",year));

  TCanvas c2;
  c2.cd();
  
  gStyle->SetOptStat(0);
  histo_after->SetTitle("");
  histo_after->GetXaxis()->SetTitle("m(K^{+}#pi^{-}#mu^{+}#mu^{-}) [GeV/c^{2}]");
  histo_after->GetYaxis()->SetTitle("m(#mu^{+}#mu^{-})^{2} [GeV^{2}/c^{4}]");
  histo_after->Draw("COLZ");
  line0->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");

  c2.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_q2bins_after_%i.gif",year));
  c2.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_q2bins_after_%i.pdf",year));

  TCanvas c3;
  c3.cd();

  gStyle->SetOptStat(0);
  histo_radiation->SetTitle("");
  histo_radiation->GetXaxis()->SetTitle("m(K^{+}#pi^{-}#mu^{+}#mu^{-}) [GeV/c^{2}]");
  histo_radiation->GetYaxis()->SetTitle("m(#mu^{+}#mu^{-})^{2} [GeV^{2}/c^{4}]");
  histo_radiation->Draw("COLZ");
  line0->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");

  c3.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_q2bins_radiation_%i.gif",year));
  c3.SaveAs(Form("~/public/UML-fit/Mass_scatter/mass_scatter_q2bins_radiation_%i.pdf",year));

  return;
}




