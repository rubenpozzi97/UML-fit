#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <RooFitResult.h>
#include <TGraphErrors.h>
#include <TGraphPainter.h>
#include <TLine.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLatex.h>
#include <RooRealVar.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <vector>
#include <TVectorT.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TGraphPainter.h>
#include <TGraph.h>
#include <TAttFill.h>
#include <map>
#include <string>

using namespace std;

void fit_stability(int q2Bin){

  std::map<int,std::vector<float>> fM_sigmas = {
    {2016, {0.024, 0.017, 0.018, 0.013, 0.005, 0.010, 0.006, 0.010}},
    {2017, {0.020, 0.015, 0.016, 0.011, 0.004, 0.009, 0.005, 0.013}},
    {2018, {0.015, 0.011, 0.013, 0.008, 0.002, 0.006, 0.003, 0.008}},
  };

  int n_year = 3;

  TString names[] = {"mean_rt", "sigma_rt", "alpha_rt1", "n_rt1", "alpha_rt2", "n_rt2", "sigma_rt2", "f_rt", "mean_difference", "sigma_wt1", "alpha_wt1", "alpha_wt2", "n_wt1", "n_wt2", "sig_yield", "CB_yield", "mFrac"};

  TString names_SF[] = {"mean_rt", "sig_yield", "CB_yield"};
  
  TString names_RT[] = {"sigma_rt", "alpha_rt1", "n_rt1", "alpha_rt2", "n_rt2", "sigma_rt2", "f_rt", "mean_difference", "mFrac"};

  TString names_WT[] = {"sigma_wt1", "alpha_wt1", "alpha_wt2", "n_wt1", "n_wt2"};

  int n_names = sizeof(names)/sizeof(names[0]);

  for(int year = 6; year <= 8; year ++){

    TString variables[] = {Form("mean_{RT}^{201%i}",year), Form("#sigma_{RT1}^{201%i}",year), Form("#alpha_{RT1}^{201%i}",year), Form("n_{RT1}^{201%i}",year), Form("#alpha_{RT2}^{201%i}",year), Form("n_{RT2}^{201%i}",year),Form("#sigma_{RT2}^{201%i}",year), Form("f^{RT201%i}",year), Form("mean_difference^{201%i}",year), Form("#sigma_{WT1}^{201%i}",year), Form("#alpha_{WT1}^{201%i}",year), Form("#alpha_{WT2}^{201%i}",year), Form("n_{WT1}^{201%i}",year), Form("n_{WT2}^{201%i}",year), Form("sig_yield^{201%i}",year), Form("CB_yield^{201%i}",year), Form("mFrac^{201%i}",year)};

    TString variables_SF[] = {Form("mean_{RT}^{201%i}",year), Form("sig_yield^{201%i}",year), Form("CB_yield^{201%i}",year)};

    TString variables_RT[] = {Form("#sigma_{RT1}^{201%i}",year), Form("#alpha_{RT1}^{201%i}",year), Form("n_{RT1}^{201%i}",year), Form("#alpha_{RT2}^{201%i}",year), Form("n_{RT2}^{201%i}",year),Form("#sigma_{RT2}^{201%i}",year), Form("f^{RT201%i}",year), Form("mean_difference^{201%i}",year), Form("mFrac^{201%i}",year)};

    TString variables_WT[] = {Form("#sigma_{WT1}^{201%i}",year), Form("#alpha_{WT1}^{201%i}",year), Form("#alpha_{WT2}^{201%i}",year), Form("n_{WT1}^{201%i}",year), Form("n_{WT2}^{201%i}",year)};

    int n_var = sizeof(variables)/sizeof(variables[0]);
    int n_var_SF = sizeof(variables_SF)/sizeof(variables_SF[0]);
    int n_var_RT = sizeof(variables_RT)/sizeof(variables_RT[0]);

    TFile* f_nominal = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c1m0_subs0CT+WT.root",year,q2Bin));
    TFile* f_sf = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c2m0_subs0CT+WT.root",year,q2Bin)); // w/ scale factor
    TFile* f_range1 = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c1m1_subs0CT+WT.root",year,q2Bin)); // mass > 5.1
    TFile* f_range2 = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c1m2_subs0CT+WT.root",year,q2Bin)); // mass < 5.5
    TFile* f_pdf1 = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c1m3_subs0CT+WT.root",year,q2Bin)); // WT fixed / exp + erf

    RooFitResult* r_nominal = (RooFitResult*)f_nominal->Get(Form("simFitResult_b%ip2c1m0subs0",q2Bin));
    RooFitResult* r_range1 = (RooFitResult*)f_range1->Get(Form("simFitResult_b%ip2c1m1subs0",q2Bin));
    RooFitResult* r_range2 = (RooFitResult*)f_range2->Get(Form("simFitResult_b%ip2c1m2subs0",q2Bin));
    RooFitResult* r_pdf1 = (RooFitResult*)f_pdf1->Get(Form("simFitResult_b%ip2c1m3subs0",q2Bin));
    RooFitResult* r_sf = (RooFitResult*)f_sf->Get(Form("simFitResult_b%ip2c2m0subs0",q2Bin));

    TGraphErrors *gr[n_year][n_var];

    TLine* line;
    TLatex* latex0; // nominal
    TLatex* latex1; // range1
    TLatex* latex2; // range2
    TLatex* latex3; // pdf1
    TLatex* latex4; // scale factor

    double x_var[5];
    double xe_var[5];

    for(int i = 0; i < n_var; i++){
      if( (q2Bin <= 3) && ( (i == 6) || (i == 7) ) ){continue;}
      if( (q2Bin == 7) && ( (i == 4) || (i == 5) ) ){continue;}

      RooRealVar* var0;
      RooRealVar* var1;
      RooRealVar* var2;
      RooRealVar* var3;
      RooRealVar* var4;

      if(q2Bin != 4){ // for non-resonant channel and psi(2S), variations 0, 1, 2 have WT shape unfixed 
        var0 = (RooRealVar*)r_nominal->floatParsFinal().find(variables[i]);
        var1 = (RooRealVar*)r_range1->floatParsFinal().find(variables[i]);
        var2 = (RooRealVar*)r_range2->floatParsFinal().find(variables[i]);
      
        x_var[0] = var0->getVal();
        x_var[1] = var1->getVal();
        x_var[2] = var2->getVal();

        xe_var[0] = var0->getError();
        xe_var[1] = var1->getError();
        xe_var[2] = var2->getError();

        latex0 = new TLatex(x_var[0]-xe_var[0],1.1,Form("#scale[0.5]{Nominal: %f #pm %f }",x_var[0],xe_var[0]));
        latex1 = new TLatex(x_var[1]-xe_var[1],2.1,Form("#scale[0.5]{Exp(5.1<m<5.6): %f #pm %f}",x_var[1],xe_var[1]));
        latex2 = new TLatex(x_var[2]-xe_var[2],3.1,Form("#scale[0.5]{Exp(5.0<m<5.5): %f #pm %f}",x_var[2],xe_var[2]));
      }
      if(q2Bin == 6){ // for psi(2S) channel, variation 3 also has WT shape unfixed
        var3 = (RooRealVar*)r_pdf1->floatParsFinal().find(variables[i]);

        x_var[3] = var3->getVal();
        xe_var[3] = var3->getError();
   
        latex3 = new TLatex(x_var[3]-xe_var[3],4.1,Form("#scale[0.5]{Exp+erf: %f #pm %f}",x_var[3],xe_var[3]));
      }
      if(q2Bin == 4){ // for JPsi channel, variations 1, 2 and 3 have WT shape fixed
        var0 = (RooRealVar*)r_nominal->floatParsFinal().find(variables[i]);
        if( (std::find(std::begin(variables_RT), std::end(variables_RT), variables[i]) != std::end(variables_RT)) || (std::find(std::begin(variables_SF), std::end(variables_SF), variables[i]) != std::end(variables_SF)) ){
          var1 = (RooRealVar*)r_range1->floatParsFinal().find(variables[i]);
          var2 =  (RooRealVar*)r_range2->floatParsFinal().find(variables[i]);
          var3 = (RooRealVar*)r_pdf1->floatParsFinal().find(variables[i]);

          x_var[1] = var1->getVal();
          x_var[2] = var2->getVal();
          x_var[3] = var3->getVal();

          xe_var[1] = var1->getError();
          xe_var[2] = var2->getError();
          xe_var[3] = var3->getError();

          latex1 = new TLatex(x_var[1]-xe_var[1],2.1,Form("#scale[0.5]{Exp(5.1<m<5.6): %f #pm %f}",x_var[1],xe_var[1]));
          latex2 = new TLatex(x_var[2]-xe_var[2],3.1,Form("#scale[0.5]{Exp(5.0<m<5.5): %f #pm %f}",x_var[2],xe_var[2]));
          latex3 = new TLatex(x_var[3]-xe_var[3],4.1,Form("#scale[0.5]{Exp+erf: %f #pm %f}",x_var[3],xe_var[3]));
        }
        x_var[0] = var0->getVal(); // nominal has WT shape unfixed
        xe_var[0] = var0->getError();

        latex0 = new TLatex(x_var[0]-xe_var[0],1.1,Form("#scale[0.5]{Nominal: %f #pm %f}",x_var[0],xe_var[0]));
      }
      if((q2Bin != 4) && (q2Bin != 6)){ // for non-resonant channel, variation 3 has WT shape fixed
        if( (std::find(std::begin(variables_RT), std::end(variables_RT), variables[i]) != std::end(variables_RT)) || (std::find(std::begin(variables_SF), std::end(variables_SF), variables[i]) != std::end(variables_SF)) ){
          var3 = (RooRealVar*)r_pdf1->floatParsFinal().find(variables[i]);

          x_var[3] = var3->getVal();
          xe_var[3] = var3->getError();

          latex3 = new TLatex(x_var[3]-xe_var[3],4.1,Form("#scale[0.5]{WT fix: %f #pm %f}",x_var[3],xe_var[3]));
        }
      }
      if(std::find(std::begin(variables_SF), std::end(variables_SF), variables[i]) != std::end(variables_SF)){ // scale factor (variation 4)
        var4 = (RooRealVar*)r_sf->floatParsFinal().find(variables[i]); 

        x_var[4] = var4->getVal();
        xe_var[4] = var4->getError();
   
        latex4 = new TLatex(x_var[4]-xe_var[4],5.1,Form("#scale[0.5]{SF: %f #pm %f}",x_var[4],xe_var[4]));
      }

      TCanvas c = new TCanvas(variables[i]);
      TFile* f = new TFile("~/public/UML-fit/Fit_stability/root_files/"+names[i]+Form("_%i_b%i.root",year,q2Bin), "RECREATE");

      // signal_yield, CB_yield, mFrac
      if(std::find(std::begin(variables_SF), std::end(variables_SF), variables[i]) != std::end(variables_SF)){
        double x_sf[] = {x_var[0], x_var[1], x_var[2], x_var[3], x_var[4]};
        double xe_sf[] = {xe_var[0], xe_var[1], xe_var[2], xe_var[3], xe_var[4]};
        double y_sf[] = {1., 2., 3., 4., 5.};

        gr[year-6][i] = new TGraphErrors(5,x_sf,y_sf,xe_sf);
        gr[year-6][i]->GetYaxis()->SetRangeUser(0.5,5.5);

        line = new TLine(x_var[0],0.5,x_var[0],5.5);
        line->SetLineStyle(2);
        line->SetLineColor(kBlue);

        gr[year-6][i]->SetLineColor(kBlack);
        gr[year-6][i]->SetTitle( variables[i] + Form(" - q2Bin = %i ",q2Bin) );
        gr[year-6][i]->GetYaxis()->SetLabelSize(0);
        gr[year-6][i]->GetYaxis()->SetTickLength(0);

        c.cd();
        gr[year-6][i]->Draw("AP");
        line->Draw("same");
        latex0->Draw("same");
        latex1->Draw("same");
        latex2->Draw("same");
        latex3->Draw("same");
        latex4->Draw("same");

        f->cd();
        gr[year-6][i]->Write();
        f->Close();

        c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.pdf",year,q2Bin));
        c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.gif",year,q2Bin));
      }

      // WT variables
      if(std::find(std::begin(variables_WT), std::end(variables_WT), variables[i]) != std::end(variables_WT)){
        if(q2Bin == 4){
          double x_WT_JPsi[] = {x_var[0]};
          double xe_WT_JPsi[] = {xe_var[0]};
          double y_WT_JPsi[] = {1.};

          gr[year-6][i] = new TGraphErrors(1,x_WT_JPsi,y_WT_JPsi,xe_WT_JPsi);
          gr[year-6][i]->GetYaxis()->SetRangeUser(0.5,1.5);

          line = new TLine(x_var[0],0.5,x_var[0],1.5);
          line->SetLineStyle(2);
          line->SetLineColor(kBlue);

          gr[year-6][i]->SetLineColor(kBlack);
          gr[year-6][i]->SetTitle( variables[i] + Form(" - q2Bin = %i ",q2Bin) );
          gr[year-6][i]->GetYaxis()->SetLabelSize(0);
          gr[year-6][i]->GetYaxis()->SetTickLength(0);
  
          c.cd();
          gr[year-6][i]->Draw("AP");
          latex3->Draw("same");

          f->cd();
          gr[year-6][i]->Write();
          f->Close();

          c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.pdf",year,q2Bin));
          c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.gif",year,q2Bin));
        }
        else if((q2Bin != 4) && (q2Bin != 6)){
          double x_WT_non_res[] = {x_var[0], x_var[1], x_var[2]};
          double xe_WT_non_res[] = {xe_var[0], xe_var[1], xe_var[2]};
          double y_WT_non_res[] = {1., 2., 3.};

          gr[year-6][i] = new TGraphErrors(3,x_WT_non_res,y_WT_non_res,xe_WT_non_res);
          gr[year-6][i]->GetYaxis()->SetRangeUser(0.5,3.5);

          line = new TLine(x_var[0],0.5,x_var[0],3.5);
          line->SetLineStyle(2);
          line->SetLineColor(kBlue);

          gr[year-6][i]->SetLineColor(kBlack);
          gr[year-6][i]->SetTitle( variables[i] + Form(" - q2Bin = %i ",q2Bin) );
          gr[year-6][i]->GetYaxis()->SetLabelSize(0);
          gr[year-6][i]->GetYaxis()->SetTickLength(0);

          c.cd();
          gr[year-6][i]->Draw("AP");
          line->Draw("same");
          latex0->Draw("same");
          latex1->Draw("same");
          latex2->Draw("same");          
          latex3->Draw("same");

          f->cd();
          gr[year-6][i]->Write();
          f->Close();          

          c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.pdf",year,q2Bin));
          c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.gif",year,q2Bin));
        }
        else if(q2Bin == 6){
          double x_WT_Psi[] = {x_var[0], x_var[1], x_var[2], x_var[3]};
          double xe_WT_Psi[] = {xe_var[0], xe_var[1], xe_var[2], xe_var[3]};
          double y_WT_Psi[] = {1., 2., 3., 4.};

          gr[year-6][i] = new TGraphErrors(4,x_WT_Psi,y_WT_Psi,xe_WT_Psi);
          gr[year-6][i]->GetYaxis()->SetRangeUser(0.5,4.5);

          line = new TLine(x_var[0],0.5,x_var[0],4.5);
          line->SetLineStyle(2);
          line->SetLineColor(kBlue);

          gr[year-6][i]->SetLineColor(kBlack);
          gr[year-6][i]->SetTitle( variables[i] + Form(" - q2Bin = %i ",q2Bin) );
          gr[year-6][i]->GetYaxis()->SetLabelSize(0);
          gr[year-6][i]->GetYaxis()->SetTickLength(0);

          c.cd();
          gr[year-6][i]->Draw("AP");
          line->Draw("same");
          latex0->Draw("same");
          latex1->Draw("same");
          latex2->Draw("same");
          latex3->Draw("same");
          latex4->Draw("same");

          f->cd();
          gr[year-6][i]->Write();
          f->Close();

          c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.pdf",year,q2Bin));
          c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.gif",year,q2Bin));
        }
      }

      // RT variables + mean_difference
      if(std::find(std::begin(variables_RT), std::end(variables_RT), variables[i]) != std::end(variables_RT)){
        double x_RT[] = {x_var[0], x_var[1], x_var[2], x_var[3]};
        double xe_RT[] = {xe_var[0], xe_var[1], xe_var[2], xe_var[3]};
        double y_RT[] = {1., 2., 3., 4.};
  
        gr[year-6][i] = new TGraphErrors(4,x_RT,y_RT,xe_RT);
        gr[year-6][i]->GetYaxis()->SetRangeUser(0.5,4.5);

        line = new TLine(x_var[0],0.5,x_var[0],4.5);
        line->SetLineStyle(2);
        line->SetLineColor(kBlue);

        gr[year-6][i]->SetLineColor(kBlack);
        gr[year-6][i]->SetTitle( variables[i] + Form(" - q2Bin = %i ",q2Bin) );
        gr[year-6][i]->GetYaxis()->SetLabelSize(0);
        gr[year-6][i]->GetYaxis()->SetTickLength(0);

        gr[year-6][i]->Draw("AP");
        line->Draw("same");
        latex0->Draw("same");
        latex1->Draw("same");
        latex2->Draw("same");
        latex3->Draw("same");
        latex4->Draw("same");

        f->cd();
        gr[year-6][i]->Write();
        f->Close();

        c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.pdf",year,q2Bin));
        c.SaveAs("~/public/UML-fit/Fit_stability/single_year/"+variables[i]+Form("_%i_b%i.gif",year,q2Bin));
      }
      c.Clear();
    }//loop over variables
  }//loop over years  

  for(int i = 0; i < n_names; i++){

    if( (q2Bin <= 3) && ( (i == 6) || (i == 7) ) ){continue;}
    if( (q2Bin == 7) && ( (i == 4) || (i == 5) ) ){continue;}

    TMultiGraph* mg = new TMultiGraph();
    
    TFile* f_6 = new TFile("~/public/UML-fit/Fit_stability/root_files/"+names[i]+Form("_%i_b%i.root",6,q2Bin));
    TFile* f_7 = new TFile("~/public/UML-fit/Fit_stability/root_files/"+names[i]+Form("_%i_b%i.root",7,q2Bin));
    TFile* f_8 = new TFile("~/public/UML-fit/Fit_stability/root_files/"+names[i]+Form("_%i_b%i.root",8,q2Bin));
 
    TFile* f_RT6 = new TFile(Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0CT.root", 2016, q2Bin));
    TFile* f_RT7 = new TFile(Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0CT.root", 2017, q2Bin));
    TFile* f_RT8 = new TFile(Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0CT.root", 2018, q2Bin));

    RooFitResult* r_RT6 = (RooFitResult*)f_RT6->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));
    RooFitResult* r_RT7 = (RooFitResult*)f_RT7->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));
    RooFitResult* r_RT8 = (RooFitResult*)f_RT8->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));

    TFile* f_WT6 = new TFile(Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0WT.root", 2016, q2Bin));
    TFile* f_WT7 = new TFile(Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0WT.root", 2017, q2Bin));
    TFile* f_WT8 = new TFile(Form("simFitMassResults/simFitResult_recoMC_fullMass%i_MCStat_b%ip2c0m0_subs0WT.root", 2018, q2Bin));

    RooFitResult* r_WT6 = (RooFitResult*)f_WT6->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));
    RooFitResult* r_WT7 = (RooFitResult*)f_WT7->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));
    RooFitResult* r_WT8 = (RooFitResult*)f_WT8->Get(Form("simFitResult_b%ip2c0m0subs0",q2Bin));

    TGraphErrors* g6 = (TGraphErrors*)f_6->Get("Graph");
    TGraphErrors* g7 = (TGraphErrors*)f_7->Get("Graph");
    TGraphErrors* g8 = (TGraphErrors*)f_8->Get("Graph");

    g6->SetLineColor(kBlack);
    g7->SetLineColor(kRed);
    g8->SetLineColor(kBlue);

    mg->Add(g6);
    mg->Add(g7);
    mg->Add(g8);

    TCanvas c1;
    c1.cd();
    mg->Draw("AP");

    mg->SetTitle( names[i] + Form(" - q2Bin = %i ",q2Bin) );
    mg->GetYaxis()->SetLabelSize(0);
    mg->GetYaxis()->SetTickLength(0);
    mg->GetYaxis()->SetRangeUser(0.5,5.5);

    TLatex* tex0;
    TLatex* tex1;
    TLatex* tex2;
    TLatex* tex3;
    TLatex* tex4;

    if(q2Bin != 4){
      tex0 = new TLatex(mg->GetXaxis()->GetXmin(), 1.1,"#scale[0.5]{Nominal}"); 
      tex1 = new TLatex(mg->GetXaxis()->GetXmin(), 2.1,"#scale[0.5]{Exp(5.1<mass<5.6)}");
      tex2 = new TLatex(mg->GetXaxis()->GetXmin(), 3.1,"#scale[0.5]{Exp(5.0<mass<5.5)}");
    }
    if(q2Bin == 6){
      tex3 = new TLatex(mg->GetXaxis()->GetXmin(), 4.1,"#scale[0.5]{Exp+erf}");
    }
    if(q2Bin == 4){
      if( (std::find(std::begin(names_RT), std::end(names_RT), names[i]) != std::end(names_RT)) || (std::find(std::begin(names_SF), std::end(names_SF), names[i]) != std::end(names_SF)) ){
        tex1 = new TLatex(mg->GetXaxis()->GetXmin(), 2.1,"#scale[0.5]{Exp(5.1<mass<5.6)}");
        tex2 = new TLatex(mg->GetXaxis()->GetXmin(), 3.1,"#scale[0.5]{Exp(5.0<mass<5.5)}");
        tex3 = new TLatex(mg->GetXaxis()->GetXmin(), 4.1,"#scale[0.5]{Exp+erf}");
      }
      tex0 = new TLatex(mg->GetXaxis()->GetXmin(), 1.1,"#scale[0.5]{Nominal}");
    }
    if((q2Bin != 4) && (q2Bin != 6)){
      if( (std::find(std::begin(names_RT), std::end(names_RT), names[i]) != std::end(names_RT)) || (std::find(std::begin(names_SF), std::end(names_SF), names[i]) != std::end(names_SF)) ){
        tex3 = new TLatex(mg->GetXaxis()->GetXmin(), 4.1,"#scale[0.5]{WT fixed}");
      }
    }
    if(std::find(std::begin(names_SF), std::end(names_SF), names[i]) != std::end(names_SF)){
        tex4 = new TLatex(mg->GetXaxis()->GetXmin(), 5.1,"#scale[0.5]{SF}");
    }

    TLegend* leg = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg->SetFillColor(0);
    leg->AddEntry(g6, "2016", "lp");
    leg->AddEntry(g7, "2017", "lp");    
    leg->AddEntry(g8, "2018", "lp");                       

    // RT vars
    TString var_RT6[] = {Form("mean_{RT}^{201%i}",6), Form("#sigma_{RT1}^{201%i}",6), Form("#alpha_{RT1}^{201%i}",6), Form("n_{RT1}^{201%i}",6), Form("#alpha_{RT2}^{201%i}",6), Form("n_{RT2}^{201%i}",6), Form("#sigma_{RT2}^{201%i}",6), Form("f^{RT201%i}",6)};
    TString var_RT7[] = {Form("mean_{RT}^{201%i}",7), Form("#sigma_{RT1}^{201%i}",7), Form("#alpha_{RT1}^{201%i}",7), Form("n_{RT1}^{201%i}",7), Form("#alpha_{RT2}^{201%i}",7), Form("n_{RT2}^{201%i}",7), Form("#sigma_{RT2}^{201%i}",7), Form("f^{RT201%i}",7)};
    TString var_RT8[] = {Form("mean_{RT}^{201%i}",8), Form("#sigma_{RT1}^{201%i}",8), Form("#alpha_{RT1}^{201%i}",8), Form("n_{RT1}^{201%i}",8), Form("#alpha_{RT2}^{201%i}",8), Form("n_{RT2}^{201%i}",8), Form("#sigma_{RT2}^{201%i}",8), Form("f^{RT201%i}",8)};

    // WT vars
    TString var_WT6[] = {Form("#sigma_{WT1}^{201%i}",6), Form("#alpha_{WT1}^{201%i}",6), Form("#alpha_{WT2}^{201%i}",6), Form("n_{WT1}^{201%i}",6), Form("n_{WT2}^{201%i}",6), Form("sig_yield^{201%i}",6), Form("CB_yield^{201%i}",6), Form("mFrac^{201%i}",6)};
    TString var_WT7[] = {Form("#sigma_{WT1}^{201%i}",7), Form("#alpha_{WT1}^{201%i}",7), Form("#alpha_{WT2}^{201%i}",7), Form("n_{WT1}^{201%i}",7), Form("n_{WT2}^{201%i}",7), Form("sig_yield^{201%i}",7), Form("CB_yield^{201%i}",7)};
    TString var_WT8[] = {Form("#sigma_{WT1}^{201%i}",8), Form("#alpha_{WT1}^{201%i}",8), Form("#alpha_{WT2}^{201%i}",8), Form("n_{WT1}^{201%i}",8), Form("n_{WT2}^{201%i}",8), Form("sig_yield^{201%i}",8), Form("CB_yield^{201%i}",8)};    

    RooRealVar* var6;
    RooRealVar* var7;
    RooRealVar* var8;

    TGraph* gr6;
    TGraph* gr7;
    TGraph* gr8;

    if(i == 8){//mean_difference

      double y1[5] = {0.5, 0.5, 4.5, 4.5, 0.5};

      RooRealVar* var_RT6 = (RooRealVar*)r_RT6->floatParsFinal().find(Form("mean_{RT}^{201%i}",6));
      RooRealVar* var_RT7 = (RooRealVar*)r_RT7->floatParsFinal().find(Form("mean_{RT}^{201%i}",7));
      RooRealVar* var_RT8 = (RooRealVar*)r_RT8->floatParsFinal().find(Form("mean_{RT}^{201%i}",8));

      RooRealVar* var_WT6 = (RooRealVar*)r_WT6->floatParsFinal().find(Form("mean_{WT}^{201%i}",6));
      RooRealVar* var_WT7 = (RooRealVar*)r_WT7->floatParsFinal().find(Form("mean_{WT}^{201%i}",7));
      RooRealVar* var_WT8 = (RooRealVar*)r_WT8->floatParsFinal().find(Form("mean_{WT}^{201%i}",8));

      double m6 = (var_WT6->getVal()) - (var_RT6->getVal());
      double m7 = (var_WT7->getVal()) - (var_RT7->getVal());
      double m8 = (var_WT8->getVal()) - (var_RT8->getVal());

      double me6 = sqrt( pow((var_RT6->getError()),2) + pow((var_WT6->getError()),2) );
      double me7 = sqrt( pow((var_RT7->getError()),2) + pow((var_WT7->getError()),2) );
      double me8 = sqrt( pow((var_RT8->getError()),2) + pow((var_WT8->getError()),2) );

      double x1_6[5] = {m6-me6, m6+me6, m6+me6, m6-me6, m6-me6};
      double x1_7[5] = {m7-me7, m7+me7, m7+me7, m7-me7, m7-me7};
      double x1_8[5] = {m8-me8, m8+me8, m8+me8, m8-me8, m8-me8};

      gr6 = new TGraph(5,x1_6, y1);
      gr6->SetLineColor(kBlack);
      gr6->SetFillColorAlpha(kBlack,0.25);

      gr7 = new TGraph(5,x1_7, y1);
      gr7->SetLineColor(kRed);
      gr7->SetFillColorAlpha(kRed,0.25);

      gr8 = new TGraph(5,x1_8, y1);
      gr8->SetLineColor(kBlue);
      gr8->SetFillColorAlpha(kBlue,0.25);

      if( x1_6[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_6[0],mg->GetXaxis()->GetXmax());}
      else if( x1_7[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_7[0],mg->GetXaxis()->GetXmax());}
      else if( x1_8[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_8[0],mg->GetXaxis()->GetXmax());}
      else if( mg->GetXaxis()->GetXmax() < x1_6[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_6[1]);}
      else if( mg->GetXaxis()->GetXmax() < x1_7[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_7[1]);}
      else if( mg->GetXaxis()->GetXmax() < x1_8[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_8[1]);}

    }
    else if( (i<=13) && (i != 8) ){// RT variables
 
      if(i <= 7){
        var6 = (RooRealVar*)r_RT6->floatParsFinal().find(var_RT6[i]);
        var7 = (RooRealVar*)r_RT7->floatParsFinal().find(var_RT7[i]);
        var8 = (RooRealVar*)r_RT8->floatParsFinal().find(var_RT8[i]);
      }
      else{
        var6 = (RooRealVar*)r_WT6->floatParsFinal().find(var_WT6[i-9]);
        var7 = (RooRealVar*)r_WT7->floatParsFinal().find(var_WT7[i-9]);
        var8 = (RooRealVar*)r_WT8->floatParsFinal().find(var_WT8[i-9]);
      }

      double y1[5] = {0.5, 0.5, 4.5, 4.5, 0.5};
      double x1_6[5] = {var6->getVal()-var6->getError(), var6->getVal()+var6->getError(), var6->getVal()+var6->getError(), var6->getVal()-var6->getError(), var6->getVal()-var6->getError()};
      double x1_7[5] = {var7->getVal()-var7->getError(), var7->getVal()+var7->getError(), var7->getVal()+var7->getError(), var7->getVal()-var7->getError(), var7->getVal()-var7->getError()};
      double x1_8[5] = {var8->getVal()-var8->getError(), var8->getVal()+var8->getError(), var8->getVal()+var8->getError(), var8->getVal()-var8->getError(), var8->getVal()-var8->getError()};

      gr6 = new TGraph(5,x1_6, y1);
      gr6->SetLineColor(kBlack);
      gr6->SetFillColorAlpha(kBlack,0.25);

      gr7 = new TGraph(5,x1_7, y1);
      gr7->SetLineColor(kRed);
      gr7->SetFillColorAlpha(kRed,0.25);

      gr8 = new TGraph(5,x1_8, y1);
      gr8->SetLineColor(kBlue);
      gr8->SetFillColorAlpha(kBlue,0.25);

      if( x1_6[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_6[0],mg->GetXaxis()->GetXmax());}
      else if( x1_7[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_7[0],mg->GetXaxis()->GetXmax());}
      else if( x1_8[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_8[0],mg->GetXaxis()->GetXmax());}
      else if( mg->GetXaxis()->GetXmax() < x1_6[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_6[1]);}
      else if( mg->GetXaxis()->GetXmax() < x1_7[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_7[1]);}
      else if( mg->GetXaxis()->GetXmax() < x1_8[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_8[1]);}

    }
    else if(i == 16){//mFrac
      RooRealVar* RT_yield_2016 = (RooRealVar*)r_RT6->floatParsFinal().find(Form("sig_yield^{201%i}",6));
      RooRealVar* RT_yield_2017 = (RooRealVar*)r_RT7->floatParsFinal().find(Form("sig_yield^{201%i}",7));
      RooRealVar* RT_yield_2018 = (RooRealVar*)r_RT8->floatParsFinal().find(Form("sig_yield^{201%i}",8));

      RooRealVar* WT_yield_2016 = (RooRealVar*)r_WT6->floatParsFinal().find(Form("sig_yield^{201%i}",6));
      RooRealVar* WT_yield_2017 = (RooRealVar*)r_WT7->floatParsFinal().find(Form("sig_yield^{201%i}",7));
      RooRealVar* WT_yield_2018 = (RooRealVar*)r_WT8->floatParsFinal().find(Form("sig_yield^{201%i}",8));

      double n_rt6 = RT_yield_2016->getVal();
      double n_rt7 = RT_yield_2017->getVal();
      double n_rt8 = RT_yield_2018->getVal();
  
      double n_wt6 = WT_yield_2016->getVal();
      double n_wt7 = WT_yield_2017->getVal();
      double n_wt8 = WT_yield_2018->getVal();
 
      double fraction6 = n_wt6/(n_rt6+n_wt6);
      double fraction7 = n_wt7/(n_rt7+n_wt7);
      double fraction8 = n_wt8/(n_rt8+n_wt8);

      double width6 = fM_sigmas[2016][q2Bin];
      double width7 = fM_sigmas[2017][q2Bin]; 
      double width8 = fM_sigmas[2018][q2Bin]; 

      double y1[5] = {0.5, 0.5, 4.5, 4.5, 0.5};
      double x1_6[5] = {fraction6-width6, fraction6+width6, fraction6+width6, fraction6-width6, fraction6-width6};
      double x1_7[5] = {fraction7-width7, fraction7+width7, fraction7+width7, fraction7-width7, fraction7-width7};
      double x1_8[5] = {fraction8-width8, fraction8+width8, fraction8+width8, fraction8-width8, fraction8-width8};

      gr6 = new TGraph(5,x1_6, y1);
      gr6->SetLineColor(kBlack);
      gr6->SetFillColorAlpha(kBlack,0.25);

      gr7 = new TGraph(5,x1_7, y1);
      gr7->SetLineColor(kRed);
      gr7->SetFillColorAlpha(kRed,0.25);

      gr8 = new TGraph(5,x1_8, y1);
      gr8->SetLineColor(kBlue);
      gr8->SetFillColorAlpha(kBlue,0.25);

      if( x1_6[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_6[0],mg->GetXaxis()->GetXmax());}
      else if( x1_7[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_7[0],mg->GetXaxis()->GetXmax());}
      else if( x1_8[0] < mg->GetXaxis()->GetXmin() ){mg->GetXaxis()->SetLimits(x1_8[0],mg->GetXaxis()->GetXmax());}
      else if( mg->GetXaxis()->GetXmax() < x1_6[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_6[1]);}
      else if( mg->GetXaxis()->GetXmax() < x1_7[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_7[1]);}
      else if( mg->GetXaxis()->GetXmax() < x1_8[1] ){mg->GetXaxis()->SetLimits(mg->GetXaxis()->GetXmin(),x1_8[1]);}

    }

    leg->Draw("same");
    if( ((i<=13) && (i!=0)) || (i==16)){
      gr6->Draw("F");
      gr7->Draw("F");
      gr8->Draw("F");
    }
 
    tex0->Draw("same");
    tex1->Draw("same");
    tex2->Draw("same");
    tex3->Draw("same");

    TString names_sf[] = {"mean_rt", "sig_yield", "CB_yield"};
    if(std::find(std::begin(names_sf), std::end(names_sf), names[i]) != std::end(names_sf)){
      tex4->Draw("same");
    }

    c1.SaveAs("~/public/UML-fit/Fit_stability/all_years/"+names[i]+Form("_b%i.pdf",q2Bin));
    c1.SaveAs("~/public/UML-fit/Fit_stability/all_years/"+names[i]+Form("_b%i.gif",q2Bin));
  }
}


