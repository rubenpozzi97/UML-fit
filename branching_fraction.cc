#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <RooFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphPainter.h>
#include <TLine.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TLatex.h>

void branching_fraction(){

  const int n_q2Bin = 8;
  double q2Bin[n_q2Bin+1] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};

  //bin width
  double bin_width[n_q2Bin];
  for(int i = 0; i < n_q2Bin; i++){
    bin_width[i] = q2Bin[i+1]-q2Bin[i];
  }

  //efficiencies + efficiency systematics
  TFile* f_eff_2016 = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",2016));
  TFile* f_eff_2017 = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",2017));
  TFile* f_eff_2018 = new TFile(Form("~/public/UML-fit/Efficiency/eff_x_acc_%i.root",2018));

  TGraphErrors* eff_2016 = (TGraphErrors*)f_eff_2016->Get("Graph");
  TGraphErrors* eff_2017 = (TGraphErrors*)f_eff_2017->Get("Graph");
  TGraphErrors* eff_2018 = (TGraphErrors*)f_eff_2018->Get("Graph");

  double *eff_x_acc_2016 = eff_2016->GetY();
  double *eff_x_acc_2017 = eff_2017->GetY();
  double *eff_x_acc_2018 = eff_2018->GetY();

  TGraph* efficiency_2016 = (TGraph*)f_eff_2016->Get("Graph;1");
  TGraph* efficiency_2017 = (TGraph*)f_eff_2017->Get("Graph;1");
  TGraph* efficiency_2018 = (TGraph*)f_eff_2018->Get("Graph;1");

  TGraph* acceptance_2016 = (TGraph*)f_eff_2016->Get("Graph;2");
  TGraph* acceptance_2017 = (TGraph*)f_eff_2017->Get("Graph;2");
  TGraph* acceptance_2018 = (TGraph*)f_eff_2018->Get("Graph;2");

  TGraph* weighted_2016 = (TGraph*)f_eff_2016->Get("Graph;3");
  TGraph* weighted_2017 = (TGraph*)f_eff_2017->Get("Graph;3");
  TGraph* weighted_2018 = (TGraph*)f_eff_2018->Get("Graph;3");

  double* acc_2016 = acceptance_2016->GetY();
  double* acc_2017 = acceptance_2017->GetY();
  double* acc_2018 = acceptance_2018->GetY();

  double* wei_2016 = weighted_2016->GetY();
  double* wei_2017 = weighted_2017->GetY();
  double* wei_2018 = weighted_2018->GetY();

  //double *eff_syst_2016 = eff_2016->GetEY();
  //double *eff_syst_2017 = eff_2017->GetEY();
  //double *eff_syst_2018 = eff_2018->GetEY();

  //yields + yield systematics + yied statistical errors
  double eff_x_acc_wei_2016[n_q2Bin];
  double eff_x_acc_wei_2017[n_q2Bin];
  double eff_x_acc_wei_2018[n_q2Bin];

  double yields_2016[n_q2Bin];
  double yields_2017[n_q2Bin];
  double yields_2018[n_q2Bin];

  double yields_errors_2016[n_q2Bin];
  double yields_errors_2017[n_q2Bin];
  double yields_errors_2018[n_q2Bin];

  double yields_average[n_q2Bin];
  double yields_average_stat[n_q2Bin];
  double yields_average_syst[n_q2Bin];

  TFile* f_yield_syst_2016 = new TFile("~/public/UML-fit/Systematics/root_files/yield_syst_2016.root");
  TFile* f_yield_syst_2017 = new TFile("~/public/UML-fit/Systematics/root_files/yield_syst_2017.root");
  TFile* f_yield_syst_2018 = new TFile("~/public/UML-fit/Systematics/root_files/yield_syst_2018.root");

  TGraphErrors* graph_2016 = (TGraphErrors*)f_yield_syst_2016->Get("Graph;2");
  TGraphErrors* graph_2017 = (TGraphErrors*)f_yield_syst_2017->Get("Graph;2");
  TGraphErrors* graph_2018 = (TGraphErrors*)f_yield_syst_2018->Get("Graph;2");

  double *yield_syst_2016 = graph_2016->GetEY();
  double *yield_syst_2017 = graph_2017->GetEY();
  double *yield_syst_2018 = graph_2018->GetEY();

  for(int i = 0; i < n_q2Bin; i++){
    eff_x_acc_wei_2016[i] = wei_2016[i]*acc_2016[i];
    eff_x_acc_wei_2017[i] = wei_2017[i]*acc_2017[i];
    eff_x_acc_wei_2018[i] = wei_2018[i]*acc_2018[i];

    TFile* f_yield_2016 = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass%i_DATA_b%ip2c1m0_subs0CT+WT.root",2016,i));
    TFile* f_yield_2017 = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass%i_DATA_b%ip2c1m0_subs0CT+WT.root",2017,i));
    TFile* f_yield_2018 = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass%i_DATA_b%ip2c1m0_subs0CT+WT.root",2018,i));

    RooFitResult* fitresult_2016 = (RooFitResult*)f_yield_2016->Get(Form("simFitResult_b%ip2c1m0subs0",i));
    RooFitResult* fitresult_2017 = (RooFitResult*)f_yield_2017->Get(Form("simFitResult_b%ip2c1m0subs0",i));
    RooFitResult* fitresult_2018 = (RooFitResult*)f_yield_2018->Get(Form("simFitResult_b%ip2c1m0subs0",i));

    RooRealVar* sig_yield_2016 = (RooRealVar*)fitresult_2016->floatParsFinal().find(Form("sig_yield^{%i}",2016));   
    RooRealVar* sig_yield_2017 = (RooRealVar*)fitresult_2017->floatParsFinal().find(Form("sig_yield^{%i}",2017));
    RooRealVar* sig_yield_2018 = (RooRealVar*)fitresult_2018->floatParsFinal().find(Form("sig_yield^{%i}",2018));

    yields_2016[i] = sig_yield_2016->getVal();
    yields_2017[i] = sig_yield_2017->getVal();
    yields_2018[i] = sig_yield_2018->getVal();

    yields_errors_2016[i] = (sig_yield_2016->getError());
    yields_errors_2017[i] = (sig_yield_2017->getError());
    yields_errors_2018[i] = (sig_yield_2018->getError());

    yields_average[i] = (1./3.)*(yields_2016[i]+yields_2017[i]+yields_2018[i]);
    yields_average_stat[i] = (1./3.)*sqrt( pow(yields_errors_2016[i],2) + pow(yields_errors_2017[i],2) + pow(yields_errors_2018[i],2) );
    yields_average_syst[i] = (1./3.)*sqrt( pow(yield_syst_2016[i],2) + pow(yield_syst_2017[i],2) + pow(yield_syst_2018[i],2) );

    delete f_yield_2016;
    delete f_yield_2017;
    delete f_yield_2018;
  }

  //resonant channels branching ratio
  double PDG_JPsi_K = 0.00127;
  double PDG_JPsi_K_error = 0.00005;

  double PDG_Psi_K = 0.00059;
  double PDG_Psi_K_error = 0.00004;

  double PDG_JPsi = 0.05961;
  double PDG_JPsi_error = 0.00033;

  double PDG_Psi = 0.0080;
  double PDG_Psi_error = 0.0006;

  double PDG_ratio = PDG_Psi_K/PDG_JPsi_K;
  double PDG_ratio_error = sqrt( pow((1/PDG_JPsi_K),2)*pow(PDG_Psi_K_error,2) + pow((PDG_Psi_K/(PDG_JPsi_K*PDG_JPsi_K)),2)*pow(PDG_JPsi_K_error,2) );

  double ratio_2016 =  (yields_2016[6]/yields_2016[4])*(eff_x_acc_2016[4]/eff_x_acc_2016[6])*(PDG_JPsi/PDG_Psi);
  double ratio_2017 =  (yields_2017[6]/yields_2017[4])*(eff_x_acc_2017[4]/eff_x_acc_2017[6])*(PDG_JPsi/PDG_Psi);
  double ratio_2018 =  (yields_2018[6]/yields_2018[4])*(eff_x_acc_2018[4]/eff_x_acc_2018[6])*(PDG_JPsi/PDG_Psi);

  double ratio_wei_2016 = (yields_2016[6]/yields_2016[4])*(eff_x_acc_wei_2016[4]/eff_x_acc_wei_2016[6])*(PDG_JPsi/PDG_Psi);
  double ratio_wei_2017 = (yields_2017[6]/yields_2017[4])*(eff_x_acc_wei_2017[4]/eff_x_acc_wei_2017[6])*(PDG_JPsi/PDG_Psi);
  double ratio_wei_2018 = (yields_2018[6]/yields_2018[4])*(eff_x_acc_wei_2018[4]/eff_x_acc_wei_2018[6])*(PDG_JPsi/PDG_Psi);

  double eff_syst1_2016 = abs(ratio_2016-ratio_wei_2016);
  double eff_syst1_2017 = abs(ratio_2017-ratio_wei_2017);
  double eff_syst1_2018 = abs(ratio_2018-ratio_wei_2018);

  double ratio_error_2016_syst = sqrt( pow(yield_syst_2016[6]/yields_2016[6],2) + pow(yield_syst_2016[4]/yields_2016[4],2) + pow(eff_syst1_2016/ratio_2016,2) + pow(PDG_JPsi_error/PDG_JPsi,2) + pow(PDG_Psi_error/PDG_Psi,2) )*ratio_2016;
  double ratio_error_2017_syst = sqrt( pow(yield_syst_2017[6]/yields_2017[6],2) + pow(yield_syst_2017[4]/yields_2017[4],2) + pow(eff_syst1_2017/ratio_wei_2017,2) + pow(PDG_JPsi_error/PDG_JPsi,2) + pow(PDG_Psi_error/PDG_Psi,2) )*ratio_2017;
  double ratio_error_2018_syst = sqrt( pow(yield_syst_2018[6]/yields_2018[6],2) + pow(yield_syst_2018[4]/yields_2018[4],2) + pow(eff_syst1_2018/ratio_wei_2018,2) + pow(PDG_JPsi_error/PDG_JPsi,2) + pow(PDG_Psi_error/PDG_Psi,2) )*ratio_2018;

  double ratio_error_2016_stat = sqrt( pow(yields_errors_2016[4]/yields_2016[4],2) + pow(yields_errors_2016[6]/yields_2016[6],2) )*ratio_2016;
  double ratio_error_2017_stat = sqrt( pow(yields_errors_2017[4]/yields_2017[4],2) + pow(yields_errors_2017[6]/yields_2017[6],2) )*ratio_2017;
  double ratio_error_2018_stat = sqrt( pow(yields_errors_2018[4]/yields_2018[4],2) + pow(yields_errors_2018[6]/yields_2018[6],2) )*ratio_2018;

  double ratio = (1./3.)*(ratio_2016+ratio_2017+ratio_2018);
  double ratio_error_syst =  (1./3.)*sqrt( pow(ratio_error_2016_syst,2) + pow(ratio_error_2017_syst,2) + pow(ratio_error_2018_syst,2) );
  double ratio_error_stat = (1./3.)*sqrt( pow(ratio_error_2016_stat,2) + pow(ratio_error_2017_stat,2) + pow(ratio_error_2018_stat,2) );

  double x[] = {PDG_ratio, ratio_2016, ratio_2017, ratio_2018, ratio};
  double y[] = {1., 2., 3., 4., 5.};
  double ex_stat[] = {0., ratio_error_2016_stat, ratio_error_2017_stat, ratio_error_2018_stat, ratio_error_stat};
  double ex_syst[] = {PDG_ratio_error, ratio_error_2016_syst, ratio_error_2017_syst, ratio_error_2018_syst, ratio_error_syst}; 

  cout << '|' << setw(15) << "B(JPsi)" << '|' << setw(15) << "Yield Syst" << '|' << setw(15) << "Eff syst" << '|' << setw(15) << "Total syst" << '|' << setw(15) << "Total stat " << '|' << endl;

  cout << '|' << setw(15) << PDG_JPsi_error/PDG_JPsi << '|' << setw(15) << yield_syst_2016[4]/yields_2016[4] << '|' << setw(15) << eff_syst1_2016/ratio_2016 << '|' << setw(15) << ratio_error_2016_syst/ratio_2016 << '|' << setw(15) << ratio_error_2016_stat/ratio_2016   << '|' << endl;
  cout << '|' << setw(15) << PDG_JPsi_error/PDG_JPsi << '|' << setw(15) << yield_syst_2017[4]/yields_2017[4] << '|' << setw(15) << eff_syst1_2017/ratio_2017 << '|' << setw(15) << ratio_error_2017_syst/ratio_2017 << '|' << setw(15) << ratio_error_2017_stat/ratio_2017   << '|' << endl;
  cout << '|' << setw(15) << PDG_JPsi_error/PDG_JPsi << '|' << setw(15) << yield_syst_2018[4]/yields_2018[4] << '|' << setw(15) << eff_syst1_2018/ratio_2018 << '|' << setw(15) << ratio_error_2018_syst/ratio_2018 << '|' << setw(15) << ratio_error_2018_stat/ratio_2018   << '|' << endl;

  TGraphErrors* gr_stat = new TGraphErrors(5,x,y,ex_stat);
  gr_stat->SetLineColor(kBlack);

  TGraphErrors* gr_syst = new TGraphErrors(5,x,y,ex_syst);
  gr_syst->SetLineColor(kRed);

  TMultiGraph* mg = new TMultiGraph();

  mg->Add(gr_syst);
  mg->Add(gr_stat);

  mg->SetTitle("");
  mg->GetYaxis()->SetLabelSize(0);
  mg->GetYaxis()->SetTickLength(0);
  mg->GetYaxis()->SetRangeUser(0.5,5.5);

  TLine* line = new TLine(PDG_ratio,0.5,PDG_ratio,5.5);
  line->SetLineStyle(2);
  line->SetLineColor(kBlue);
 
  TLatex* latex1 = new TLatex(PDG_ratio-PDG_ratio_error,1.1,Form("#scale[0.5]{PDG: %f #pm %f }",PDG_ratio,PDG_ratio_error));
  TLatex* latex2 = new TLatex(ratio_2016-ratio_error_2016_syst,2.1,Form("#scale[0.5]{2016: %f #pm %f #pm %f}",ratio_2016,ratio_error_2016_stat,ratio_error_2016_syst));
  TLatex* latex3 = new TLatex(ratio_2017-ratio_error_2017_syst,3.1,Form("#scale[0.5]{2017: %f #pm %f #pm %f}",ratio_2017,ratio_error_2017_stat,ratio_error_2017_syst));
  TLatex* latex4 = new TLatex(ratio_2018-ratio_error_2018_syst,4.1,Form("#scale[0.5]{2018: %f #pm %f #pm %f}",ratio_2018,ratio_error_2018_stat,ratio_error_2018_syst));
  TLatex* latex5 = new TLatex(ratio-ratio_error_syst,5.1,Form("#scale[0.5]{Year average: %f #pm %f #pm %f}",ratio,ratio_error_stat,ratio_error_syst));

  TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(gr_stat, "Statistical Uncertainty", "lp");
  leg->AddEntry(gr_syst, "Systematic Uncertainty", "lp");

  TCanvas graph;
  graph.cd();

  mg->Draw("AP");
  line->Draw("same");
  latex1->Draw("same");
  latex2->Draw("same");
  latex3->Draw("same");
  latex4->Draw("same");
  latex5->Draw("same");
  leg->Draw("same");

  graph.SaveAs("~/public/UML-fit/Branching_fraction/resonance_ratio.gif");
  graph.SaveAs("~/public/UML-fit/Branching_fraction/resonance_ratio.pdf");

  // differential measurement
  double dB_dq_2016[n_q2Bin];
  double dB_dq_2016_wei[n_q2Bin];
  double dB_dq_syst_2016[n_q2Bin];
  double dB_dq_stat_2016[n_q2Bin];

  double dB_dq_2017[n_q2Bin];
  double dB_dq_2017_wei[n_q2Bin];
  double dB_dq_syst_2017[n_q2Bin];
  double dB_dq_stat_2017[n_q2Bin];

  double dB_dq_2018[n_q2Bin];
  double dB_dq_2018_wei[n_q2Bin];
  double dB_dq_syst_2018[n_q2Bin];
  double dB_dq_stat_2018[n_q2Bin];

  double eff_syst_2016[n_q2Bin];
  double eff_syst_2017[n_q2Bin];
  double eff_syst_2018[n_q2Bin];

  for(int i = 0; i < n_q2Bin; i++){
    if((i != 4) && (i != 6)){  
      dB_dq_2016[i] = (yields_2016[i]/yields_2016[4])*(eff_x_acc_2016[4]/eff_x_acc_2016[i])*(PDG_JPsi_K*PDG_JPsi/bin_width[i]);
      dB_dq_2017[i] = (yields_2017[i]/yields_2017[4])*(eff_x_acc_2017[4]/eff_x_acc_2017[i])*(PDG_JPsi_K*PDG_JPsi/bin_width[i]);
      dB_dq_2018[i] = (yields_2018[i]/yields_2018[4])*(eff_x_acc_2018[4]/eff_x_acc_2018[i])*(PDG_JPsi_K*PDG_JPsi/bin_width[i]);

      dB_dq_2016_wei[i] = (yields_2016[i]/yields_2016[4])*(eff_x_acc_wei_2016[4]/eff_x_acc_wei_2016[i])*(PDG_JPsi_K*PDG_JPsi/bin_width[i]);    
      dB_dq_2017_wei[i] = (yields_2017[i]/yields_2017[4])*(eff_x_acc_wei_2017[4]/eff_x_acc_wei_2017[i])*(PDG_JPsi_K*PDG_JPsi/bin_width[i]);
      dB_dq_2018_wei[i] = (yields_2018[i]/yields_2018[4])*(eff_x_acc_wei_2018[4]/eff_x_acc_wei_2018[i])*(PDG_JPsi_K*PDG_JPsi/bin_width[i]);     

      eff_syst_2016[i] = abs(dB_dq_2016[i]-dB_dq_2016_wei[i]);
      eff_syst_2017[i] = abs(dB_dq_2017[i]-dB_dq_2017_wei[i]);
      eff_syst_2018[i] = abs(dB_dq_2018[i]-dB_dq_2018_wei[i]);

      dB_dq_syst_2016[i] = sqrt( pow(yield_syst_2016[i]/yields_2016[i],2) + pow(yield_syst_2016[4]/yields_2016[4],2) + pow(eff_syst_2016[i]/dB_dq_2016[i],2) + pow(eff_syst_2016[4]/eff_x_acc_2016[4],2) + pow(PDG_JPsi_error/PDG_JPsi,2) + pow(PDG_JPsi_K_error/PDG_JPsi_K,2) )*dB_dq_2016[i];
      dB_dq_syst_2017[i] = sqrt( pow(yield_syst_2017[i]/yields_2017[i],2) + pow(yield_syst_2017[4]/yields_2017[4],2) + pow(eff_syst_2017[i]/dB_dq_2017[i],2) + pow(eff_syst_2017[4]/eff_x_acc_2017[4],2) + pow(PDG_JPsi_error/PDG_JPsi,2) + pow(PDG_JPsi_K_error/PDG_JPsi_K,2) )*dB_dq_2017[i];
      dB_dq_syst_2018[i] = sqrt( pow(yield_syst_2018[i]/yields_2018[i],2) + pow(yield_syst_2018[4]/yields_2018[4],2) + pow(eff_syst_2018[i]/dB_dq_2018[i],2) + pow(eff_syst_2018[4]/eff_x_acc_2018[4],2) + pow(PDG_JPsi_error/PDG_JPsi,2) + pow(PDG_JPsi_K_error/PDG_JPsi_K,2) )*dB_dq_2018[i];

      dB_dq_stat_2016[i] = sqrt( pow(yields_errors_2016[i]/yields_2016[i],2) + pow(yields_errors_2016[4]/yields_2016[4],2) )*dB_dq_2016[i];   
      dB_dq_stat_2017[i] = sqrt( pow(yields_errors_2017[i]/yields_2017[i],2) + pow(yields_errors_2017[4]/yields_2017[4],2) )*dB_dq_2017[i];
      dB_dq_stat_2018[i] = sqrt( pow(yields_errors_2018[i]/yields_2018[i],2) + pow(yields_errors_2018[4]/yields_2018[4],2) )*dB_dq_2018[i];

      cout << '|' << setw(15) << "dB/dq^{2}" << '|' << setw(15) << "Yield syst" << '|' << setw(15) << "Eff syst" << '|' << setw(15) << "B_JPsi_K syst" << '|' << setw(15) << "B_JPsi syst " << '|' << setw(15) << "Total syst" << '|' << setw(15) << "Total stat" << '|' << endl;
      cout << '|' << setw(15) << dB_dq_2016[i] << '|' << setw(15) << yield_syst_2016[i]/yields_2016[i]  << '|' << setw(15) << eff_syst_2016[i]/dB_dq_2016[i] << '|' << setw(15) << PDG_JPsi_K_error/PDG_JPsi_K << '|' << setw(15) << PDG_JPsi_error/PDG_JPsi << '|' << setw(15) << dB_dq_syst_2016[i]/dB_dq_2016[i] << '|' << setw(15) << dB_dq_stat_2016[i]/dB_dq_2016[i] << '|' << endl;
      cout << '|' << setw(15) << dB_dq_2017[i] << '|' << setw(15) << yield_syst_2017[i]/yields_2017[i]  << '|' << setw(15) << eff_syst_2017[i]/dB_dq_2017[i] << '|' << setw(15) << PDG_JPsi_K_error/PDG_JPsi_K << '|' << setw(15) << PDG_JPsi_error/PDG_JPsi << '|' << setw(15) << dB_dq_syst_2017[i]/dB_dq_2017[i] << '|' << setw(15) << dB_dq_stat_2017[i]/dB_dq_2017[i] << '|' << endl;
      cout << '|' << setw(15) << dB_dq_2018[i] << '|' << setw(15) << yield_syst_2018[i]/yields_2018[i]  << '|' << setw(15) << eff_syst_2018[i]/dB_dq_2018[i] << '|' << setw(15) << PDG_JPsi_K_error/PDG_JPsi_K << '|' << setw(15) << PDG_JPsi_error/PDG_JPsi << '|' << setw(15) << dB_dq_syst_2018[i]/dB_dq_2018[i] << '|' << setw(15) << dB_dq_stat_2018[i]/dB_dq_2018[i] << '|' << endl;
    }
  }

  double q2Bins_half[n_q2Bin];
  double q2Bins_err[n_q2Bin];

  for(int i = 0; i < n_q2Bin; i++){
    if((i != 4) && (i != 6)){
      q2Bins_half[i] = q2Bin[i] + 0.5* (q2Bin[i+1]-q2Bin[i]);
      q2Bins_err[i] = 0.5* (q2Bin[i+1]-q2Bin[i]);
    }
    else{
      q2Bins_half[i] = 0;
      q2Bins_err[i] = 0;
    }
  }

  double JPsi_mean = (q2Bin[5]-q2Bin[4])/2 + q2Bin[4];
  double Psi_mean = (q2Bin[7]-q2Bin[6])/2 + q2Bin[6];

  double JPsi[5] = {JPsi_mean-bin_width[4]/2, JPsi_mean+bin_width[4]/2, JPsi_mean+bin_width[4]/2, JPsi_mean-bin_width[4]/2, JPsi_mean-bin_width[4]/2};
  double Psi[5] = {Psi_mean-bin_width[6]/2, Psi_mean+bin_width[6]/2, Psi_mean+bin_width[6]/2, Psi_mean-bin_width[6]/2, Psi_mean-bin_width[6]/2};

  // 8 TeV CMS results
  double dB_dq_8TeV[] = {4.6*pow(10,-8), 3.3*pow(10,-8), 3.4*pow(10,-8), 4.7*pow(10,-8), 0., 6.2*pow(10,-8), 0., 6.7*pow(10,-8)};
  double dB_dq_8TeV_stat[] = {0.7*pow(10,-8), 0.5*pow(10,-8), 0.5*pow(10,-8), 0.4*pow(10,-8), 0., 0.4*pow(10,-8), 0., 0.6*pow(10,-8)}; 
  double dB_dq_8TeV_syst[] = {0.3*pow(10,-8), 0.2*pow(10,-8), 0.3*pow(10,-8), 0.3*pow(10,-8), 0., 0.5*pow(10,-8), 0., 0.5*pow(10,-8)};

  TGraphErrors* gr_stat_8TeV = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_8TeV,q2Bins_err,dB_dq_8TeV_stat);
  TGraphErrors* gr_syst_8TeV = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_8TeV,q2Bins_err,dB_dq_8TeV_syst);
  gr_stat_8TeV->SetLineColor(kGreen);
  gr_syst_8TeV->SetLineColor(kBlue);
  gr_stat_8TeV->SetMarkerColor(kGreen);
  gr_syst_8TeV->SetMarkerColor(kBlue);

  // LHCb results
  double lhcb_bins[] = {0.1, 2.0, 4.3, 8.68, 10.09, 12.86, 14.18, 16., 19.};
  double n_lhcb_bins = 8.;  

  double q2Bins_half_lhcb[n_q2Bin];
  double q2Bins_err_lhcb[n_q2Bin];

  for(int i = 0; i < n_lhcb_bins; i++){
    if((i != 3) && (i != 5)){
      q2Bins_half_lhcb[i] = lhcb_bins[i] + 0.5* (lhcb_bins[i+1]-lhcb_bins[i]);
      q2Bins_err_lhcb[i] = 0.5* (lhcb_bins[i+1]-lhcb_bins[i]);
    }
    else{
      q2Bins_half_lhcb[i] = 0;
      q2Bins_err_lhcb[i] = 0;
    }
  }

  double dB_dq_lhcb[] = {0.6*pow(10,-7), 0.3*pow(10,-7), 0.49*pow(10,-7), 0.43*pow(10,-7), 0.56*pow(10,-7), 0.41*pow(10,-7)};
  double dB_dq_lhcb_stat[] = {0.06*pow(10,-7), 0.03*pow(10,-7), 0.04*pow(10,-7), 0.04*pow(10,-7), 0.06*pow(10,-7), 0.04*pow(10,-7)};
  double dB_dq_lhcb_syst[] = {0.05*pow(10,-7), 0.03*pow(10,-7), 0.04*pow(10,-7), 0.04*pow(10,-7), 0.04*pow(10,-7), 0.04*pow(10,-7)};

  TGraphErrors* gr_stat_lhcb = new TGraphErrors(n_lhcb_bins,q2Bins_half_lhcb,dB_dq_lhcb,q2Bins_err_lhcb,dB_dq_lhcb_stat);
  TGraphErrors* gr_syst_lhcb = new TGraphErrors(n_lhcb_bins,q2Bins_half_lhcb,dB_dq_lhcb,q2Bins_err_lhcb,dB_dq_lhcb_syst);
  gr_stat_lhcb->SetLineColor(kCyan);
  gr_syst_lhcb->SetLineColor(kMagenta);
  gr_stat_lhcb->SetMarkerColor(kCyan);
  gr_syst_lhcb->SetMarkerColor(kMagenta);

  TMultiGraph* mg_2016 = new TMultiGraph();

  TGraphErrors* gr_stat_2016 = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_2016,q2Bins_err,dB_dq_stat_2016);
  gr_stat_2016->SetLineColor(kBlack);

  TGraphErrors* gr_syst_2016 = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_2016,q2Bins_err,dB_dq_syst_2016);
  gr_syst_2016->SetLineColor(kRed);

  mg_2016->Add(gr_syst_2016);
  mg_2016->Add(gr_stat_2016);
  mg_2016->Add(gr_stat_8TeV);
  mg_2016->Add(gr_syst_8TeV);
  mg_2016->Add(gr_stat_lhcb);
  mg_2016->Add(gr_syst_lhcb);

  TLegend *leg_2016 = new TLegend(0.1, 0.1, 0.4, 0.3);
  leg_2016->SetFillColor(0);
  leg_2016->AddEntry(gr_stat_2016, "Stat. Error (CMS Run 2)", "lp");
  leg_2016->AddEntry(gr_syst_2016, "Syst. Error (CMS Run 2)", "lp");
  leg_2016->AddEntry(gr_stat_8TeV, "Stat. Error (CMS Run 1)", "lp");
  leg_2016->AddEntry(gr_syst_8TeV, "Syst. Error (CMS Run 1)", "lp");
  leg_2016->AddEntry(gr_stat_lhcb, "Stat. Error (LHCb Run 1)", "lp");
  leg_2016->AddEntry(gr_syst_lhcb, "Syst. Error (LHCb Run 1)", "lp");

  TCanvas c_2016;
  c_2016.cd();

  mg_2016->Draw("AP");
  leg_2016->Draw("same");
  mg_2016->SetTitle("Year 2016");
  mg_2016->GetYaxis()->SetTitle("dB (B^{0} #rightarrow K^{*0} #mu^{+} #mu^{-})/ dq^{2}");
  mg_2016->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  double y_2016[5] = {mg_2016->GetYaxis()->GetXmin(), mg_2016->GetYaxis()->GetXmin(), mg_2016->GetYaxis()->GetXmax(), mg_2016->GetYaxis()->GetXmax(), mg_2016->GetYaxis()->GetXmin()};

  TGraph* gr6J = new TGraph(5,JPsi, y_2016);
  gr6J->SetLineColor(kBlack);
  gr6J->SetFillColorAlpha(kBlack,0.25);

  TGraph* gr6P = new TGraph(5,Psi, y_2016);
  gr6P->SetLineColor(kBlack);
  gr6P->SetFillColorAlpha(kBlack,0.25);

  gr6J->Draw("F");
  gr6P->Draw("F");

  c_2016.SaveAs("~/public/UML-fit/Branching_fraction/diff_b_2016.gif");
  c_2016.SaveAs("~/public/UML-fit/Branching_fraction/diff_b_2016.pdf");

  ///

  TMultiGraph* mg_2017 = new TMultiGraph();

  TGraphErrors* gr_stat_2017 = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_2017,q2Bins_err,dB_dq_stat_2017);
  gr_stat_2017->SetLineColor(kBlack);

  TGraphErrors* gr_syst_2017 = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_2017,q2Bins_err,dB_dq_syst_2017);
  gr_syst_2017->SetLineColor(kRed);

  mg_2017->Add(gr_syst_2017);
  mg_2017->Add(gr_stat_2017);
  mg_2017->Add(gr_stat_8TeV);
  mg_2017->Add(gr_syst_8TeV);
  mg_2017->Add(gr_stat_lhcb);
  mg_2017->Add(gr_syst_lhcb);

  TLegend *leg_2017 = new TLegend(0.1, 0.1, 0.4, 0.3);
  leg_2017->SetFillColor(0);
  leg_2017->AddEntry(gr_stat_2017, "Stat. Error (CMS Run 2)", "lp");
  leg_2017->AddEntry(gr_syst_2017, "Syst. Error (CMS Run 2)", "lp");
  leg_2017->AddEntry(gr_stat_8TeV, "Stat. Error (CMS Run 1)", "lp");
  leg_2017->AddEntry(gr_syst_8TeV, "Syst. Error (CMS Run 1)", "lp");
  leg_2017->AddEntry(gr_stat_lhcb, "Stat. Error (LHCb Run 1)", "lp");
  leg_2017->AddEntry(gr_syst_lhcb, "Syst. Error (LHCb Run 1)", "lp");

  TCanvas c_2017;
  c_2017.cd();

  mg_2017->Draw("AP");
  leg_2017->Draw("same");
  mg_2017->SetTitle("Year 2017");
  mg_2017->GetYaxis()->SetTitle("dB (B^{0} #rightarrow K^{*0} #mu^{+} #mu^{-})/ dq^{2}");
  mg_2017->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  double y_2017[5] = {mg_2017->GetYaxis()->GetXmin(), mg_2017->GetYaxis()->GetXmin(), mg_2017->GetYaxis()->GetXmax(), mg_2017->GetYaxis()->GetXmax(), mg_2017->GetYaxis()->GetXmin()};

  TGraph* gr7J = new TGraph(5,JPsi, y_2017);
  gr7J->SetLineColor(kBlack);
  gr7J->SetFillColorAlpha(kBlack,0.25);

  TGraph* gr7P = new TGraph(5,Psi, y_2017);
  gr7P->SetLineColor(kBlack);
  gr7P->SetFillColorAlpha(kBlack,0.25);

  gr7J->Draw("F");
  gr7P->Draw("F");

  c_2017.SaveAs("~/public/UML-fit/Branching_fraction/diff_b_2017.gif");
  c_2017.SaveAs("~/public/UML-fit/Branching_fraction/diff_b_2017.pdf");

  ///

  TMultiGraph* mg_2018 = new TMultiGraph();

  TGraphErrors* gr_stat_2018 = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_2018,q2Bins_err,dB_dq_stat_2018);
  gr_stat_2018->SetLineColor(kBlack);

  TGraphErrors* gr_syst_2018 = new TGraphErrors(n_q2Bin,q2Bins_half,dB_dq_2018,q2Bins_err,dB_dq_syst_2018);
  gr_syst_2018->SetLineColor(kRed);

  mg_2018->Add(gr_syst_2018);
  mg_2018->Add(gr_stat_2018);
  mg_2018->Add(gr_stat_8TeV);
  mg_2018->Add(gr_syst_8TeV);
  mg_2018->Add(gr_stat_lhcb);
  mg_2018->Add(gr_syst_lhcb);

  TLegend *leg_2018 = new TLegend(0.1, 0.1, 0.4, 0.3);
  leg_2018->SetFillColor(0);
  leg_2018->AddEntry(gr_stat_2018, "Stat. Error (CMS Run 2)", "lp");
  leg_2018->AddEntry(gr_syst_2018, "Syst. Error (CMS Run 2)", "lp");
  leg_2018->AddEntry(gr_stat_8TeV, "Stat. Error (CMS Run 1)", "lp");
  leg_2018->AddEntry(gr_syst_8TeV, "Syst. Error (CMS Run 1)", "lp");
  leg_2018->AddEntry(gr_stat_lhcb, "Stat. Error (LHCb Run 1)", "lp");
  leg_2018->AddEntry(gr_syst_lhcb, "Syst. Error (LHCb Run 1)", "lp");

  TCanvas c_2018;
  c_2018.cd();

  mg_2018->Draw("AP");
  leg_2018->Draw("same");
  mg_2018->SetTitle("Year 2018");
  mg_2018->GetYaxis()->SetTitle("dB (B^{0} #rightarrow K^{*0} #mu^{+} #mu^{-})/ dq^{2}");
  mg_2018->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  double y_2018[5] = {mg_2018->GetYaxis()->GetXmin(), mg_2018->GetYaxis()->GetXmin(), mg_2018->GetYaxis()->GetXmax(), mg_2018->GetYaxis()->GetXmax(), mg_2018->GetYaxis()->GetXmin()};

  TGraph* gr8J = new TGraph(5,JPsi, y_2018);
  gr8J->SetLineColor(kBlack);
  gr8J->SetFillColorAlpha(kBlack,0.25);

  TGraph* gr8P = new TGraph(5,Psi, y_2018);
  gr8P->SetLineColor(kBlack);
  gr8P->SetFillColorAlpha(kBlack,0.25);

  gr8J->Draw("F");
  gr8P->Draw("F");

  c_2018.SaveAs("~/public/UML-fit/Branching_fraction/diff_b_2018.gif");
  c_2018.SaveAs("~/public/UML-fit/Branching_fraction/diff_b_2018.pdf");

  double yields_CMS[] = {84., 145., 117., 254., 163583., 362., 10508., 225.};
  double yields_CMS_stat[] = {11., 16., 15., 21., 503., 25., 135., 18.};

  TMultiGraph* mg1 = new TMultiGraph();
  
  double bins[] = {0, 1, 2, 3, 4, 5, 6, 7};

  double relative_2016[n_q2Bin];
  double relative_2017[n_q2Bin];
  double relative_2018[n_q2Bin];
  double relative[n_q2Bin]; 
  double CMS_relative[n_q2Bin];

  double relative_syst_2016[n_q2Bin];
  double relative_syst_2017[n_q2Bin];
  double relative_syst_2018[n_q2Bin];
  double relative_syst[n_q2Bin];

  double total_2016[n_q2Bin];
  double total_2017[n_q2Bin];
  double total_2018[n_q2Bin];
  double total[n_q2Bin];

  for(int i = 0; i < n_q2Bin; i++){
    relative_2016[i] = yields_errors_2016[i]/yields_2016[i];
    relative_2017[i] = yields_errors_2017[i]/yields_2017[i];
    relative_2018[i] = yields_errors_2018[i]/yields_2018[i];
    relative[i] = yields_average_stat[i]/yields_average[i];
    CMS_relative[i] = yields_CMS_stat[i]/yields_CMS[i];

    relative_syst_2016[i] = yield_syst_2016[i]/yields_2016[i];
    relative_syst_2017[i] = yield_syst_2017[i]/yields_2017[i];
    relative_syst_2018[i] = yield_syst_2018[i]/yields_2018[i];
    relative_syst[i] = yields_average_syst[i]/yields_average[i];
  }
 
  TGraph* year_2016 = new TGraph(n_q2Bin, bins, relative_2016);
  TGraph* year_2017 = new TGraph(n_q2Bin, bins, relative_2017);
  TGraph* year_2018 = new TGraph(n_q2Bin, bins, relative_2018);
  TGraph* year_average = new TGraph(n_q2Bin, bins, relative);
  TGraph* year_CMS = new TGraph(n_q2Bin, bins, CMS_relative);

  year_2016->SetLineColor(kBlue);
  year_2017->SetLineColor(kOrange+7);
  year_2018->SetLineColor(kGreen+3);
  year_average->SetLineColor(kBlack);
  year_CMS->SetLineColor(kRed);

  year_2016->SetLineWidth(2); 
  year_2017->SetLineWidth(2);
  year_2018->SetLineWidth(2);
  year_average->SetLineWidth(2);
  year_CMS->SetLineWidth(2);

  year_2016->SetMarkerStyle(21);
  year_2017->SetMarkerStyle(22);
  year_2018->SetMarkerStyle(20);
  year_average->SetMarkerStyle(23);
  year_CMS->SetMarkerStyle(33);

  year_2016->SetMarkerColor(kBlue);
  year_2017->SetMarkerColor(kOrange+7);
  year_2018->SetMarkerColor(kGreen+3);
  year_average->SetMarkerColor(kBlack);
  year_CMS->SetMarkerColor(kRed);

  mg1->Add(year_2016);
  mg1->Add(year_2017);
  mg1->Add(year_2018);
  mg1->Add(year_average);
  mg1->Add(year_CMS);

  TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg2->SetFillColor(0);
  leg2->AddEntry(year_2016, "2016", "lp");
  leg2->AddEntry(year_2017, "2017", "lp");
  leg2->AddEntry(year_2018, "2018", "lp"); 
  leg2->AddEntry(year_average, "2016+2017+2018", "lp");     
  leg2->AddEntry(year_CMS, "2012", "lp");     

  TCanvas c2;
  c2.cd();
  c2.SetGrid();
  mg1->Draw("ALP"); 
  mg1->SetTitle("Yields statistical error (relative)");
  mg1->GetXaxis()->SetTitle("q^{2} bins");
  mg1->GetYaxis()->SetTitle("Relative statistical error");
  mg1->GetYaxis()->SetRangeUser(0.,0.2);
  leg2->Draw("same");

  c2.SaveAs("~/public/UML-fit/Branching_fraction/stat_errors.gif");
  c2.SaveAs("~/public/UML-fit/Branching_fraction/stat_errors.pdf");

  TMultiGraph* mg2_2016 = new TMultiGraph();
  TMultiGraph* mg2_2017 = new TMultiGraph();
  TMultiGraph* mg2_2018 = new TMultiGraph();
  TMultiGraph* mg2 = new TMultiGraph();

  TGraph* Y_2016_stat = new TGraph(n_q2Bin, bins, relative_2016);
  TGraph* Y_2017_stat = new TGraph(n_q2Bin, bins, relative_2017);
  TGraph* Y_2018_stat = new TGraph(n_q2Bin, bins, relative_2018);
  TGraph* Y_stat = new TGraph(n_q2Bin, bins, relative);

  TGraph* Y_2016_syst = new TGraph(n_q2Bin, bins, relative_syst_2016);
  TGraph* Y_2017_syst = new TGraph(n_q2Bin, bins, relative_syst_2017);
  TGraph* Y_2018_syst = new TGraph(n_q2Bin, bins, relative_syst_2018);
  TGraph* Y_syst = new TGraph(n_q2Bin, bins, relative_syst);

  Y_2016_stat->SetLineColor(kBlack);
  Y_2016_syst->SetLineColor(kRed);
  Y_2016_stat->SetMarkerColor(kBlack);
  Y_2016_syst->SetMarkerColor(kRed);
  
  Y_2017_stat->SetLineColor(kBlack);
  Y_2017_syst->SetLineColor(kRed);
  Y_2017_stat->SetMarkerColor(kBlack);
  Y_2017_syst->SetMarkerColor(kRed); 

  Y_2018_stat->SetLineColor(kBlack);
  Y_2018_syst->SetLineColor(kRed);
  Y_2018_stat->SetMarkerColor(kBlack);
  Y_2018_syst->SetMarkerColor(kRed);

  Y_stat->SetLineColor(kBlack);
  Y_syst->SetLineColor(kRed);
  Y_stat->SetMarkerColor(kBlack);
  Y_syst->SetMarkerColor(kRed);

  Y_2016_stat->SetLineWidth(2);
  Y_2017_stat->SetLineWidth(2);
  Y_2018_stat->SetLineWidth(2);
  Y_stat->SetLineWidth(2);

  Y_2016_syst->SetLineWidth(2);
  Y_2017_syst->SetLineWidth(2);
  Y_2018_syst->SetLineWidth(2);
  Y_syst->SetLineWidth(2);

  mg2_2016->Add(Y_2016_stat);
  mg2_2017->Add(Y_2017_stat);
  mg2_2018->Add(Y_2018_stat);
  mg2->Add(Y_stat);

  mg2_2016->Add(Y_2016_syst);
  mg2_2017->Add(Y_2017_syst);
  mg2_2018->Add(Y_2018_syst);
  mg2->Add(Y_syst);

  TLegend *leg3_2016 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg3_2016->SetFillColor(0);
  leg3_2016->AddEntry(Y_2016_stat, "Stat. error", "lp");
  leg3_2016->AddEntry(Y_2016_syst, "Syst. error", "lp");

  TLegend *leg3_2017 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg3_2017->SetFillColor(0);
  leg3_2017->AddEntry(Y_2017_stat, "Stat. error", "lp");
  leg3_2017->AddEntry(Y_2017_syst, "Syst. error", "lp");

  TLegend *leg3_2018 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg3_2018->SetFillColor(0);
  leg3_2018->AddEntry(Y_2018_stat, "Stat. error", "lp");
  leg3_2018->AddEntry(Y_2018_syst, "Syst. error", "lp");

  TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg3->SetFillColor(0);
  leg3->AddEntry(Y_stat, "Stat. error", "lp");
  leg3->AddEntry(Y_syst, "Syst. error", "lp");

  TLatex* latex_2016_0 = new TLatex(Y_2016_stat->GetX() [0], Y_2016_stat->GetY() [0], Form("%.1f %%",relative_2016[0]*100));
  TLatex* latex_2016_1 = new TLatex(Y_2016_stat->GetX() [1], Y_2016_stat->GetY() [1], Form("%.1f %%",relative_2016[1]*100));
  TLatex* latex_2016_2 = new TLatex(Y_2016_stat->GetX() [2], Y_2016_stat->GetY() [2], Form("%.1f %%",relative_2016[2]*100));
  TLatex* latex_2016_3 = new TLatex(Y_2016_stat->GetX() [3], Y_2016_stat->GetY() [3], Form("%.1f %%",relative_2016[3]*100));
  TLatex* latex_2016_4 = new TLatex(Y_2016_stat->GetX() [4], Y_2016_stat->GetY() [4], Form("%.1f %%",relative_2016[4]*100));
  TLatex* latex_2016_5 = new TLatex(Y_2016_stat->GetX() [5], Y_2016_stat->GetY() [5], Form("%.1f %%",relative_2016[5]*100));
  TLatex* latex_2016_6 = new TLatex(Y_2016_stat->GetX() [6], Y_2016_stat->GetY() [6], Form("%.1f %%",relative_2016[6]*100));
  TLatex* latex_2016_7 = new TLatex(Y_2016_stat->GetX() [7], Y_2016_stat->GetY() [7], Form("%.1f %%",relative_2016[7]*100));

  TLatex* latex_2016_0s = new TLatex(Y_2016_syst->GetX() [0], Y_2016_syst->GetY() [0], Form("%.1f %%",relative_syst_2016[0]*100));
  TLatex* latex_2016_1s = new TLatex(Y_2016_syst->GetX() [1], Y_2016_syst->GetY() [1], Form("%.1f %%",relative_syst_2016[1]*100));
  TLatex* latex_2016_2s = new TLatex(Y_2016_syst->GetX() [2], Y_2016_syst->GetY() [2], Form("%.1f %%",relative_syst_2016[2]*100));
  TLatex* latex_2016_3s = new TLatex(Y_2016_syst->GetX() [3], Y_2016_syst->GetY() [3], Form("%.1f %%",relative_syst_2016[3]*100));
  TLatex* latex_2016_4s = new TLatex(Y_2016_syst->GetX() [4], Y_2016_syst->GetY() [4], Form("%.1f %%",relative_syst_2016[4]*100));
  TLatex* latex_2016_5s = new TLatex(Y_2016_syst->GetX() [5], Y_2016_syst->GetY() [5], Form("%.1f %%",relative_syst_2016[5]*100));
  TLatex* latex_2016_6s = new TLatex(Y_2016_syst->GetX() [6], Y_2016_syst->GetY() [6], Form("%.1f %%",relative_syst_2016[6]*100));
  TLatex* latex_2016_7s = new TLatex(Y_2016_syst->GetX() [7], Y_2016_syst->GetY() [7], Form("%.1f %%",relative_syst_2016[7]*100));

  TLatex* latex_2017_0 = new TLatex(Y_2017_stat->GetX() [0], Y_2017_stat->GetY() [0], Form("%.1f %%",relative_2017[0]*100));
  TLatex* latex_2017_1 = new TLatex(Y_2017_stat->GetX() [1], Y_2017_stat->GetY() [1], Form("%.1f %%",relative_2017[1]*100));
  TLatex* latex_2017_2 = new TLatex(Y_2017_stat->GetX() [2], Y_2017_stat->GetY() [2], Form("%.1f %%",relative_2017[2]*100));
  TLatex* latex_2017_3 = new TLatex(Y_2017_stat->GetX() [3], Y_2017_stat->GetY() [3], Form("%.1f %%",relative_2017[3]*100));
  TLatex* latex_2017_4 = new TLatex(Y_2017_stat->GetX() [4], Y_2017_stat->GetY() [4], Form("%.1f %%",relative_2017[4]*100));
  TLatex* latex_2017_5 = new TLatex(Y_2017_stat->GetX() [5], Y_2017_stat->GetY() [5], Form("%.1f %%",relative_2017[5]*100));
  TLatex* latex_2017_6 = new TLatex(Y_2017_stat->GetX() [6], Y_2017_stat->GetY() [6], Form("%.1f %%",relative_2017[6]*100));
  TLatex* latex_2017_7 = new TLatex(Y_2017_stat->GetX() [7], Y_2017_stat->GetY() [7], Form("%.1f %%",relative_2017[7]*100));

  TLatex* latex_2017_0s = new TLatex(Y_2017_syst->GetX() [0], Y_2017_syst->GetY() [0], Form("%.1f %%",relative_syst_2017[0]*100));
  TLatex* latex_2017_1s = new TLatex(Y_2017_syst->GetX() [1], Y_2017_syst->GetY() [1], Form("%.1f %%",relative_syst_2017[1]*100));
  TLatex* latex_2017_2s = new TLatex(Y_2017_syst->GetX() [2], Y_2017_syst->GetY() [2], Form("%.1f %%",relative_syst_2017[2]*100));
  TLatex* latex_2017_3s = new TLatex(Y_2017_syst->GetX() [3], Y_2017_syst->GetY() [3], Form("%.1f %%",relative_syst_2017[3]*100));
  TLatex* latex_2017_4s = new TLatex(Y_2017_syst->GetX() [4], Y_2017_syst->GetY() [4], Form("%.1f %%",relative_syst_2017[4]*100));
  TLatex* latex_2017_5s = new TLatex(Y_2017_syst->GetX() [5], Y_2017_syst->GetY() [5], Form("%.1f %%",relative_syst_2017[5]*100));
  TLatex* latex_2017_6s = new TLatex(Y_2017_syst->GetX() [6], Y_2017_syst->GetY() [6], Form("%.1f %%",relative_syst_2017[6]*100));
  TLatex* latex_2017_7s = new TLatex(Y_2017_syst->GetX() [7], Y_2017_syst->GetY() [7], Form("%.1f %%",relative_syst_2017[7]*100));

  TLatex* latex_2018_0 = new TLatex(Y_2018_stat->GetX() [0], Y_2018_stat->GetY() [0], Form("%.1f %%",relative_2018[0]*100));
  TLatex* latex_2018_1 = new TLatex(Y_2018_stat->GetX() [1], Y_2018_stat->GetY() [1], Form("%.1f %%",relative_2018[1]*100));
  TLatex* latex_2018_2 = new TLatex(Y_2018_stat->GetX() [2], Y_2018_stat->GetY() [2], Form("%.1f %%",relative_2018[2]*100));
  TLatex* latex_2018_3 = new TLatex(Y_2018_stat->GetX() [3], Y_2018_stat->GetY() [3], Form("%.1f %%",relative_2018[3]*100));
  TLatex* latex_2018_4 = new TLatex(Y_2018_stat->GetX() [4], Y_2018_stat->GetY() [4], Form("%.1f %%",relative_2018[4]*100));
  TLatex* latex_2018_5 = new TLatex(Y_2018_stat->GetX() [5], Y_2018_stat->GetY() [5], Form("%.1f %%",relative_2018[5]*100));
  TLatex* latex_2018_6 = new TLatex(Y_2018_stat->GetX() [6], Y_2018_stat->GetY() [6], Form("%.1f %%",relative_2018[6]*100));
  TLatex* latex_2018_7 = new TLatex(Y_2018_stat->GetX() [7], Y_2018_stat->GetY() [7], Form("%.1f %%",relative_2018[7]*100));

  TLatex* latex_2018_0s = new TLatex(Y_2018_syst->GetX() [0], Y_2018_syst->GetY() [0], Form("%.1f %%",relative_syst_2018[0]*100));
  TLatex* latex_2018_1s = new TLatex(Y_2018_syst->GetX() [1], Y_2018_syst->GetY() [1], Form("%.1f %%",relative_syst_2018[1]*100));
  TLatex* latex_2018_2s = new TLatex(Y_2018_syst->GetX() [2], Y_2018_syst->GetY() [2], Form("%.1f %%",relative_syst_2018[2]*100));
  TLatex* latex_2018_3s = new TLatex(Y_2018_syst->GetX() [3], Y_2018_syst->GetY() [3], Form("%.1f %%",relative_syst_2018[3]*100));
  TLatex* latex_2018_4s = new TLatex(Y_2018_syst->GetX() [4], Y_2018_syst->GetY() [4], Form("%.1f %%",relative_syst_2018[4]*100));
  TLatex* latex_2018_5s = new TLatex(Y_2018_syst->GetX() [5], Y_2018_syst->GetY() [5], Form("%.1f %%",relative_syst_2018[5]*100));
  TLatex* latex_2018_6s = new TLatex(Y_2018_syst->GetX() [6], Y_2018_syst->GetY() [6], Form("%.1f %%",relative_syst_2018[6]*100));
  TLatex* latex_2018_7s = new TLatex(Y_2018_syst->GetX() [7], Y_2018_syst->GetY() [7], Form("%.1f %%",relative_syst_2018[7]*100));

  TLatex* latex_0 = new TLatex(Y_stat->GetX() [0], Y_stat->GetY() [0], Form("%.1f %%",relative[0]*100));
  TLatex* latex_1 = new TLatex(Y_stat->GetX() [1], Y_stat->GetY() [1], Form("%.1f %%",relative[1]*100));
  TLatex* latex_2 = new TLatex(Y_stat->GetX() [2], Y_stat->GetY() [2], Form("%.1f %%",relative[2]*100));
  TLatex* latex_3 = new TLatex(Y_stat->GetX() [3], Y_stat->GetY() [3], Form("%.1f %%",relative[3]*100));
  TLatex* latex_4 = new TLatex(Y_stat->GetX() [4], Y_stat->GetY() [4], Form("%.1f %%",relative[4]*100));
  TLatex* latex_5 = new TLatex(Y_stat->GetX() [5], Y_stat->GetY() [5], Form("%.1f %%",relative[5]*100));
  TLatex* latex_6 = new TLatex(Y_stat->GetX() [6], Y_stat->GetY() [6], Form("%.1f %%",relative[6]*100));
  TLatex* latex_7 = new TLatex(Y_stat->GetX() [7], Y_stat->GetY() [7], Form("%.1f %%",relative[7]*100));

  TLatex* latex_0s = new TLatex(Y_syst->GetX() [0], Y_syst->GetY() [0], Form("%.1f %%",relative_syst[0]*100));
  TLatex* latex_1s = new TLatex(Y_syst->GetX() [1], Y_syst->GetY() [1], Form("%.1f %%",relative_syst[1]*100));
  TLatex* latex_2s = new TLatex(Y_syst->GetX() [2], Y_syst->GetY() [2], Form("%.1f %%",relative_syst[2]*100));
  TLatex* latex_3s = new TLatex(Y_syst->GetX() [3], Y_syst->GetY() [3], Form("%.1f %%",relative_syst[3]*100));
  TLatex* latex_4s = new TLatex(Y_syst->GetX() [4], Y_syst->GetY() [4], Form("%.1f %%",relative_syst[4]*100));
  TLatex* latex_5s = new TLatex(Y_syst->GetX() [5], Y_syst->GetY() [5], Form("%.1f %%",relative_syst[5]*100));
  TLatex* latex_6s = new TLatex(Y_syst->GetX() [6], Y_syst->GetY() [6], Form("%.1f %%",relative_syst[6]*100));
  TLatex* latex_7s = new TLatex(Y_syst->GetX() [7], Y_syst->GetY() [7], Form("%.1f %%",relative_syst[7]*100));

  latex_2016_0->SetTextSize(0.03); 
  latex_2016_1->SetTextSize(0.03);
  latex_2016_2->SetTextSize(0.03);
  latex_2016_3->SetTextSize(0.03);
  latex_2016_4->SetTextSize(0.03);
  latex_2016_5->SetTextSize(0.03);
  latex_2016_6->SetTextSize(0.03);
  latex_2016_7->SetTextSize(0.03);

  latex_2016_0s->SetTextSize(0.03);
  latex_2016_1s->SetTextSize(0.03);
  latex_2016_2s->SetTextSize(0.03);
  latex_2016_3s->SetTextSize(0.03);
  latex_2016_4s->SetTextSize(0.03);
  latex_2016_5s->SetTextSize(0.03);
  latex_2016_6s->SetTextSize(0.03);
  latex_2016_7s->SetTextSize(0.03);

  latex_2017_0->SetTextSize(0.03);
  latex_2017_1->SetTextSize(0.03);
  latex_2017_2->SetTextSize(0.03);
  latex_2017_3->SetTextSize(0.03);
  latex_2017_4->SetTextSize(0.03);
  latex_2017_5->SetTextSize(0.03);
  latex_2017_6->SetTextSize(0.03);
  latex_2017_7->SetTextSize(0.03);

  latex_2017_0s->SetTextSize(0.03);
  latex_2017_1s->SetTextSize(0.03);
  latex_2017_2s->SetTextSize(0.03);
  latex_2017_3s->SetTextSize(0.03);
  latex_2017_4s->SetTextSize(0.03);
  latex_2017_5s->SetTextSize(0.03);
  latex_2017_6s->SetTextSize(0.03);
  latex_2017_7s->SetTextSize(0.03);

  latex_2018_0->SetTextSize(0.03);
  latex_2018_1->SetTextSize(0.03);
  latex_2018_2->SetTextSize(0.03);
  latex_2018_3->SetTextSize(0.03);
  latex_2018_4->SetTextSize(0.03);
  latex_2018_5->SetTextSize(0.03);
  latex_2018_6->SetTextSize(0.03);
  latex_2018_7->SetTextSize(0.03);

  latex_2018_0s->SetTextSize(0.03);
  latex_2018_1s->SetTextSize(0.03);
  latex_2018_2s->SetTextSize(0.03);
  latex_2018_3s->SetTextSize(0.03);
  latex_2018_4s->SetTextSize(0.03);
  latex_2018_5s->SetTextSize(0.03);
  latex_2018_6s->SetTextSize(0.03);
  latex_2018_7s->SetTextSize(0.03);

  latex_0->SetTextSize(0.03);
  latex_1->SetTextSize(0.03);
  latex_2->SetTextSize(0.03);
  latex_3->SetTextSize(0.03);
  latex_4->SetTextSize(0.03);
  latex_5->SetTextSize(0.03);
  latex_6->SetTextSize(0.03);
  latex_7->SetTextSize(0.03);

  latex_0s->SetTextSize(0.03);
  latex_1s->SetTextSize(0.03);
  latex_2s->SetTextSize(0.03);
  latex_3s->SetTextSize(0.03);
  latex_4s->SetTextSize(0.03);
  latex_5s->SetTextSize(0.03);
  latex_6s->SetTextSize(0.03);
  latex_7s->SetTextSize(0.03);

  latex_2016_0->SetTextColor(kBlack);
  latex_2016_1->SetTextColor(kBlack);
  latex_2016_2->SetTextColor(kBlack);
  latex_2016_3->SetTextColor(kBlack);
  latex_2016_4->SetTextColor(kBlack);
  latex_2016_5->SetTextColor(kBlack);
  latex_2016_6->SetTextColor(kBlack);
  latex_2016_7->SetTextColor(kBlack);

  latex_2016_0s->SetTextColor(kRed);
  latex_2016_1s->SetTextColor(kRed);
  latex_2016_2s->SetTextColor(kRed);
  latex_2016_3s->SetTextColor(kRed);
  latex_2016_4s->SetTextColor(kRed);
  latex_2016_5s->SetTextColor(kRed);
  latex_2016_6s->SetTextColor(kRed);
  latex_2016_7s->SetTextColor(kRed);

  latex_2017_0->SetTextColor(kBlack);
  latex_2017_1->SetTextColor(kBlack);
  latex_2017_2->SetTextColor(kBlack);
  latex_2017_3->SetTextColor(kBlack);
  latex_2017_4->SetTextColor(kBlack);
  latex_2017_5->SetTextColor(kBlack);
  latex_2017_6->SetTextColor(kBlack);
  latex_2017_7->SetTextColor(kBlack);

  latex_2017_0s->SetTextColor(kRed);
  latex_2017_1s->SetTextColor(kRed);
  latex_2017_2s->SetTextColor(kRed);
  latex_2017_3s->SetTextColor(kRed);
  latex_2017_4s->SetTextColor(kRed);
  latex_2017_5s->SetTextColor(kRed);
  latex_2017_6s->SetTextColor(kRed);
  latex_2017_7s->SetTextColor(kRed);

  latex_2018_0->SetTextColor(kBlack);
  latex_2018_1->SetTextColor(kBlack);
  latex_2018_2->SetTextColor(kBlack);
  latex_2018_3->SetTextColor(kBlack);
  latex_2018_4->SetTextColor(kBlack);
  latex_2018_5->SetTextColor(kBlack);
  latex_2018_6->SetTextColor(kBlack);
  latex_2018_7->SetTextColor(kBlack);

  latex_2018_0s->SetTextColor(kRed);
  latex_2018_1s->SetTextColor(kRed);
  latex_2018_2s->SetTextColor(kRed);
  latex_2018_3s->SetTextColor(kRed);
  latex_2018_4s->SetTextColor(kRed);
  latex_2018_5s->SetTextColor(kRed);
  latex_2018_6s->SetTextColor(kRed);
  latex_2018_7s->SetTextColor(kRed);

  latex_0->SetTextColor(kBlack);
  latex_1->SetTextColor(kBlack);
  latex_2->SetTextColor(kBlack);
  latex_3->SetTextColor(kBlack);
  latex_4->SetTextColor(kBlack);
  latex_5->SetTextColor(kBlack);
  latex_6->SetTextColor(kBlack);
  latex_7->SetTextColor(kBlack);

  latex_0s->SetTextColor(kRed);
  latex_1s->SetTextColor(kRed);
  latex_2s->SetTextColor(kRed);
  latex_3s->SetTextColor(kRed);
  latex_4s->SetTextColor(kRed);
  latex_5s->SetTextColor(kRed);
  latex_6s->SetTextColor(kRed);
  latex_7s->SetTextColor(kRed);

  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_0);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_1);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_2);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_3);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_4);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_5);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_6);
  Y_2016_stat->GetListOfFunctions()->Add(latex_2016_7);

  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_0s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_1s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_2s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_3s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_4s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_5s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_6s);
  Y_2016_syst->GetListOfFunctions()->Add(latex_2016_7s);

  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_0);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_1);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_2);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_3);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_4);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_5);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_6);
  Y_2017_stat->GetListOfFunctions()->Add(latex_2017_7);

  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_0s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_1s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_2s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_3s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_4s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_5s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_6s);
  Y_2017_syst->GetListOfFunctions()->Add(latex_2017_7s);

  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_0);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_1);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_2);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_3);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_4);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_5);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_6);
  Y_2018_stat->GetListOfFunctions()->Add(latex_2018_7);

  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_0s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_1s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_2s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_3s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_4s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_5s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_6s);
  Y_2018_syst->GetListOfFunctions()->Add(latex_2018_7s);

  Y_stat->GetListOfFunctions()->Add(latex_0);
  Y_stat->GetListOfFunctions()->Add(latex_1);
  Y_stat->GetListOfFunctions()->Add(latex_2);
  Y_stat->GetListOfFunctions()->Add(latex_3);
  Y_stat->GetListOfFunctions()->Add(latex_4);
  Y_stat->GetListOfFunctions()->Add(latex_5);
  Y_stat->GetListOfFunctions()->Add(latex_6);
  Y_stat->GetListOfFunctions()->Add(latex_7);

  Y_syst->GetListOfFunctions()->Add(latex_0s);
  Y_syst->GetListOfFunctions()->Add(latex_1s);
  Y_syst->GetListOfFunctions()->Add(latex_2s);
  Y_syst->GetListOfFunctions()->Add(latex_3s);
  Y_syst->GetListOfFunctions()->Add(latex_4s);
  Y_syst->GetListOfFunctions()->Add(latex_5s);
  Y_syst->GetListOfFunctions()->Add(latex_6s);
  Y_syst->GetListOfFunctions()->Add(latex_7s);

  TCanvas c3_2016;
  c3_2016.cd();
  c3_2016.SetGrid();
  mg2_2016->Draw("AL");
  mg2_2016->SetTitle("Yields (relative) uncertainties - 2016");
  mg2_2016->GetXaxis()->SetTitle("q^{2} bins");
  leg3_2016->Draw("same");
  c3_2016.SaveAs("~/public/UML-fit/Branching_fraction/yields_2016.gif");
  c3_2016.SaveAs("~/public/UML-fit/Branching_fraction/yields_2016.pdf");

  TCanvas c3_2017;
  c3_2017.cd();
  c3_2017.SetGrid();
  mg2_2017->Draw("AL");
  mg2_2017->SetTitle("Yields (relative) uncertainties - 2017");
  mg2_2017->GetXaxis()->SetTitle("q^{2} bins");
  leg3_2017->Draw("same");
  c3_2017.SaveAs("~/public/UML-fit/Branching_fraction/yields_2017.gif");
  c3_2017.SaveAs("~/public/UML-fit/Branching_fraction/yields_2017.pdf");

  TCanvas c3_2018;
  c3_2018.cd();
  c3_2018.SetGrid();
  mg2_2018->Draw("AL");
  mg2_2018->SetTitle("Yields (relative) uncertainties - 2018");
  mg2_2018->GetXaxis()->SetTitle("q^{2} bins");
  leg3_2018->Draw("same");
  c3_2018.SaveAs("~/public/UML-fit/Branching_fraction/yields_2018.gif");
  c3_2018.SaveAs("~/public/UML-fit/Branching_fraction/yields_2018.pdf");

  TCanvas c3;
  c3.cd();
  c3.SetGrid();
  mg2->Draw("AL");
  mg2->SetTitle("Yields (relative) uncertainties - 2016+2017+2018");
  mg2->GetXaxis()->SetTitle("q^{2} bins");
  leg3->Draw("same");
  c3.SaveAs("~/public/UML-fit/Branching_fraction/yields.gif");
  c3.SaveAs("~/public/UML-fit/Branching_fraction/yields.pdf");

}



