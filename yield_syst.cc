#include <RooFitResult.h>
#include <RooRealVar.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <TAttMarker.h>
#include <vector>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TGraphPainter.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

std::ofstream file_stat;
std::ofstream file_syst;

double getMax(double list[5], int n_pdf);

void yield_syst(){

  int n_q2Bin = 8;
  double q2Bins[] = {1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16};
  double q2Bins_half[n_q2Bin];
  double q2Bins_err[n_q2Bin];

  for(int i = 0; i < n_q2Bin; i++){
    q2Bins_half[i] = q2Bins[i] + 0.5* (q2Bins[i+1]-q2Bins[i]);
    q2Bins_err[i] = 0.5* (q2Bins[i+1]-q2Bins[i]);
  }

  double yields_2016[n_q2Bin];
  double yields_2017[n_q2Bin];
  double yields_2018[n_q2Bin];

  double yields_stat_2016[n_q2Bin];
  double yields_stat_2017[n_q2Bin];
  double yields_stat_2018[n_q2Bin];

  double yields_syst_2016[n_q2Bin];
  double yields_syst_2017[n_q2Bin];
  double yields_syst_2018[n_q2Bin];

  cout << '|' << setw(15) << "Year" << '|' << setw(15) << "q2Bin" << '|' << setw(15) << "PDF" << '|' << setw(15) << "Yield" << '|' << setw(15) << "Yield Error" << '|' << setw(15) << "Rel. Stat. (%)" << '|' << setw(15) << "Rel. Diff. (%)" << '|' << endl;

  int n_pdf = 5;

  double signal_yield_2016[n_q2Bin][n_pdf];
  double signal_yield_2017[n_q2Bin][n_pdf];
  double signal_yield_2018[n_q2Bin][n_pdf];

  double stat_yield_2016[n_q2Bin][n_pdf];
  double stat_yield_2017[n_q2Bin][n_pdf];
  double stat_yield_2018[n_q2Bin][n_pdf];

  double syst_yield_2016[n_q2Bin][n_pdf];
  double syst_yield_2017[n_q2Bin][n_pdf];
  double syst_yield_2018[n_q2Bin][n_pdf];

  for(int q2Bin = 0; q2Bin < n_q2Bin; q2Bin++){


    for(int year = 6; year < 9; year++){

      for(int pdf = 0; pdf < n_pdf; pdf++){   
        TFile* f; 
        RooFitResult* fitresult;

        if(pdf == 4){
          f = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c2m0_subs0CT+WT.root",year,q2Bin));
          fitresult = (RooFitResult*)f->Get(Form("simFitResult_b%ip2c2m0subs0",q2Bin));
        }
        else{
          f = new TFile(Form("~/public/UML-fit/simFitMassResults/simFitResult_recoMC_fullMass201%i_DATA_b%ip2c1m%i_subs0CT+WT.root",year,q2Bin,pdf));
          fitresult = (RooFitResult*)f->Get(Form("simFitResult_b%ip2c1m%isubs0",q2Bin,pdf));
        }

        RooRealVar* signal_yield = (RooRealVar*)fitresult->floatParsFinal().find(Form("sig_yield^{201%i}",year));

        if(year == 6){
          signal_yield_2016[q2Bin][pdf] = signal_yield->getVal();
          stat_yield_2016[q2Bin][pdf] = signal_yield->getError();
          syst_yield_2016[q2Bin][pdf] = abs(signal_yield_2016[q2Bin][0] - signal_yield_2016[q2Bin][pdf]);

          cout << '|' << setw(15) << 2016 << '|' << setw(15) << q2Bin << '|' << setw(15) << pdf << '|' << setw(15) <<  signal_yield_2016[q2Bin][pdf]  << '|' << setw(15) << stat_yield_2016[q2Bin][pdf]  << '|' << setw(15) << (stat_yield_2016[q2Bin][pdf]/signal_yield_2016[q2Bin][pdf])*100 << '|' << setw(15) << (syst_yield_2016[q2Bin][pdf]/signal_yield_2016[q2Bin][0])*100 << '|' << endl;
        }
        else if(year == 7){
          signal_yield_2017[q2Bin][pdf] = signal_yield->getVal();
          stat_yield_2017[q2Bin][pdf] = signal_yield->getError();
          syst_yield_2017[q2Bin][pdf] = abs(signal_yield_2017[q2Bin][0] - signal_yield_2017[q2Bin][pdf]);

          cout << '|' << setw(15) << 2017 << '|' << setw(15) << q2Bin << '|' << setw(15) << pdf << '|' << setw(15) << signal_yield_2017[q2Bin][pdf] << '|' << setw(15) << stat_yield_2017[q2Bin][pdf]  << '|' << setw(15) << (stat_yield_2017[q2Bin][pdf]/signal_yield_2017[q2Bin][pdf])*100 << '|' << setw(15) << (syst_yield_2017[q2Bin][pdf]/signal_yield_2017[q2Bin][0])*100 << '|' << endl;
        }
        else if(year == 8){
          signal_yield_2018[q2Bin][pdf] = signal_yield->getVal();
          stat_yield_2018[q2Bin][pdf] = signal_yield->getError();
          syst_yield_2018[q2Bin][pdf] = abs(signal_yield_2018[q2Bin][0] - signal_yield_2018[q2Bin][pdf]);

          cout << '|' << setw(15) << 2018 << '|' << setw(15) << q2Bin << '|' << setw(15) << pdf  << '|' << setw(15) << signal_yield_2018[q2Bin][pdf]  << '|' << setw(15) << stat_yield_2018[q2Bin][pdf] << '|' << setw(15) << (stat_yield_2018[q2Bin][pdf]/signal_yield_2018[q2Bin][pdf])*100 << '|' << setw(15) <<  (syst_yield_2018[q2Bin][pdf]/signal_yield_2018[q2Bin][0])*100 << '|' << endl;
        }        
        delete f;
      }       
 
      if(year == 6){// from nominal fit
        yields_2016[q2Bin] = signal_yield_2016[q2Bin][0];
        yields_stat_2016[q2Bin] = stat_yield_2016[q2Bin][0];

        yields_syst_2016[q2Bin] = getMax(syst_yield_2016[q2Bin],n_pdf);
      }
      else if(year == 7){
        yields_2017[q2Bin] = signal_yield_2017[q2Bin][0];
        yields_stat_2017[q2Bin] = stat_yield_2017[q2Bin][0];
        yields_syst_2017[q2Bin] = getMax(syst_yield_2017[q2Bin],n_pdf);
      } 
      else if(year == 8){
        yields_2018[q2Bin] = signal_yield_2018[q2Bin][0];
        yields_stat_2018[q2Bin] = stat_yield_2018[q2Bin][0];
        yields_syst_2018[q2Bin] = getMax(syst_yield_2018[q2Bin],n_pdf);
     }
    }//ends loop over years  
  }//ends loop over q2Bins


  cout << '|' << setw(15) << "Year" << '|' << setw(15) << "q2Bin" << '|' << setw(15) << "Yield (nominal)" << '|' << setw(15) << "Yield Error" << '|' << setw(15) << "Rel. Stat. (%)" << '|' << setw(15) << "Rel. Syst. (%)" << '|' << endl;

  for(int q2Bin = 0; q2Bin < n_q2Bin; q2Bin++){
    cout << '|' << setw(15) << 2016 << '|' << setw(15) << q2Bin << '|' << setw(15) << yields_2016[q2Bin] << '|' << setw(15) << yields_stat_2016[q2Bin] << '|' << setw(15) << (yields_stat_2016[q2Bin]/yields_2016[q2Bin])*100 << '|' << setw(15) << (yields_syst_2016[q2Bin]/yields_2016[q2Bin])*100 << '|' << endl;
    cout << '|' << setw(15) << 2017 << '|' << setw(15) << q2Bin << '|' << setw(15) << yields_2017[q2Bin] << '|' << setw(15) << yields_stat_2017[q2Bin] << '|' << setw(15) << (yields_stat_2017[q2Bin]/yields_2017[q2Bin])*100 << '|' << setw(15) << (yields_syst_2017[q2Bin]/yields_2017[q2Bin])*100 << '|' << endl;
    cout << '|' << setw(15) << 2018 << '|' << setw(15) << q2Bin << '|' << setw(15) << yields_2018[q2Bin] << '|' << setw(15) << yields_stat_2018[q2Bin] << '|' << setw(15) << (yields_stat_2018[q2Bin]/yields_2018[q2Bin])*100 << '|' << setw(15) << (yields_syst_2018[q2Bin]/yields_2018[q2Bin])*100 << '|' << endl;
  } 

  TFile* f_output_2016 = new TFile("~/public/UML-fit/Systematics/root_files/yield_syst_2016.root", "UPDATE");
  TFile* f_output_2017 = new TFile("~/public/UML-fit/Systematics/root_files/yield_syst_2017.root", "UPDATE");
  TFile* f_output_2018 = new TFile("~/public/UML-fit/Systematics/root_files/yield_syst_2018.root", "UPDATE");

  // 2016 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas c_2016;
  c_2016.SetLogy();
  f_output_2016->cd();
  TMultiGraph* mg_2016 = new TMultiGraph();

  TGraphErrors* g_stat_2016 = new TGraphErrors(n_q2Bin,q2Bins_half,yields_2016,q2Bins_err,yields_stat_2016);
  g_stat_2016->SetMarkerColor(4);
  g_stat_2016->SetMarkerStyle(1);
  g_stat_2016->SetLineColor(1);
  g_stat_2016->Write();

  TGraphErrors* g_syst_2016 = new TGraphErrors(n_q2Bin,q2Bins_half,yields_2016,q2Bins_err,yields_syst_2016);
  g_syst_2016->SetMarkerColor(4);
  g_syst_2016->SetMarkerStyle(1);
  g_syst_2016->SetLineColor(2);
  g_syst_2016->Write();
  
  f_output_2016->Close();
  
  mg_2016->Add(g_syst_2016);
  mg_2016->Add(g_stat_2016);

  TLegend *leg_2016 = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg_2016->SetFillColor(0);
  leg_2016->AddEntry(g_stat_2016, "Statistical Uncertainty", "lp");
  leg_2016->AddEntry(g_syst_2016, "Systematic Uncertainty", "lp");

  mg_2016->Draw("AP");
  leg_2016->Draw("same");
  mg_2016->SetTitle("Signal Yield - 2016");
  mg_2016->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c_2016.SaveAs("~/public/UML-fit/Systematics/plots/yield_syst_2016.gif");
  c_2016.SaveAs("~/public/UML-fit/Systematics/plots/yield_syst_2016.pdf");

  // 2017 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas c_2017;
  c_2017.SetLogy();
  f_output_2017->cd();
  TMultiGraph* mg_2017 = new TMultiGraph();

  TGraphErrors* g_stat_2017 = new TGraphErrors(n_q2Bin,q2Bins_half,yields_2017,q2Bins_err,yields_stat_2017);
  g_stat_2017->SetMarkerColor(4);
  g_stat_2017->SetMarkerStyle(1);
  g_stat_2017->SetLineColor(1);
  g_stat_2017->Write();

  TGraphErrors* g_syst_2017 = new TGraphErrors(n_q2Bin,q2Bins_half,yields_2017,q2Bins_err,yields_syst_2017);
  g_syst_2017->SetMarkerColor(4);
  g_syst_2017->SetMarkerStyle(1);
  g_syst_2017->SetLineColor(2);
  g_syst_2017->Write();

  f_output_2017->Close();

  mg_2017->Add(g_syst_2017);
  mg_2017->Add(g_stat_2017);

  TLegend *leg_2017 = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg_2017->SetFillColor(0);
  leg_2017->AddEntry(g_stat_2017, "Statistical Uncertainty", "lp");
  leg_2017->AddEntry(g_syst_2017, "Systematic Uncertainty", "lp");

  mg_2017->Draw("AP");
  leg_2017->Draw("same");
  mg_2017->SetTitle("Signal Yield - 2017");
  mg_2017->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c_2017.SaveAs("~/public/UML-fit/Systematics/plots/yield_syst_2017.gif");
  c_2017.SaveAs("~/public/UML-fit/Systematics/plots/yield_syst_2017.pdf");

  // 2018 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TCanvas c_2018;
  c_2018.SetLogy(); 
  f_output_2018->cd();
  TMultiGraph* mg_2018 = new TMultiGraph();

  TGraphErrors* g_stat_2018 = new TGraphErrors(n_q2Bin,q2Bins_half,yields_2018,q2Bins_err,yields_stat_2018);
  g_stat_2018->SetMarkerColor(4);
  g_stat_2018->SetMarkerStyle(1);
  g_stat_2018->SetLineColor(1);
  g_stat_2018->Write();

  TGraphErrors* g_syst_2018 = new TGraphErrors(n_q2Bin,q2Bins_half,yields_2018,q2Bins_err,yields_syst_2018);
  g_syst_2018->SetMarkerColor(4);
  g_syst_2018->SetMarkerStyle(1);
  g_syst_2018->SetLineColor(2);
  g_syst_2018->Write();

  f_output_2018->Close();

  mg_2018->Add(g_syst_2018);
  mg_2018->Add(g_stat_2018);

  TLegend *leg_2018 = new TLegend(0.7, 0.8, 0.9, 0.9);
  leg_2018->SetFillColor(0);
  leg_2018->AddEntry(g_stat_2018, "Statistical Uncertainty", "lp");
  leg_2018->AddEntry(g_syst_2018, "Systematic Uncertainty", "lp");

  mg_2018->Draw("AP");
  leg_2018->Draw("same");
  mg_2018->SetTitle("Signal Yield - 2018");
  mg_2018->GetXaxis()->SetTitle("q^{2} [GeV^{2}]");

  c_2018.SaveAs("~/public/UML-fit/Systematics/plots/yield_syst_2018.gif");
  c_2018.SaveAs("~/public/UML-fit/Systematics/plots/yield_syst_2018.pdf");

 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // STAT TABLE
  std::string filename_stat = "/afs/cern.ch/user/m/mfaria/public/UML-fit/Systematics/tables/yield_stat_table.tex";
  file_stat.open(filename_stat);

  if(!file_stat.is_open())
      cout << "Output file_stat not opened!";

  file_stat << "\\begin{table}" << std::endl;
  file_stat << "\\caption{Yields with respective statistical uncertainties for each fit variation.}" << std::endl;
  file_stat << "\\tiny" << std::endl;
  file_stat << "\\centering" << std::endl;
  file_stat << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << std::endl;
  file_stat << "\\hline" << std::endl;

  std::vector<std::string> col_name_stat = {"$q^2$ bin", "Year", "Nominal", "Exp($5.1<m<5.6$)", "Exp($5.0<m<5.5$)", "WT fixed/Exp+erf", "Scale factor"};
  for(int c = 0; c < 6; c++){
    file_stat << col_name_stat[c] << " & ";
  }
  file_stat << col_name_stat[6] << " \\\\ " << std::endl;
  file_stat << "\\hline" << std::endl;

  for(int i = 0; i < n_q2Bin; i++){
    file_stat << " & " << " 2016 & " << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2016[i][0], stat_yield_2016[i][0], stat_yield_2016[i][0]/signal_yield_2016[i][0]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2016[i][1], stat_yield_2016[i][1], stat_yield_2016[i][1]/signal_yield_2016[i][1]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2016[i][2], stat_yield_2016[i][2], stat_yield_2016[i][2]/signal_yield_2016[i][2]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2016[i][3], stat_yield_2016[i][3], stat_yield_2016[i][3]/signal_yield_2016[i][3]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2016[i][4], stat_yield_2016[i][4], stat_yield_2016[i][4]/signal_yield_2016[i][4]) << "\\\\" << std::endl;

    file_stat << Form("%i & ",i) << " 2017 & " << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2017[i][0], stat_yield_2017[i][0], stat_yield_2017[i][0]/signal_yield_2017[i][0]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2017[i][1], stat_yield_2017[i][1], stat_yield_2017[i][1]/signal_yield_2017[i][1]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2017[i][2], stat_yield_2017[i][2], stat_yield_2017[i][2]/signal_yield_2017[i][2]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2017[i][3], stat_yield_2017[i][3], stat_yield_2017[i][3]/signal_yield_2017[i][3]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2017[i][4], stat_yield_2017[i][4], stat_yield_2017[i][4]/signal_yield_2017[i][4]) << "\\\\" << std::endl;

    file_stat << " & " << " 2018 & " << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2018[i][0], stat_yield_2018[i][0], stat_yield_2018[i][0]/signal_yield_2018[i][0]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2018[i][1], stat_yield_2018[i][1], stat_yield_2018[i][1]/signal_yield_2018[i][1]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2018[i][2], stat_yield_2018[i][2], stat_yield_2018[i][2]/signal_yield_2018[i][2]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2018[i][3], stat_yield_2018[i][3], stat_yield_2018[i][3]/signal_yield_2018[i][3]) << Form("%.0lf $\\pm$ %.0lf (%.3lf) & ", signal_yield_2018[i][4], stat_yield_2018[i][4], stat_yield_2018[i][4]/signal_yield_2018[i][4]) << "\\\\ \\hline" << std::endl;
  }

  file_stat << "\\end{tabular}" << std::endl;
  file_stat << "\\end{table}" << std::endl;
  file_stat.close();

  //SYST TABLE
  std::string filename_syst = "/afs/cern.ch/user/m/mfaria/public/UML-fit/Systematics/tables/yield_syst_table.tex";
  file_syst.open(filename_syst);

  if(!file_syst.is_open())
      cout << "Output file_syst not opened!";

  file_syst << "\\begin{table}" << std::endl;
  file_syst << "\\caption{Yields with respective systematic uncertainties for each fit variation.}" << std::endl;
  file_syst << "\\scriptsize" << std::endl;
  file_syst << "\\centering" << std::endl;
  file_syst << "\\begin{tabular}{|c|c|c|c|c|c|c|c|}" << std::endl;
  file_syst << "\\hline" << std::endl;

  std::vector<std::string> col_name_syst = {"$q^2$ bin", "Year", "Exp($5.1<m<5.6$)", "Exp($5.0<m<5.5$)", "WT fixed/Exp+erf", "Scale factor", "Total"};
  for(int c = 0; c < 6; c++){
    file_syst << col_name_syst[c] << " & ";
  }
  file_syst << col_name_syst[6] << " \\\\ " << std::endl;
  file_syst << "\\hline" << std::endl;

  for(int i = 0; i < n_q2Bin; i++){

    file_syst << " & " << " 2016 & " << Form("%.2lf & ", (syst_yield_2016[i][1]/signal_yield_2016[i][1])*100.) << Form("%.2lf & ", (syst_yield_2016[i][2]/signal_yield_2016[i][2])*100.) << Form("%.2lf & ", (syst_yield_2016[i][3]/signal_yield_2016[i][3])*100.) << Form("%.2lf & ", (syst_yield_2016[i][4]/signal_yield_2016[i][4])*100.) << Form("%.2lf ", (yields_syst_2016[i]/yields_2016[i])*100.) << "\\\\" << std::endl;

    file_syst << Form("%i & ",i) << " 2017 & " << Form("%.2lf & ", (syst_yield_2017[i][1]/signal_yield_2017[i][1])*100.) << Form("%.2lf & ", (syst_yield_2017[i][2]/signal_yield_2017[i][2])*100.) << Form("%.2lf & ", (syst_yield_2017[i][3]/signal_yield_2017[i][3])*100.) << Form("%.2lf & ", (syst_yield_2017[i][4]/signal_yield_2017[i][4])*100.) << Form("%.2lf ", (yields_syst_2017[i]/yields_2017[i])*100.) << "\\\\" << std::endl;

    file_syst << " & " << " 2018 & " << Form("%.2lf & ", (syst_yield_2018[i][1]/signal_yield_2018[i][1])*100.) << Form("%.2lf & ", (syst_yield_2018[i][2]/signal_yield_2018[i][2])*100.) << Form("%.2lf & ", (syst_yield_2018[i][3]/signal_yield_2018[i][3])*100.) << Form("%.2lf & ", (syst_yield_2018[i][4]/signal_yield_2018[i][4])*100.) << Form("%.2lf ", (yields_syst_2018[i]/yields_2018[i])*100.) << "\\\\ \\hline" << std::endl;
  }

  file_syst << "\\end{tabular}" << std::endl;
  file_syst << "\\end{table}" << std::endl;
  file_syst.close();

}

double getMax(double list[5], int n_pdf){

  double max = list[0];

  for(int i = 0; i <= n_pdf; i++){
    if(max < list[i]) {
      max = list[i];
    }   
  }
  return max; 
}
