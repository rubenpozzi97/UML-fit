#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TColor.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooMinimizer.h>

#include "DecayRate.h"
#include "DecayRate_Pen.h"
#include "PdfSigAng.h"
#include "PdfSigAng_Pen.h"
#include "Penalty.h"
#include "BoundCheck.h"

using namespace RooFit;
using namespace std;

int nPlotBins = 40;
int nPlotBinsZoom = 20;

// int nSigZoom = 60;

TCanvas* c;
TCanvas* cZ;
TCanvas* c2;

void plotBound(TH1* h2Bound, TVirtualPad* pad, int color, float size);

void plotMultiResultsBin(int q2Bin, int parity, bool do2016, bool do2017, bool do2018, bool isGEN)
{

  gStyle->SetOptStat(0);

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

  // Load results' container
  string yearString = ""; 
  string fileName = "";
  if (isGEN) {
    if (do2016 && do2017 && do2018) yearString = "0";
    else if (do2016 && !do2017 && !do2018) yearString = "1";
    else if (!do2016 && do2017 && !do2018) yearString = "2";
    else if (!do2016 && !do2017 && do2018) yearString = "3";
    else {
      cout<<"Error: GEN fits not executed on two-year samples"<<endl;
      return;
    }
    fileName = "fitResult_genMC_penalty";
  } else {
    if (do2016) yearString = "2016";
    if (do2017) yearString = yearString + "2017";
    if (do2018) yearString = yearString + "2018";
    fileName = Form(("simFitResult_recoMC_fullAngular"+yearString+"_dataStat_b%i").c_str(),q2Bin);
  }

  TFile* fin_data = TFile::Open( ("simFitResults/"+fileName+".root").c_str() );
  if ( !fin_data || !fin_data->IsOpen() ) {
    cout<<"File not found: simFitResults/"<<fileName<<".root"<<endl;
    return;
  }
  string wsName = "";
  if (isGEN) wsName = "wsMulti_" + shortString + "_s" + yearString + (q2Bin==7?"_n100_pow1.0":"_n500_pow1.0");
  else wsName = "ws_" + shortString;
  RooWorkspace* wsp = (RooWorkspace*)fin_data->Get(wsName.c_str());
  if ( !wsp || wsp->IsZombie() ) {
    cout<<"Workspace "<<wsName<<" not found in file: simFitResults/"<<fileName<<".root"<<endl;
    return;
  }
  RooRealVar* Fl  = wsp->var("Fl");
  RooRealVar* P1  = wsp->var("P1");
  RooRealVar* P2  = wsp->var("P2");
  RooRealVar* P3  = wsp->var("P3");
  RooRealVar* P4p = wsp->var("P4p");
  RooRealVar* P5p = wsp->var("P5p");
  RooRealVar* P6p = wsp->var("P6p");
  RooRealVar* P8p = wsp->var("P8p");

  RooDataSet* subResults = (RooDataSet*)wsp->data("subResults");
  if ( !subResults || subResults->IsZombie() ) {
    cout<<"Dataset subResults not found in file: simFitResults/"<<fileName<<".root"<<endl;
    return;
  }
  RooDataSet* subPosConv = 0;
  if (isGEN) subPosConv = (RooDataSet*)subResults->reduce(Name("subPosConv"),Cut("resStatus==0"));
  else subPosConv = (RooDataSet*)wsp->data("subPosConv");
  if ( !subPosConv || subPosConv->IsZombie() ) {
    cout<<"Dataset subPosConv not found in file: simFitResults/"<<fileName<<".root"<<endl;
    return;
  }
  cout<<"Plotting dataset with "<<subPosConv->numEntries()<<" results"<<endl;

  // import NLL and fit result of high-stat MC sample
  string fileNameRes = "";
  if (isGEN) fileNameRes = "fitResult_genMC_penalty";
  else fileNameRes = Form(("simFitResult_recoMC_fullAngular"+yearString+"_MCStat_b%i").c_str(),q2Bin);
  TFile* fin_res = TFile::Open( ("simFitResults/"+fileNameRes+".root").c_str() );
  if ( !fin_res || !fin_res->IsOpen() ) {
    cout<<"File not found: simFitResults/"<<fileName<<".root"<<endl;
    return;
  }
  string wsNameRes = "";
  if (isGEN) wsNameRes = Form("ws_%s_s0_pow1.0",shortString.c_str());
  else wsNameRes = "ws_" + shortString;
  RooWorkspace* wspRes = (RooWorkspace*)fin_res->Get(wsNameRes.c_str());
  if ( !wspRes || wspRes->IsZombie() ) {
    cout<<"Workspace "<<wsNameRes<<" not found in file: simFitResults/"<<fileNameRes<<".root"<<endl;
    return;
  }
  string NllName = "";
  if (isGEN) NllName = "nll_PDF_sig_ang_decayRate_"+shortString+"_data_genDen_"+(parity==1?"ev":"od")+Form("_b%i",q2Bin);
  else NllName = "nll_simPdf_allcombData";
  RooAbsReal* nll = wspRes->function(NllName.c_str());
  if ( !nll || nll->IsZombie() ) {
    cout<<"NLL "<<NllName<<" not found in workspace "<<wsNameRes<<" in file: simFitResults/"<<fileNameRes<<".root"<<endl;
    return;
  }
  RooFitResult* fitResult = 0;
  if (isGEN) {
    fitResult = (RooFitResult*)wspRes->obj("penFitResult");
    if (!fitResult || fitResult->IsZombie()) {
      fitResult = (RooFitResult*)wspRes->obj("fitResult");
      if (!fitResult || fitResult->IsZombie()) {
	cout<<"No fit result found in workspace "<<wsNameRes<<" in file: simFitResults/"<<fileNameRes<<".root"<<endl;
	return;
      }
    }
    if (fitResult->status()!=0 || fitResult->covQual()!=3) {
      cout<<"Not valid fit result in workspace "<<wsNameRes<<" in file: simFitResults/"<<fileNameRes<<".root"<<endl;
      return;
    }
  } else {
    string fitResultName = "simFitResult_" + shortString + "subs0";
    fitResult = (RooFitResult*)fin_res->Get(fitResultName.c_str());
    if (!fitResult || fitResult->IsZombie() || fitResult->status()!=0 || fitResult->covQual()!=3) {
      cout<<"Not valid fit result: "<<fitResultName<<" in file simFitResults/"<<fileNameRes<<".root"<<endl;
      return;
    }
  }

  // Define boundary check (returning 0 in physical region and 1 outside)
  BoundCheck* boundary = new BoundCheck("bound","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  // BoundCheck* boundary4 = new BoundCheck("bound4","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,true,-1);
  // BoundCheck* boundary15 = new BoundCheck("bound15","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false);

  // vector <BoundCheck*> boundary15v (0);
  // for (int iPhi=4; iPhi<=100; iPhi+=4)
  //   boundary15v.push_back(new BoundCheck(Form("bound15%i",iPhi),"Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false,0.01*iPhi*TMath::Pi()));

  // BoundCheck* boundary15_1 = new BoundCheck("bound15_1","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false,0.10*TMath::Pi());
  // BoundCheck* boundary15_2 = new BoundCheck("bound15_2","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false,0.15*TMath::Pi());
  // BoundCheck* boundary15_3 = new BoundCheck("bound15_3","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false,0.85*TMath::Pi());
  // BoundCheck* boundary15_4 = new BoundCheck("bound15_4","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false,0.90*TMath::Pi());

  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  nll->Print();
  for (int iPar = 0; iPar < pars.getSize(); ++iPar)
    ((RooRealVar*)pars.at(iPar))->setVal(((RooRealVar*)fitResult->floatParsFinal().find(pars.at(iPar)->GetName()))->getValV());

  c  = new TCanvas (("c_"+shortString).c_str(),("c_"+shortString).c_str(),1800,1800);
  cZ = new TCanvas (("cZ_"+shortString).c_str(),("cZ_"+shortString).c_str(),1800,1800);
  c2 = new TCanvas (("c2_"+shortString).c_str(),("c2_"+shortString).c_str(),3600,3000);
  c ->Divide(3,3);
  cZ->Divide(3,3);
  c2->Divide(6,5);

  RooPlot* frame [8];
  RooPlot* frame2D [28];

  vector <double> lowRange  (0);
  vector <double> highRange (0);

  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)pars.at(iPar);

    auto hRes = subResults->createHistogram("hRes",*par,Binning(nPlotBins));
    auto hResGood = subPosConv->createHistogram("hResGood",*par,Binning(nPlotBins));

    hResGood->SetLineColor(8);
    
    c->cd(iPar+1);
    hRes->Draw();
    hResGood->Draw("same");

    int iBin = 1;
    for (iBin=1; hRes->GetBinContent(iBin)==0; ++iBin);
    lowRange.push_back (max(par->getMin(),
			    hRes->GetBinCenter(iBin) - 0.5*hRes->GetBinWidth(iBin) - 0.25 * ( hRes->GetMean() - hRes->GetBinCenter(iBin) )) );
    for (iBin=hRes->GetNbinsX(); hRes->GetBinContent(iBin)==0; --iBin);
    highRange.push_back(min(par->getMax(),
			    hRes->GetBinCenter(iBin) + 0.5*hRes->GetBinWidth(iBin) - 0.25 * ( hRes->GetMean() - hRes->GetBinCenter(iBin) )) );

  }

  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {
  cout<<iPar<<endl;
    RooRealVar* par = (RooRealVar*)pars.at(iPar);

    // plots with zoom around the best-fit value
    // double halfRange = ((RooRealVar*)fitResult->floatParsFinal().find(pars.at(iPar)->GetName()))->getError();
    // double lowRange  = max(par->getMin(),par->getValV()-nSigZoom*halfRange);
    // double highRange = min(par->getMax(),par->getValV()+nSigZoom*halfRange);

    frame[iPar] = par->frame(Name(Form("f%s",par->GetName())),Title(Form("Result distribution of %s",par->GetTitle())),
			     Range(lowRange[iPar],highRange[iPar])) ;

    nll->plotOn(frame[iPar],
		PrintEvalErrors(-1),
		ShiftToZero(),
		EvalErrorValue(nll->getVal()+10),
		LineColor(kRed),
		LineWidth(2)) ;
  cout<<iPar<<endl;
    // double hMax = frame[iPar]->GetMaximum();
    // double hZoom = 0.5 * nSigZoom*nSigZoom;
    double parSigma = ((RooRealVar*)fitResult->floatParsFinal().find(pars.at(iPar)->GetName()))->getError();
    double hZoom = 0.5 * pow( (highRange[iPar]-lowRange[iPar])/2/parSigma, 2 );
    cout<<hZoom<<endl;
    if (iPar>0)
      boundary->plotOn(frame[iPar],
		       LineColor(13),
		       FillColor(13),
		       FillStyle(3545),
		       Normalization(1.1*hZoom,RooAbsReal::Raw),
		       DrawOption("LF"),
		       VLines(),
		       LineWidth(2));

    frame[iPar]->SetMaximum(hZoom);

    auto hzRes = subResults->createHistogram("hzRes",*par,Binning(nPlotBinsZoom,lowRange[iPar],highRange[iPar]));
    auto hzResGood = subPosConv->createHistogram("hzResGood",*par,Binning(nPlotBinsZoom,lowRange[iPar],highRange[iPar]));

    double histScaling = 0.8*hZoom/hzRes->GetMaximum();
    hzRes->Scale(histScaling);
    hzResGood->Scale(histScaling);
    hzResGood->SetLineColor(8);
    
    cZ->cd(iPar+1);
    frame[iPar]->Draw();
    hzRes->Draw("same");
    hzResGood->Draw("same");

    for (int iPar2 = 0; iPar2 < iPar; ++iPar2) {

      int indx = iPar2 + iPar*(iPar-1)/2;
      RooRealVar* par2 = (RooRealVar*)pars.at(iPar2);

      // double halfRange2 = ((RooRealVar*)fitResult->floatParsFinal().find(pars.at(iPar2)->GetName()))->getError();
      // double lowRange2  = max(par2->getMin(),par2->getValV()-nSigZoom*halfRange2);
      // double highRange2 = min(par2->getMax(),par2->getValV()+nSigZoom*halfRange2);

      auto h2ResGood = subPosConv->createHistogram("h2ResGood",*par2,Binning(50,lowRange[iPar2],highRange[iPar2]),
      						   YVar(*par,Binning(50,lowRange[iPar],highRange[iPar])));
      h2ResGood->SetMaximum(7);

      // frame2D[indx] = new RooPlot(*par2,*par);

      auto h2Bound = boundary->createHistogram("h2Bound",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
					       YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
					       Extended(false),Scaling(false));
      // auto h2Bound4 = boundary4->createHistogram("h2Bound4",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 						 YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 						 Extended(false),Scaling(false));
      // auto h2Bound15 = boundary15->createHistogram("h2Bound15",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 						   YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 						   Extended(false),Scaling(false));
      // vector <TH1*> h2Bound15v (0);
      // for (uint iBound=0; iBound<boundary15v.size(); ++iBound)
      // 	h2Bound15v.push_back(boundary15v[iBound]->createHistogram(Form("h2Bound15%u",iBound),*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 								  YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 								  Extended(false),Scaling(false)));

      // auto h2Bound15_1 = boundary15_1->createHistogram("h2Bound15_1",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 						       YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 						       Extended(false),Scaling(false));
      // auto h2Bound15_2 = boundary15_2->createHistogram("h2Bound15_2",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 						       YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 						       Extended(false),Scaling(false));
      // auto h2Bound15_3 = boundary15_3->createHistogram("h2Bound15_3",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 						       YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 						       Extended(false),Scaling(false));
      // auto h2Bound15_4 = boundary15_4->createHistogram("h2Bound15_4",*par2,Binning(500,lowRange[iPar2],highRange[iPar2]),
      // 						       YVar(*par,Binning(500,lowRange[iPar],highRange[iPar])),
      // 						       Extended(false),Scaling(false));

      TGraph* grBestFit = new TGraph(1);
      grBestFit->SetPoint(0,par2->getValV(),par->getValV());
      grBestFit->SetMarkerColor(2);
      grBestFit->SetMarkerStyle(29);
      grBestFit->SetMarkerSize(5);

      TVirtualPad* pad = c2->cd(indx+1);
      // frame2D[indx]->Draw();
      h2ResGood->Draw("COLZ");


      // const Int_t Number = 3;
      // Double_t Red[Number]    = { 1.00, 0.00, 0.00};
      // Double_t Green[Number]  = { 0.00, 1.00, 0.00};
      // Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
      // Double_t Length[Number] = { 0.00, 0.50, 1.00 };
      // Int_t nb=50;
      // TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
 

      plotBound(h2Bound,pad,13,0.5);
      // for (uint iHist=0; iHist<h2Bound15v.size(); ++iHist)
      // 	plotBound(h2Bound15v[iHist],pad,0,0.4);
      // plotBound(h2Bound4,pad,797,0.4);
      // plotBound(h2Bound15,pad,393,0.4);

      // plotBound(h2Bound15_1,pad,4,0.4);
      // plotBound(h2Bound15_2,pad,3,0.4);
      // plotBound(h2Bound15_3,pad,6,0.4);
      // plotBound(h2Bound15_4,pad,1,0.4);

      grBestFit->Draw("P");

    }

  }

  string plotName = (isGEN?"multiFitPlot_genMC_":"multiFitPlot_recoMC_") + shortString;
  if (do2016) plotName = plotName + "_2016";
  if (do2017) plotName = plotName + "_2017";
  if (do2018) plotName = plotName + "_2018";

  c->Update();
  c->SaveAs( ("plotSimFit_d/"+plotName+".pdf").c_str() );

  cZ->Update();
  cZ->SaveAs( ("plotSimFit_d/"+plotName+"_zoom.pdf").c_str() );
 
  c2->Update();
  c2->SaveAs( ("plotSimFit_d/"+plotName+"_2D.pdf").c_str() );

}

void plotBound(TH1* h2Bound, TVirtualPad* pad, int color, float size)
{

  double cont[1] = {0.5};
  h2Bound->SetContour(1,cont);

  TCanvas* cDummy = new TCanvas("cDummy","",600,600);
  cDummy->cd();
  h2Bound->Draw("CONT LIST");
  gPad->Update();

  TObjArray* contours = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  pad->cd();
  if (contours->GetSize()>0) {
    TList *list = (TList*)contours->At(0);
    TGraph* gr = (TGraph*)list->First();
    TGraph* gc = 0;
    for (int iGraph=0; iGraph<list->GetSize(); ++iGraph) {
      if (color>0) gr->SetMarkerColor(color);
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(size);
      gc = (TGraph*)gr->Clone();
      if (color>0) gc->Draw("P");
      else gc->Draw("P PMC");
      gr = (TGraph*)list->After(gr);
    }
  }

  delete cDummy;

  return;
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency

  int q2Bin  = -1;
  int parity = 1;

  if ( argc > 1 ) q2Bin  = atoi(argv[1]);
  if ( argc > 2 ) parity = atoi(argv[2]);

  bool do2016 = true;
  bool do2017 = true;
  bool do2018 = true;

  if ( argc > 3 && atoi(argv[3]) == 0 ) do2016 = false;
  if ( argc > 4 && atoi(argv[4]) == 0 ) do2017 = false;
  if ( argc > 5 && atoi(argv[5]) == 0 ) do2018 = false;

  bool doGen = false;

  if ( argc > 6 && atoi(argv[6]) > 0 ) doGen = true;

  if ( q2Bin  < -1 || q2Bin  > 7 ) return 1;
  if ( parity <  0 || parity > 1 ) return 1;

  if ( q2Bin==-1 ) {
    cout<<"Running all the q2 bins"<<endl;
    for (q2Bin=0; q2Bin<8; ++q2Bin) {
      if ( q2Bin==4 || q2Bin==6 ) continue;
      plotMultiResultsBin(q2Bin, parity, do2016, do2017, do2018, doGen);
    }
  } else
    plotMultiResultsBin(q2Bin, parity, do2016, do2017, do2018, doGen);

  return 0;

}
