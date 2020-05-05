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
#include "PdfSigAng.h"
#include "Penalty.h"
#include "BoundCheck.h"

using namespace RooFit;
using namespace std;

int nBins = 500;
int nPlotBins = 40;
int nPlotBinsZoom = 20;

// int nSigZoom = 60;

TCanvas* c;
TCanvas* cZ;
TCanvas* c2d1;
TCanvas* c2d2;

void plotBound(TH1* h2Bound, TVirtualPad* pad, int color, float size);

void plotMultiResultsBin(int q2Bin, int parity, bool do2016, bool do2017, bool do2018, bool isGEN, int boundAnalysis)
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
  if (isGEN) wsName = "wsMulti_" + shortString + "_s" + yearString + "_n200_pow1.0";
  else wsName = "wsMulti_" + shortString + "_s100_pow1.0";
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

  // import fit result of high-stat MC sample, and data+PDF to create NLL
  string fileNameRes = "";
  if (isGEN) fileNameRes = "fitResult_genMC_penalty";
  else fileNameRes = Form(("simFitResult_recoMC_fullAngular"+yearString+"_MCStat_b%i").c_str(),q2Bin);
  TFile* fin_res = TFile::Open( ("simFitResults/"+fileNameRes+".root").c_str() );
  if ( !fin_res || !fin_res->IsOpen() ) {
    cout<<"File not found: simFitResults/"<<fileName<<".root"<<endl;
    return;
  }
  string wsNameRes = "ws_" + shortString + "_s0_pow1.0";
  RooWorkspace* wspRes = (RooWorkspace*)fin_res->Get(wsNameRes.c_str());
  if ( !wspRes || wspRes->IsZombie() ) {
    cout<<"Workspace "<<wsNameRes<<" not found in file: simFitResults/"<<fileNameRes<<".root"<<endl;
    return;
  }
  RooDataSet* data = (RooDataSet*)wspRes->data("data");
  if ( !data || data->IsZombie() ) {
    cout<<"Dataset not found in workspace "<<wsNameRes<<" in file: simFitResults/"<<fileNameRes<<".root"<<endl;
    return;
  }
  RooAbsPdf* pdf = wspRes->pdf("pdf");
  if ( !pdf || pdf->IsZombie() ) {
    cout<<"Free pdf not found in workspace "<<wsNameRes<<" in file: simFitResults/"<<fileNameRes<<".root"<<endl;
    return;
  }
  RooAbsReal* nll = pdf->createNLL(*data,Extended(kFALSE),NumCPU(1));
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

  // when boundary analysis is selected, single component of the boundary are added in the 2D plots
  BoundCheck* boundary4  = 0;
  BoundCheck* boundary15 = 0;
  vector <BoundCheck*> boundary15v (0);
  if ( boundAnalysis>0 ) {
    // CTL4 component
    boundary4 = new BoundCheck("bound4","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,true,-1);
    if ( boundAnalysis>1 && boundAnalysis<=100 ) {
      // here a set of CTL15 components are created for individual phi values
      int nScan = 100 / boundAnalysis;
      for (int iPhi=nScan; iPhi<=100; iPhi+=nScan)
	boundary15v.push_back(new BoundCheck(Form("bound15%i",iPhi),"Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false,0.01*iPhi*TMath::Pi()));
    } else 
      // here the global CTL15 component is created
      boundary15 = new BoundCheck("bound15","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,false);
  }

  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  for (int iPar = 0; iPar < pars.getSize(); ++iPar)
    ((RooRealVar*)pars.at(iPar))->setVal(((RooRealVar*)fitResult->floatParsFinal().find(pars.at(iPar)->GetName()))->getValV());

  if (boundAnalysis==0) {
    c  = new TCanvas (("c_"+shortString).c_str(),("c_"+shortString).c_str(),1800,1800);
    cZ = new TCanvas (("cZ_"+shortString).c_str(),("cZ_"+shortString).c_str(),1800,1800);
    c ->Divide(3,3);
    cZ->Divide(3,3);
  }
  c2d1 = new TCanvas (("c2d1_"+shortString).c_str(),("c2d1_"+shortString).c_str(),3000,1800);
  c2d2 = new TCanvas (("c2d2_"+shortString).c_str(),("c2d2_"+shortString).c_str(),3000,1800);
  c2d1->Divide(5,3);
  c2d2->Divide(5,3);

  RooPlot* frame [8];

  vector <double> lowRange  (0);
  vector <double> highRange (0);

  for (int iPar = 0; iPar < pars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)pars.at(iPar);

    auto hRes = subResults->createHistogram("hRes",*par,Binning(nPlotBins));
    auto hResGood = subPosConv->createHistogram("hResGood",*par,Binning(nPlotBins));

    hResGood->SetLineColor(8);
    
    if (boundAnalysis==0) {
      c->cd(iPar+1);
      hRes->Draw();
      hResGood->Draw("same");
    }

    // compute range for zommed plots
    int iBin = 1;
    for (iBin=1; hRes->GetBinContent(iBin)==0; ++iBin);
    lowRange.push_back (max(par->getMin(),
    			    hRes->GetBinCenter(iBin) - 0.5*hRes->GetBinWidth(iBin) - 0.25 * ( hRes->GetMean() - hRes->GetBinCenter(iBin) )) );
    for (iBin=hRes->GetNbinsX(); hRes->GetBinContent(iBin)==0; --iBin);
    highRange.push_back(min(par->getMax(),
    			    hRes->GetBinCenter(iBin) + 0.5*hRes->GetBinWidth(iBin) - 0.25 * ( hRes->GetMean() - hRes->GetBinCenter(iBin) )) );

    if (boundAnalysis==0) {
      frame[iPar] = par->frame(Name(Form("f%s",par->GetName())),Title(Form("Result distribution of %s",par->GetTitle())),
			       Range(lowRange[iPar],highRange[iPar])) ;
      
      // NLL function of high-stat fit
      nll->plotOn(frame[iPar],
		  PrintEvalErrors(-1),
		  ShiftToZero(),
		  EvalErrorValue(nll->getVal()+1),
		  LineColor(kRed),
		  LineWidth(2)) ;
      
      // compute y range needed to contain the likelihood, expecting a shape similar to a parabolic function
      double parSigma = ((RooRealVar*)fitResult->floatParsFinal().find(pars.at(iPar)->GetName()))->getError();
      double hZoom = 0.5 * pow( (highRange[iPar]-lowRange[iPar])/2/parSigma, 2 );
      
      // boundary represented by shading the non-physical region
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
      
      // 1D histos with fit results (hzRes is visible only if there are not-good fit results)
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
      
    }

    // Loop over a second parameter to produce 2D plots
    for (int iPar2 = 0; iPar2 < iPar; ++iPar2) {

      int indx = iPar2 + iPar*(iPar-1)/2;
      RooRealVar* par2 = (RooRealVar*)pars.at(iPar2);

      // Histogram with the results
      TH1* h2ResGood = subPosConv->createHistogram("h2ResGood",*par2,Binning(50,lowRange[iPar2],highRange[iPar2]),
      						   YVar(*par,Binning(50,lowRange[iPar],highRange[iPar])));
      h2ResGood->SetMaximum(7);

      // Boundaries are represented as 2D histo, then will be plot as TGraph with their countours by the plotBound function
      TH1* h2Bound = boundary->createHistogram("h2Bound",*par2,Binning(nBins,lowRange[iPar2],highRange[iPar2]),
					       YVar(*par,Binning(nBins,lowRange[iPar],highRange[iPar])),
					       Extended(false),Scaling(false));

      TH1* h2Bound4 = 0;
      TH1* h2Bound15 = 0;
      vector <TH1*> h2Bound15v (0);
      if ( boundAnalysis>0 ) {
	h2Bound4 = boundary4->createHistogram("h2Bound4",*par2,Binning(nBins,lowRange[iPar2],highRange[iPar2]),
					      YVar(*par,Binning(nBins,lowRange[iPar],highRange[iPar])),
					      Extended(false),Scaling(false));
 	if ( boundAnalysis>1 && boundAnalysis<=100 ) {
	  for (uint iBound=0; iBound<boundary15v.size(); ++iBound)
	    h2Bound15v.push_back(boundary15v[iBound]->createHistogram(Form("h2Bound15%u",iBound),*par2,Binning(nBins,lowRange[iPar2],highRange[iPar2]),
								      YVar(*par,Binning(nBins,lowRange[iPar],highRange[iPar])),
								      Extended(false),Scaling(false)));
	} else
	  h2Bound15 = boundary15->createHistogram("h2Bound15",*par2,Binning(nBins,lowRange[iPar2],highRange[iPar2]),
						  YVar(*par,Binning(nBins,lowRange[iPar],highRange[iPar])),
						  Extended(false),Scaling(false));
      }

      // One-point graph with high-stat fit result
      TGraph* grBestFit = new TGraph(1);
      grBestFit->SetPoint(0,par2->getValV(),par->getValV());
      grBestFit->SetMarkerColor(2);
      grBestFit->SetMarkerStyle(29);
      grBestFit->SetMarkerSize(5);

      // Prepare alternative palette to be used when boundAnalysis>1
      const Int_t Number = 3;
      Double_t Red[Number]    = { 1.00, 0.00, 0.00};
      Double_t Green[Number]  = { 0.00, 1.00, 0.00};
      Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
      Double_t Length[Number] = { 0.00, 0.50, 1.00 };
      Int_t nb=100;

      // Plot 2D objects
      TVirtualPad* pad = 0;
      if (indx<13) pad = c2d1->cd(indx+1);
      else pad = c2d2->cd(indx-12);
      gPad->SetLeftMargin(0.15);

      h2ResGood->Draw("COLZ");
      
      if ( boundAnalysis>0 ) {
	plotBound(h2Bound,pad,2,0.8);
	plotBound(h2Bound4,pad,797,0.5);
 	if ( boundAnalysis>1 && boundAnalysis<=100 ) {
	  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
	  for (uint iHist=0; iHist<h2Bound15v.size(); ++iHist)
	    plotBound(h2Bound15v[iHist],pad,0,0.4);
	} else
	  plotBound(h2Bound15,pad,4,0.4);
      } else
	plotBound(h2Bound,pad,13,0.7);

      grBestFit->Draw("P");

    }

  }

  string plotName = (isGEN?"multiFitPlot_genMC_":"multiFitPlot_recoMC_") + shortString;
  if (do2016) plotName = plotName + "_2016";
  if (do2017) plotName = plotName + "_2017";
  if (do2018) plotName = plotName + "_2018";

  if (boundAnalysis==0) {

    c->Update();
    c->SaveAs( ("plotSimFit_d/"+plotName+".pdf").c_str() );
    
    cZ->Update();
    cZ->SaveAs( ("plotSimFit_d/"+plotName+"_zoom.pdf").c_str() );

    c2d1->Update();
    c2d1->SaveAs( ("plotSimFit_d/"+plotName+"_2D.pdf").c_str() );
 
    c2d2->Update();
    c2d2->SaveAs( ("plotSimFit_d/"+plotName+"_2Dbis.pdf").c_str() );
 
  } else {

  c2d1->Update();
  c2d1->SaveAs( ("plotSimFit_d/"+plotName+Form("_boundAnalysis%i.pdf",boundAnalysis)).c_str() );

  c2d2->Update();
  c2d2->SaveAs( ("plotSimFit_d/"+plotName+Form("_boundAnalysis%i_bis.pdf",boundAnalysis)).c_str() );

  }

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

  int boundAnalysis = 0;

  if ( argc > 7 ) boundAnalysis = atoi(argv[7]);

  if ( q2Bin  < -1 || q2Bin  > 7 ) return 1;
  if ( parity <  0 || parity > 1 ) return 1;

  if ( q2Bin==-1 ) {
    cout<<"Running all the q2 bins"<<endl;
    for (q2Bin=0; q2Bin<8; ++q2Bin) {
      if ( q2Bin==4 || q2Bin==6 ) continue;
      plotMultiResultsBin(q2Bin, parity, do2016, do2017, do2018, doGen, boundAnalysis);
    }
  } else
    plotMultiResultsBin(q2Bin, parity, do2016, do2017, do2018, doGen, boundAnalysis);

  return 0;

}
