#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>
#include <TStopwatch.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>
#include <RooMinimizer.h>

#include "PdfSigAng.h"
#include "ParBound.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cnll;
TCanvas* c;

void simfit_recoMC_fullAngularBin(int q2Bin, int parity, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data, double bound_power, double bound_shift1, double bound_shift5, double bound_coeff, double bound_int, int nSample)
{

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;
  cout<<"Bound: pow="<<bound_power<<" shift1="<<bound_shift1<<" shift5="<<bound_shift5<<" coeff="<<bound_coeff<<" integral="<<bound_int<<endl;

  // Load variables and dataset
  // importing the complementary dataset, to fit with statistically uncorrelated efficiency

  string effCString = Form("effCHist_b%ip%i",q2Bin,parity);
  string effWString = Form("effWHist_b%ip%i",q2Bin,parity);
  string intCHistString = "MCint_"+shortString + "t1";
  string intWHistString = "MCint_"+shortString + "t0";
  string all_years = "";
  string year = ""; 
  string stat = datalike ? "_dataStat":"_MCStat";

  std::vector<TFile*> fin_data, fin_eff;
  std::vector<RooWorkspace*> wsp;
  std::vector<RooDataSet*> data;
  std::vector<RooAbsReal*> effC, effW;
  std::vector<TH3D*> effCHist, effWHist;
  std::vector<TH1D*> intCHist, intWHist;
  std::vector< std::vector<double> > intCVec(years.size(), std::vector<double>(0));
  std::vector< std::vector<double> > intWVec(years.size(), std::vector<double>(0));
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular (0);

  //// from https://root-forum.cern.ch/t/combining-roodatasets-using-std-map-in-pyroot/16471/20
  gInterpreter->GenerateDictionary("std::pair<std::string, RooDataSet*>", "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::map<std::string, RooDataSet*>",  "map;string;RooDataSet.h");
  gInterpreter->GenerateDictionary("std::pair<std::map<string,RooDataSet*>::iterator, bool>", "map;string;RooDataSet.h");
  std::map<std::string, RooDataSet*> map;

  RooRealVar* ctK = new RooRealVar("ctK", "ctK", -1  , 1  );
  RooRealVar* ctL = new RooRealVar("ctL", "ctL", -1  , 1  );
  RooRealVar* phi = new RooRealVar("phi", "phi", -3.14159, 3.14159  );
  RooArgList vars (* ctK,* ctL,* phi);
  RooRealVar* rand = new RooRealVar("rand", "rand", 0,1);
  RooArgSet reco_vars (*ctK, *ctL, *phi, *rand);

  // define angular parameters with ranges from positiveness requirements on the decay rate
  RooRealVar* Fl    = new RooRealVar("Fl","F_{L}",0.5,0,1);
  RooRealVar* P1    = new RooRealVar("P1","P_{1}",0,-1,1);
  RooRealVar* P2    = new RooRealVar("P2","P_{2}",0,-0.5,0.5);
  RooRealVar* P3    = new RooRealVar("P3","P_{3}",0,-0.5,0.5);
  RooRealVar* R45   = new RooRealVar("R45","R45",0.1,0,1);
  RooRealVar* R68   = new RooRealVar("R68","R68",0.1,0,1);
  RooRealVar* phi45 = new RooRealVar("phi45","phi45",0.5*TMath::Pi(),-0.5*TMath::Pi(),1.5*TMath::Pi());
  RooRealVar* phi68 = new RooRealVar("phi68","phi68",0.5*TMath::Pi(),-0.5*TMath::Pi(),1.5*TMath::Pi());
  RooRealVar* mFrac = new RooRealVar("mFrac","mistag fraction",1, 0.95, 1.05);
  mFrac->setConstant();

  RooFormulaVar* P4p = new RooFormulaVar("P4p","P'_{4}","R45*cos(phi45)*sqrt(0.5+P2)+R68*cos(phi68)*sqrt(0.5-P2)",RooArgList(*R45,*R68,*phi45,*phi68,*P2));
  RooFormulaVar* P5p = new RooFormulaVar("P5p","P'_{5}","R45*cos(phi45)*sqrt(0.5+P2)-R68*cos(phi68)*sqrt(0.5-P2)",RooArgList(*R45,*R68,*phi45,*phi68,*P2));
  RooFormulaVar* P6p = new RooFormulaVar("P6p","P'_{6}","R68*sin(phi68)*sqrt(0.5-P2)-R45*sin(phi45)*sqrt(0.5+P2)",RooArgList(*R45,*R68,*phi45,*phi68,*P2));
  RooFormulaVar* P8p = new RooFormulaVar("P8p","P'_{8}","R68*sin(phi68)*sqrt(0.5-P2)+R45*sin(phi45)*sqrt(0.5+P2)",RooArgList(*R45,*R68,*phi45,*phi68,*P2));

  RooCategory sample ("sample", "sample");
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    sample.defineType(("data"+year).c_str());
    all_years += year;
  }
  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);

  // loop on the various datasets
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
    string filename_data = Form("/eos/cms/store/user/fiorendi/p5prime/effKDE/%i/lmnr/recoMCDataset_b%i_%i.root", years[iy], q2Bin, years[iy]); 

    // import data (or MC as data proxy)
    fin_data.push_back( TFile::Open( filename_data.c_str() ) );
    if ( !fin_data[iy] || !fin_data[iy]->IsOpen() ) {
      cout << "File not found: " << filename_data << endl;
      return;
    }
  
    wsp.push_back( (RooWorkspace*)fin_data[iy]->Get(Form("ws_b%ip%i", q2Bin, 1-parity ) ) );
    if ( !wsp[iy] || wsp[iy]->IsZombie() ) {
      cout<<"Workspace not found in file: "<<filename_data<<endl;
      return;
    }
  
    // create roodataset (in case data-like option is selected, only import the correct % of data)
    RooDataSet* dataCT, *dataWT;
    if (datalike){
      if (nSample*scale_to_data[years[iy]]>1) {
	cout<<"ERROR: number of sample requested is too high, not enough MC stat in "<<years[iy]<<endl;
	return;
      }
      dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
        ->reduce( RooArgSet(reco_vars), Form("rand > %1.6f && rand < %1.6f", (nSample-1)*scale_to_data[years[iy]], nSample*scale_to_data[years[iy]] )) ;
      dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
        ->reduce( RooArgSet(reco_vars), Form("rand > %1.6f && rand < %1.6f", (nSample-1)*scale_to_data[years[iy]], nSample*scale_to_data[years[iy]] )) ;
    }
    else{
      dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
    }

    RooDataSet* datatmp = new RooDataSet(*dataCT,("data_"+shortString).c_str());
    datatmp->append(*dataWT);
    data.push_back( datatmp) ;
    if ( !data[iy] || data[iy]->IsZombie() ) {
      cout<<"Dataset not found in file: "<<filename_data<<endl;
      return;
    }

    // import KDE efficiency histograms and partial integral histograms
    string filename = "/eos/cms/store/user/fiorendi/p5prime/effKDE/";
    filename = filename + Form((parity==0 ? "%i/lmnr/KDEeff_b%i_ev_%i.root" : "%i/lmnr/KDEeff_b%i_od_%i.root"),years[iy],q2Bin,years[iy]);
    fin_eff.push_back( new TFile( filename.c_str(), "READ" ));
    if ( !fin_eff[iy] || !fin_eff[iy]->IsOpen() ) {
      cout<<"File not found: "<<filename<<endl;
      return;
    }

    effCHist.push_back( (TH3D*)fin_eff[iy]->Get(effCString.c_str()));
    effWHist.push_back( (TH3D*)fin_eff[iy]->Get(effWString.c_str()));
    if ( !effCHist[iy] || effCHist[iy]->IsZombie() || !effWHist[iy] || effWHist[iy]->IsZombie() ) {
      cout<<"Efficiency histogram "<< effCString <<" or " << effWString << " not found in file: "<< filename <<endl;
      return;
    }

    // create efficiency functions
    RooDataHist* effCData = new RooDataHist(("effCData_"+shortString+"_"+year).c_str(),"effCData",vars,effCHist[iy]);
    RooDataHist* effWData = new RooDataHist(("effWData_"+shortString+"_"+year).c_str(),"effWData",vars,effWHist[iy]);
    effC.push_back( new RooHistFunc(("effC_"+shortString+"_"+year).c_str(),
                                    ("effC"+year).c_str() ,
                                    vars,
                                    *effCData,
                                    1));
    effW.push_back( new RooHistFunc(("effW_"+shortString+"_"+year).c_str(),
                                    ("effW"+year).c_str() ,
                                    vars,
                                    *effWData,
                                    1));

    // import precomputed integrals and fill a std::vector
    intCHist.push_back( (TH1D*)fin_eff[iy]->Get(intCHistString.c_str()));
    intWHist.push_back( (TH1D*)fin_eff[iy]->Get(intWHistString.c_str()));
    intCVec.push_back (vector<double> (0));
    intWVec.push_back (vector<double> (0));
    if ( !intCHist[iy] || intCHist[iy]->IsZombie() || !intWHist[iy] || intWHist[iy]->IsZombie() ) {
      cout << "Integral histogram " << intCHistString <<" or " << intWHistString << " not found in file: "<< filename << endl << "Abort" << endl;
      return;
    } else if ( strcmp( intCHist[iy]->GetTitle(), effCHist[iy]->GetTitle() ) || strcmp( intWHist[iy]->GetTitle(), effWHist[iy]->GetTitle() )) {
    // if the eff_config tag is different between efficiency and precomputed-integral means that they are inconsistent
      cout << "Integral histograms are incoherent with efficiency in file: " << filename << endl;
      cout << "Efficiency (CT) conf: " << effCHist[iy]->GetTitle() <<endl;
      cout << "Integral (CT) conf: "   << intCHist[iy]->GetTitle() <<endl;
      cout << "Efficiency (WT) conf: " << effWHist[iy]->GetTitle() <<endl;
      cout << "Integral (WT) conf: "   << intWHist[iy]->GetTitle() <<endl;
      cout << "Abort"<<endl;
      return;
    } 
    else {
      for (int i=1; i<=intCHist[iy]->GetNbinsX(); ++i) {
        intCVec[iy].push_back(intCHist[iy]->GetBinContent(i));
      }
      for (int i=1; i<=intWHist[iy]->GetNbinsX(); ++i) {
        intWVec[iy].push_back(intWHist[iy]->GetBinContent(i));
      }
    }

    // define angular PDF for signal and only one tag component, using the custom class
    // efficiency function and integral values are passed as arguments
    PDF_sig_ang_fullAngular.push_back( new PdfSigAng(("PDF_sig_ang_fullAngular_"+shortString+"_"+year).c_str(),
                                                     ("PDF_sig_ang_fullAngular_"+year).c_str(),
      		                                     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
      		                                     *effC[iy], *effW[iy], intCVec[iy],intWVec[iy]));
    
    // insert sample in the category map, to be imported in the combined dataset
    map.insert(map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year).c_str(), data[iy]));
    // associate model with the data
    simPdf->addPdf(*PDF_sig_ang_fullAngular[iy], ("data"+year).c_str());
  }

  for(auto it = map.cbegin(); it != map.cend(); ++it)
    std::cout << "dataset: " << it->first << ", with n entries: " << it->second->sumEntries() << "\n";

  // to start the fit, parameters are restored to the center of the parameter space
  Fl ->setVal(0.5);
  P1 ->setVal(0);
  P2 ->setVal(0);
  P3 ->setVal(0);
  R45->setVal(0.1);
  R68->setVal(0.1);
  phi45->setVal(0.5*TMath::Pi());
  phi68->setVal(0.5*TMath::Pi());

  // Construct combined dataset in (x,sample)
  RooDataSet combData ("combData", "combined data", 
                         vars,
                         Index(sample),
                         Import(map)); 

  // Define the PDF defining the physical region, used as exexternal contraint in the fit
  RooAbsPdf* PDF_phys_bound = new ParBound(("PDF_phys_bound_"+shortString).c_str(),"PDF_phys_bound",*P1,*P2,*P3,*R45,*R68,*phi45,*phi68,
					   bound_power,bound_shift1,bound_shift5,bound_coeff,bound_int,true);

  // perform fit in two steps:
  // first with strategy=0 and no MINOS, to get close to the likilihood maximum
  RooAbsReal* nll = simPdf->createNLL(combData,
                                      RooFit::Extended(kFALSE),
                                      RooFit::ExternalConstraints(*PDF_phys_bound),
                                      RooFit::NumCPU(1)
                                      );
    
  RooMinimizer m (*nll) ;
  m.optimizeConst(kTRUE);
  m.setOffsetting(kTRUE);
  m.setVerbose(kTRUE);
  m.setMinimizerType("Minuit2");

  m.setStrategy(0);
  m.migrad() ;
  m.hesse() ;   

  // second with full accuracy (strategy=2) and MINOS, to get the best result
  m.setStrategy(2);
  m.migrad() ;
  m.hesse() ;   
  // m.minos() ;

  // The two step fit seemed necessary to get convergence in all q2 bins
  // if one observes that the convergence is obtained just running the second step, it is better to avoid the first step (to gain computing time)
  // Offset(true) parameter make the FCN value to be defined at 0 at the begin of the minimization
  // it is needed to get convergence in all q2 bins, and avoid the "machine precision reached" errors
  // which are due to too high FCN absolute values compared to the minimisation steps

  RooFitResult* fitResult = m.save(("Fit results "+shortString).c_str(),("Fit results "+shortString).c_str());

  RooPlot* frameFl    = Fl   ->frame(Title("-log(L) scan vs Fl")) ;
  RooPlot* frameP1    = P1   ->frame(Title("-log(L) scan vs P1")) ;
  RooPlot* frameP2    = P2   ->frame(Title("-log(L) scan vs P2")) ;
  RooPlot* frameP3    = P3   ->frame(Title("-log(L) scan vs P3")) ;
  RooPlot* frameR45   = R45  ->frame(Title("-log(L) scan vs R45")) ;
  RooPlot* frameR68   = R68  ->frame(Title("-log(L) scan vs R68")) ;
  RooPlot* framephi45 = phi45->frame(Title("-log(L) scan vs phi45")) ;
  RooPlot* framephi68 = phi68->frame(Title("-log(L) scan vs phi68")) ;
  nll->plotOn(frameFl   ,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(frameP1   ,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(frameP2   ,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(frameP3   ,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(frameR45  ,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(frameR68  ,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(framephi45,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  nll->plotOn(framephi68,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;

  PDF_phys_bound->plotOn(frameP1   ,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;
  PDF_phys_bound->plotOn(frameP2   ,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;
  PDF_phys_bound->plotOn(frameP3   ,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;
  PDF_phys_bound->plotOn(frameR45  ,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;
  PDF_phys_bound->plotOn(frameR68  ,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;
  PDF_phys_bound->plotOn(framephi45,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;
  PDF_phys_bound->plotOn(framephi68,LineColor(kBlue),Normalization(1,RooAbsReal::Raw)) ;

  cnll = new TCanvas (("cnll_"+shortString).c_str(),("cnll_"+shortString).c_str(),1800,1800);
  cnll->Divide(3,3);
  cnll->cd(1); frameFl   ->Draw();
  cnll->cd(2); frameP1   ->Draw();
  cnll->cd(3); frameP2   ->Draw();
  cnll->cd(4); frameP3   ->Draw();
  cnll->cd(5); frameR45  ->Draw();
  cnll->cd(6); frameR68  ->Draw();
  cnll->cd(7); framephi45->Draw();
  cnll->cd(8); framephi68->Draw();
  cnll->SaveAs( ("plotSimFit_d/NLL_scan_" + shortString + "_" + all_years + Form("_drawBound_%.0f_%.0f_%.0f_%.0f.pdf",bound_power,bound_shift1,bound_shift5,bound_coeff)).c_str() );

  frameFl   ->SetMaximum(2.0);
  frameP1   ->SetMaximum(2.0);
  frameP2   ->SetMaximum(2.0);
  frameP3   ->SetMaximum(2.0);
  frameR45  ->SetMaximum(2.0);
  frameR68  ->SetMaximum(2.0);
  framephi45->SetMaximum(2.0);
  framephi68->SetMaximum(2.0);
  cnll->Update();
  cnll->SaveAs( ("plotSimFit_d/NLL_scan_" + shortString + "_" + all_years + Form("_drawBound_%.0f_%.0f_%.0f_%.0f_zoom.pdf",bound_power,bound_shift1,bound_shift5,bound_coeff)).c_str() );

  fitResult->Print("v");

  cout<<endl<<"========== RESULT ==========="<<endl;
  double finalContraintVal = PDF_phys_bound->getValV();
  cout<<"Constraint value: "<<finalContraintVal<<endl;
  cout<<"Fit status: "<<fitResult->status()<<" CovMatrix status (3 is good): "<<fitResult->covQual()<<endl;

  if (save) {
    // Save fit results in file
    if (datalike) {
      TFile* fout = new TFile(("simFitResults/simFitResult_recoMC_fullAngular" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"UPDATE");
      fitResult->Write((Form("simFitResult_s%i_",nSample)+shortString).c_str(),TObject::kWriteDelete);
      fout->Close();
    } else {
      TFile* fout = new TFile(("simFitResults/simFitResult_recoMC_fullAngular" + all_years + ".root").c_str(),"UPDATE");
      fitResult->Write(("simFitResult_"+shortString).c_str(),TObject::kWriteDelete);
      fout->Close();
    }
  }

  if (!plot) return;

  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,1400);
  c->Divide(3, years.size());
  
  for (unsigned int iy = 0; iy < years.size(); iy++) {
    year.clear(); year.assign(Form("%i",years[iy]));
  
    RooPlot* xframe = ctK->frame(Title((longString+year).c_str()));
    RooPlot* yframe = ctL->frame(Title((longString+year).c_str()));
    RooPlot* zframe = phi->frame(Title((longString+year).c_str()));
    xframe->GetYaxis()->SetTitleOffset(1.8);
    yframe->GetYaxis()->SetTitleOffset(1.8);
    zframe->GetYaxis()->SetTitleOffset(1.8);
    xframe->SetMaximum(xframe->GetMaximum()*1.15);
    yframe->SetMaximum(yframe->GetMaximum()*1.15);
    zframe->SetMaximum(zframe->GetMaximum()*1.15);
    xframe->SetMinimum(0);
    yframe->SetMinimum(0);
    zframe->SetMinimum(0);
    TLegend* leg = new TLegend (0.25,0.8,0.9,0.9);

    combData.plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year).c_str()), Name(("plData"+year).c_str()));
    combData.plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year).c_str()));
    combData.plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year).c_str()));

    simPdf->plotOn(xframe,Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1),Name(("plPDF"+year).c_str()));
    simPdf->plotOn(yframe,Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1));
    simPdf->plotOn(zframe,Slice(sample, ("data"+year).c_str()), ProjWData(RooArgSet(sample), combData), LineWidth(1));

    c->cd(iy*3+1);
    gPad->SetLeftMargin(0.19); 
    xframe->Draw();
    leg->Draw("same");
    c->cd(iy*3+2);
    gPad->SetLeftMargin(0.19); 
    yframe->Draw();
    leg->Draw("same");
    c->cd(iy*3+3);
    gPad->SetLeftMargin(0.19); 
    zframe->Draw();
    leg->SetTextSize(0.03);
    leg->AddEntry(xframe->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
    leg->AddEntry(xframe->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
    leg->Draw("same");
  }

  c->SaveAs( ("plotSimFit_d/simFitResult_recoMC_fullAngular_" + shortString + "_" + all_years + stat + ".pdf").c_str() );

}


int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = 0;
  int parity  = 1; 

  double bound_power  = 6.0;
  double bound_shift1 = 1000.0;
  double bound_shift5 = 100.0;
  double bound_coeff  = 50.0;
  double bound_int   = 0.0;

  int nSample = 1;

  bool plot = false;
  bool save = false;
  bool datalike = true;

  std::vector<int> years;

  if ( argc >= 2 ) q2Bin = atoi(argv[1]);

  if ( argc >= 3 ) bound_power  = atoi(argv[2]);
  if ( argc >= 4 ) bound_shift1 = atof(argv[3]);
  if ( argc >= 5 ) bound_shift5 = atof(argv[4]);
  if ( argc >= 6 ) bound_coeff  = atof(argv[5]);
  if ( argc >= 7 ) bound_int    = atof(argv[6]);

  if ( argc >= 8 && atoi(argv[7]) != 0 ) years.push_back(atoi(argv[7]));
  else {
    cout << " no specific years selected, using default: 2016";
    years.push_back(2016);
  }
  if ( argc >= 9  && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));
  if ( argc >= 10 && atoi(argv[9]) != 0 ) years.push_back(atoi(argv[9]));

  if ( argc >= 11 ) nSample = atoi(argv[10]);

  if ( argc >= 12 ) parity = atoi(argv[11]);

  if ( argc >= 13 && atoi(argv[12]) == 0 ) datalike = false;
  if ( argc >= 14 && atoi(argv[13]) == 1 ) plot     = true;
  if ( argc >= 15 && atoi(argv[14]) == 1 ) save     = true;

  if ( q2Bin   < 0 || q2Bin   >= nBins ) return 1;
  if ( parity  < 0 || parity  > 1      ) return 1;

  if ( datalike )    cout << "Considering data-like statistics" << endl;

  std::map<int,float> scale_to_data;
  scale_to_data.insert(std::make_pair(2016, 0.01*2 /2.5/2)); // *2 since we are using only odd/even events 
  scale_to_data.insert(std::make_pair(2017, 0.01*2 /2.05/2));
  scale_to_data.insert(std::make_pair(2018, 0.015*2 /2.05/2));

  simfit_recoMC_fullAngularBin(q2Bin, parity, plot, save, datalike, years, scale_to_data, bound_power, bound_shift1, bound_shift5, bound_coeff, bound_int, nSample);

  return 0;

}
