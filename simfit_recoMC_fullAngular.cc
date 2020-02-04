#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMath.h>
#include <TH3D.h>

#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
 #include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>

#include "PdfSigAng.h"
#include "ParBound.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* c [4*nBins];

void simfit_recoMC_fullAngularBin(int q2Bin, int parity, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data)
{

  string shortString = Form("b%ip%i",q2Bin,parity);
  cout<<"Conf: "<<shortString<<endl;

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
  RooRealVar* P4p   = new RooRealVar("P4p","P'_{4}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P5p   = new RooRealVar("P5p","P'_{5}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P6p   = new RooRealVar("P6p","P'_{6}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* P8p   = new RooRealVar("P8p","P'_{8}",0,-1*sqrt(2),sqrt(2));
  RooRealVar* mFrac = new RooRealVar("mFrac","mistag fraction",1, 0, 2);
  mFrac->setConstant();

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
      dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
        ->reduce( RooArgSet(reco_vars), Form("rand < %f", scale_to_data[years[iy]] )) ;
      dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
        ->reduce( RooArgSet(reco_vars), Form("rand < %f", scale_to_data[years[iy]] )) ;
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
    string filename = "../eff-KDE/files_forSimFit/KDEeff_b";
    filename = filename + Form((parity==0 ? "%i_ev_%i.root" : "%i_od_%i.root"),q2Bin,years[iy]);
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
  P4p->setVal(0);
  P5p->setVal(0);
  P6p->setVal(0);
  P8p->setVal(0);

  // Construct combined dataset in (x,sample)
  RooDataSet combData ("combData", "combined data", 
                         vars,
                         Index(sample),
                         Import(map)); 

  // Define the PDF defining the physical region, used as exexternal contraint in the fit
  RooAbsPdf* PDF_phys_bound = new ParBound(("PDF_phys_bound_"+shortString).c_str(),"PDF_phys_bound",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  ((ParBound*)PDF_phys_bound)->Q2Bin=q2Bin;

  RooAbsReal* nll = simPdf->createNLL(combData,
                                      RooFit::Extended(kFALSE),
//                                       RooFit::ExternalConstraints(*PDF_phys_bound),
                                      RooFit::NumCPU(1)
                                      );
    
  RooMinimizer m(*nll) ;
  m.optimizeConst (kTRUE); // do not recalculate constant terms
  m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
  m.setVerbose(kTRUE);
//  Minuit2.setEps(1e-16) ;
  m.setMinimizerType("Minuit2");
  m.setStrategy(0);
//   m.setEvalErrorWall(false);
  m.migrad() ;
  m.hesse() ;   
  std::cout << std::endl;
  std::cout << "######################### now strategy 2 #########################"<< std::endl;
  m.setStrategy(2);
  m.migrad() ;
  m.hesse() ;   
  m.minos() ;
    
  RooFitResult* fitResult = m.save() ; 
  fitResult->Print("v");
  

  // perform fit in two steps:
  // first with strategy=0 and no MINOS, to get close to the likilihood maximum
//   simPdf->fitTo( combData,
// 		 ExternalConstraints(*PDF_phys_bound),
//                  Minimizer("Minuit2","migrad"), 
//                  Extended(false), 
//                  Timer(true),
//                  NumCPU(1),
//                  Hesse(true),
//                  Strategy(0),
//                  Offset(true));
//   // second with full accuracy (strategy=2) and MINOS, to get the best result
//   RooFitResult* fitResult = simPdf->fitTo( combData,
// 					   ExternalConstraints(*PDF_phys_bound),
//                                            Minimizer("Minuit2","migrad"),
//                                            Extended(false), 
//                                            Save(true),
// //                                            Timer(true),
//                                            // NumCPU(1),
//                                            Hesse(true),
//                                            Strategy(2),
//                                            Minos(true),
//                                            Offset(true)
//                                          );
  // The two step fit seemed necessary to get convergence in all q2 bins
  // if one observes that the convergence is obtained just running the second step, it is better to avoid the first step (to gain computing time)
  // Offset(true) parameter make the FCN value to be defined at 0 at the begin of the minimization
  // it is needed to get convergence in all q2 bins, and avoid the "machine precision reached" errors
  // which are due to too high FCN absolute values compared to the minimisation steps

//   fitResult->Print("v");
  ((ParBound*)PDF_phys_bound)->verbose=true;
//  
//  Boundary Checks [0 is ok, ret=final penality, ret_local4="numerically computed" penality]
  cout<<" Boundary Checks [0 is ok] = "<<PDF_phys_bound->getVal()<<endl;

  if (save) {
    // Save fit results in file
    TFile* fout = new TFile(("simFitResultsExtConstr/simFitResult_recoMC_fullAngular" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"UPDATE");
    fitResult->Write(("simFitResult_"+shortString).c_str(),TObject::kWriteDelete);
    fout->Close(); 
  }

  RooPlot* frameFl = Fl->frame(Title("-log(L) scan vs Fl")) ;
  nll->plotOn(frameFl,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP1 = P1->frame(Title("-log(L) scan vs P1")) ;
  nll->plotOn(frameP1,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP2 = P2->frame(Title("-log(L) scan vs P2")) ;
  nll->plotOn(frameP2,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP3 = P3->frame(Title("-log(L) scan vs P3")) ;
  nll->plotOn(frameP3,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP4p = P4p->frame(Title("-log(L) scan vs P4p")) ;
  nll->plotOn(frameP4p,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP5p = P5p->frame(Title("-log(L) scan vs P5p")) ;
  nll->plotOn(frameP5p,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP6p = P6p->frame(Title("-log(L) scan vs P6p")) ;
  nll->plotOn(frameP6p,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
  RooPlot* frameP8p = P8p->frame(Title("-log(L) scan vs P8p")) ;
  nll->plotOn(frameP8p,PrintEvalErrors(-1),ShiftToZero(),EvalErrorValue(nll->getVal()+10),LineColor(kRed)) ;
//   frameFl->SetMaximum(15) ;
//   frameFl->SetMinimum(0) ;

  TCanvas* cSara = new TCanvas("nllerrorhandling","nllerrorhandling",1200,900) ;
  cSara->Divide(3,3) ;
  cSara->cd(1) ; gPad->SetLeftMargin(0.15) ; frameFl->GetYaxis()->SetTitleOffset(1.4)  ; frameFl->Draw() ;
  cSara->cd(2) ; gPad->SetLeftMargin(0.15) ; frameP1->GetYaxis()->SetTitleOffset(1.4)  ; frameP1->Draw() ;
  cSara->cd(3) ; gPad->SetLeftMargin(0.15) ; frameP2->GetYaxis()->SetTitleOffset(1.4)  ; frameP2->Draw() ;
  cSara->cd(4) ; gPad->SetLeftMargin(0.15) ; frameP3->GetYaxis()->SetTitleOffset(1.4)  ; frameP3->Draw() ;
  cSara->cd(5) ; gPad->SetLeftMargin(0.15) ; frameP4p->GetYaxis()->SetTitleOffset(1.4) ; frameP4p->Draw() ;
  cSara->cd(6) ; gPad->SetLeftMargin(0.15) ; frameP5p->GetYaxis()->SetTitleOffset(1.4) ; frameP5p->Draw() ;
  cSara->cd(7) ; gPad->SetLeftMargin(0.15) ; frameP6p->GetYaxis()->SetTitleOffset(1.4) ; frameP6p->Draw() ;
  cSara->cd(8) ; gPad->SetLeftMargin(0.15) ; frameP8p->GetYaxis()->SetTitleOffset(1.4) ; frameP8p->Draw() ;
  cSara->SaveAs(("test_nll_Constr_" + all_years + stat + Form("_b%i_Paolo.pdf", q2Bin)).c_str());

  if (!plot) return;

  int confIndex = 2*nBins*parity  + q2Bin;
  string longString  = "Fit to reconstructed events";
  longString = longString + Form(parity==1?" (q2-bin %i even)":" (q2-bin %i odd)",q2Bin);

  // plot fit projections
  c[confIndex] = new TCanvas (("c_"+shortString).c_str(),("Fit to RECO-level MC - "+longString).c_str(),2000,1400);
  c[confIndex]->Divide(3, years.size());
  
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

    c[confIndex]->cd(iy*3+1);
    gPad->SetLeftMargin(0.19); 
    xframe->Draw();
    leg->Draw("same");
    c[confIndex]->cd(iy*3+2);
    gPad->SetLeftMargin(0.19); 
    yframe->Draw();
    leg->Draw("same");
    c[confIndex]->cd(iy*3+3);
    gPad->SetLeftMargin(0.19); 
    zframe->Draw();
    leg->SetTextSize(0.03);
    leg->AddEntry(xframe->findObject(("plData"+year).c_str()),("Post-selection distribution "+year).c_str() ,"lep");
    leg->AddEntry(xframe->findObject(("plPDF"+year ).c_str()),("Decay rate x efficiency "+year).c_str(),"l");
    leg->Draw("same");
  }

  
  c[confIndex]->SaveAs( ("plotSimFit_d/simFitResult_recoMC_fullAngular_" + shortString + "_" + all_years + stat + ".pdf").c_str() );
}


void simfit_recoMC_fullAngularBin1(int q2Bin, int parity, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullAngularBin(q2Bin, parity, plot, save, datalike,  years, scale_to_data);
  else
    simfit_recoMC_fullAngularBin(q2Bin, parity, plot, save, datalike, years, scale_to_data);
}

// void simfit_recoMC_fullAngularBin1(int q2Bin, int parity, int tagFlag, bool plot, bool save, bool datalike, std::vector<int> years, std::map<int,float> scale_to_data)
// {
//   if ( parity==-1 )
//     for (parity=0; parity<2; ++parity)
//       simfit_recoMC_fullAngularBin2(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);
//   else
//     simfit_recoMC_fullAngularBin2(q2Bin, parity, tagFlag, plot, save, datalike, years, scale_to_data);
// }

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively
  // tag format: [0] mistagged
  //             [1] correctly-tagged
  //             [-1] for each tag recursively

  int q2Bin   = -1;
  int parity  = -1; 
//   int tagFlag = -1;

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);
//   if ( argc >= 4 ) tagFlag = atoi(argv[3]);

  bool plot = true;
  bool save = true;
  bool datalike = true;

  if ( argc >= 4 && atoi(argv[3]) == 0 ) plot     = false;
  if ( argc >= 5 && atoi(argv[4]) == 0 ) save     = false;
  if ( argc >= 6 && atoi(argv[5]) == 0 ) datalike = false;

  std::vector<int> years;
  if ( argc >= 7 && atoi(argv[6]) != 0 ) years.push_back(atoi(argv[6]));
  else {
    cout << " no specific years selected, using default: 2016";
    years.push_back(2016);
  }
  if ( argc >= 8  && atoi(argv[7]) != 0 ) years.push_back(atoi(argv[7]));
  if ( argc >= 9  && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;
//   if ( tagFlag < -1 || tagFlag > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;
//   if ( tagFlag==-1 ) cout << "Running both the tag conditions"  << endl;
  if ( datalike )    cout << "Considering data-like statistics" << endl;

  std::map<int,float> scale_to_data;
  scale_to_data.insert(std::make_pair(2016, 0.01*2/4)); // *2 since we are using only odd/even events 
  scale_to_data.insert(std::make_pair(2017, 0.01*2/4));
  scale_to_data.insert(std::make_pair(2018, 0.015*2/4));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_recoMC_fullAngularBin1(q2Bin, parity, plot, save, datalike, years, scale_to_data);
  else
    simfit_recoMC_fullAngularBin1(q2Bin, parity, plot, save, datalike, years, scale_to_data);

  return 0;

}
