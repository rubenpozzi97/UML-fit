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
#include <RooSuperCategory.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMinimizer.h>
#include <RooPlot.h>
#include <RooHistFunc.h>
#include <RooDataHist.h>
#include <RooSimultaneous.h>
#include <RooNumIntConfig.h>

#include "PdfSigAng.h"
#include "PdfSigAng_Pen.h"
#include "BoundCheck.h"
#include "Penalty.h"

using namespace RooFit;
using namespace std;

static const int nBins = 9;

TCanvas* cnll;
TCanvas* cZoom;
TCanvas* c [4*nBins];

double fac1 = 10;
double fac4 = 10;
double fac5 = 10;
double base1 = 3;
double base4 = 3;
double base5 = 3;
double maxCoeff = 1e8;

void simfit_recoMC_fullAngularBin(int q2Bin, int parity, bool plot, bool save, int datalike, std::vector<int> years, std::map<int,float> scale_to_data, double power)
{

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

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
  string isample = ""; 
  string stat = datalike > 0 ? "_dataStat":"_MCStat";
  unsigned int nSamples = datalike>1 ? datalike : 1;
  
  std::vector<TFile*> fin_data, fin_eff;
  std::vector<RooWorkspace*> wsp;
  std::vector<std::vector<RooDataSet*>> data;
  std::vector<RooAbsReal*> effC, effW;
  std::vector<TH3D*> effCHist, effWHist;
  std::vector<TH1D*> intCHist, intWHist;
  std::vector< std::vector<double> > intCVec(years.size(), std::vector<double>(0));
  std::vector< std::vector<double> > intWVec(years.size(), std::vector<double>(0));
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular (0);
  std::vector<RooAbsPdf*> PDF_sig_ang_fullAngular_penalty (0);

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
    all_years += year;
    for (unsigned int is = 0; is < nSamples; is++) {
      isample.clear(); isample.assign( Form("%i",is) );
      sample.defineType(("data"+year+"_subs"+isample).c_str());
    }
  }
  // Construct a simultaneous pdf using category sample as index
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf", "simultaneous pdf", sample);
  RooSimultaneous* simPdf_penalty = new RooSimultaneous("simPdf_penalty", "simultaneous pdf with penalty term", sample);

  // Define boundary check (returning 0 in physical region and 1 outside)
  BoundCheck* boundary = new BoundCheck("bound","Physical region",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);

  // Define penalty term (parameters set to zero and will be set sample-by-sample)
  Penalty* penTerm = new Penalty("penTerm","Penalty term",*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,0,0,0,0);

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


    // create roodataset (in case data-like option is selected, only import the correct % of data)
    RooDataSet* dataCT, *dataWT;
    std::vector<RooDataSet*> data_isample;

    if (datalike>0){
      for (unsigned int is = 0; is < nSamples; is++) {
        dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[years[iy]], (is+1)*scale_to_data[years[iy]] )) ;
        dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin))
          ->reduce( RooArgSet(reco_vars), Form("rand > %f && rand < %f", is*scale_to_data[years[iy]], (is+1)*scale_to_data[years[iy]] )) ;

        RooDataSet* datatmp = new RooDataSet(*dataCT,("data_"+shortString + Form("_subs%i", is)).c_str());
        datatmp->append(*dataWT);
        datatmp->Print();
        data_isample.push_back (datatmp);
      }
    }
    else{
      dataCT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_ctRECO_ev_b%i":"data_ctRECO_od_b%i"),q2Bin)) ;
      dataWT = (RooDataSet*)wsp[iy]->data(Form((parity==1?"data_wtRECO_ev_b%i":"data_wtRECO_od_b%i"),q2Bin)) ;
    
      RooDataSet* datatmp = new RooDataSet(*dataCT,("data_"+shortString + "_subs0").c_str());
      datatmp->append(*dataWT);
      data_isample.push_back (datatmp);
    }

    data.push_back(data_isample) ;
    for (unsigned int is = 0; is < nSamples; is++) {
      if ( !data[iy][is] || data[iy][is]->IsZombie() ) {
        cout<<"Dataset " << is  << "not found in file: "<<filename_data<<endl;
        return;
      }
    }

    // define angular PDF for signal, using the custom class
    // efficiency function and integral values are passed as arguments
    PDF_sig_ang_fullAngular.push_back( new PdfSigAng(("PDF_sig_ang_fullAngular_"+shortString+"_"+year).c_str(),
                                                     ("PDF_sig_ang_fullAngular_"+year).c_str(),
      		                                     *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
      		                                     *effC[iy], *effW[iy], intCVec[iy],intWVec[iy]));
    // define PDF with penalty term
    PDF_sig_ang_fullAngular_penalty.push_back( new PdfSigAng_Pen(("PDF_sig_ang_fullAngular_penalty_"+shortString+"_"+year).c_str(),
								 ("PDF_sig_ang_fullAngular_penalty_"+year).c_str(),
								 *ctK,*ctL,*phi,*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p,*mFrac,
								 *effC[iy], *effW[iy], intCVec[iy],intWVec[iy],*penTerm));

    for (unsigned int is = 0; is < nSamples; is++) {
      // insert sample in the category map, to be imported in the combined dataset
      map.insert( map.cbegin(), std::pair<const string,RooDataSet*>(("data"+year+Form("_subs%d",is)).c_str(), data[iy][is]) );
      // associate model with the data
      simPdf->addPdf(*PDF_sig_ang_fullAngular[iy], ("data"+year+Form("_subs%d",is)).c_str());
      simPdf_penalty->addPdf(*PDF_sig_ang_fullAngular_penalty[iy], ("data"+year+Form("_subs%d",is)).c_str());
    }
  }


  TFile* fout = new TFile(("simFitResults/simFitResult_recoMC_fullAngular" + all_years + stat + Form("_b%i.root", q2Bin)).c_str(),"RECREATE");

  // Construct combined dataset in (x,sample)
  RooDataSet allcombData ("allcombData", "combined data", 
                            vars,
                            Index(sample), 
                            Import(map)); 
  RooDataSet * combData; 
  RooAbsReal* nll;                         
  RooAbsReal* nll_penalty;                         

  // Results' containers
  RooRealVar* fitTime = new RooRealVar("fitTime","fit time",0,"s");
  RooRealVar* co1 = new RooRealVar("co1","Coefficient 1",0);
  RooRealVar* co4 = new RooRealVar("co4","Coefficient 4",0);
  RooRealVar* co5 = new RooRealVar("co5","Coefficient 5",0);
  RooArgList pars (*Fl,*P1,*P2,*P3,*P4p,*P5p,*P6p,*P8p);
  RooArgSet savePars (*co1,*co4,*co5,*fitTime);
  savePars.add(pars);
  RooDataSet* subResults = new RooDataSet("subResults","subResults",savePars);
  RooDataSet* subPosConv = new RooDataSet("subPosConv","subPosConv",savePars);
  RooDataSet* subPosNotc = new RooDataSet("subPosNotc","subPosNotc",savePars);
  RooDataSet* subNegConv = new RooDataSet("subNegConv","subNegConv",savePars);
  RooDataSet* subNegNotc = new RooDataSet("subNegNotc","subNegNotc",savePars);

  // Timer for fitting time
  TStopwatch subTime;

  // counters to monitor results' status
  int cnt[9];
  for (int iCnt=0; iCnt<9; ++iCnt) cnt[iCnt] = 0;
                         
  for (unsigned int is = 0; is < nSamples; is++) {
    string the_cut = Form("sample==sample::data%d_subs%d", years[0], is);
    if (years.size() > 1){
      for (unsigned int iy=1; iy < years.size(); iy++){
        the_cut = the_cut + Form("|| sample==sample::data%d_subs%d", years[iy], is);
      }
    }

    combData = (RooDataSet*)allcombData.reduce(Cut(the_cut.c_str()));
//     for(auto it = map.cbegin(); it != map.cend(); ++it)
//       std::cout << "dataset: " << it->first << ", with n entries: " << it->second->sumEntries() << "\n";

    // to start the fit, parameters are restored to the center of the parameter space
    Fl ->setVal(0.5);
    P1 ->setVal(0);
    P2 ->setVal(0);
    P3 ->setVal(0);
    P4p->setVal(0);
    P5p->setVal(0);
    P6p->setVal(0);
    P8p->setVal(0);

    // set penalty term power parameter
    penTerm->setPower(power/combData->numEntries());

    // define variables needed for adaptive procedure
    bool isPhysical = false;
    double coeff1 = 0;
    double coeff4 = 0;
    double coeff5 = 0;
    double sinPh1;
    int totCoeff, phase1, phase2, totCoeffPro;

    nll = simPdf->createNLL(*combData,
                            RooFit::Extended(kFALSE),
                            RooFit::NumCPU(1)
                            );
         
    RooMinimizer m(*nll) ;
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    // m.setVerbose(kTRUE);
    m.setPrintLevel(-1);
  //  Minuit2.setEps(1e-16) ;
    m.setMinimizerType("Minuit2");

    subTime.Start(true);

    m.setStrategy(0);
  //   m.setEvalErrorWall(false);
    m.migrad() ;
    m.hesse() ;
    // std::cout << std::endl;
    // std::cout << "######################### now strategy 2 #########################"<< std::endl;
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;
    // m.minos() ;
    
    RooFitResult* fitResult = m.save(("result_" + shortString + Form("subs%d",is)).c_str()) ; 
    // fitResult->Print("v");

    RooFitResult* fitResult_penalty = 0;
    bool usedPenalty = false;
    if ( fitResult->status()!=0 || fitResult->covQual()!=3 || boundary->getValV() > 0 ) {
      usedPenalty = true;

      for (totCoeff=0; fac1*pow(base1,totCoeff)<=maxCoeff; ++totCoeff) {

	for (phase1=0; phase1<=totCoeff; ++phase1) {
	  coeff1 = fac1 * pow(base1,totCoeff*cos(0.5*TMath::Pi()*phase1/(totCoeff>0?totCoeff:1)));
	  penTerm->setCoefficient(1,coeff1);

	  sinPh1 = sin(0.5*TMath::Pi()*phase1/(totCoeff>0?totCoeff:1));
	  totCoeffPro = (int)(totCoeff*sinPh1);
	  for (phase2=0; phase2<=totCoeffPro; ++phase2) {
	    coeff4 = fac4 * pow(base4,totCoeff*sinPh1*sin(0.5*TMath::Pi()*phase2/(totCoeffPro>0?totCoeffPro:1)));
	    coeff5 = fac5 * pow(base5,totCoeff*sinPh1*cos(0.5*TMath::Pi()*phase2/(totCoeffPro>0?totCoeffPro:1)));
	    // if (coeff4>2000) continue;
	    // if (coeff5>20000) continue;

	    penTerm->setCoefficient(4,coeff4);
	    penTerm->setCoefficient(5,coeff5);

	    nll_penalty = simPdf_penalty->createNLL(*combData,
						    RooFit::Extended(kFALSE),
						    RooFit::NumCPU(1)
						    );

	    RooMinimizer m_penalty (*nll_penalty) ;
	    m_penalty.optimizeConst(kTRUE);
	    m_penalty.setOffsetting(kTRUE);
	    // m_penalty.setVerbose(kTRUE);
	    m_penalty.setMinimizerType("Minuit2");
	    // m_penalty.setProfile(kTRUE);
	    m_penalty.setPrintLevel(-1);
	    m_penalty.setStrategy(2);
    
	    m_penalty.migrad() ;
	    m_penalty.hesse() ;
	    fitResult_penalty = m_penalty.save(Form("subRes_%s_%i",shortString.c_str(),is),Form("subRes_%s_%i",shortString.c_str(),is));
	    
	    // cout<<penTerm->getCoefficient(1)<<"\t"<<penTerm->getCoefficient(5)<<"\t"<<P5p->getValV()<<endl;
	    // fitResult_penalty->Print("v");
	    
	    if ( fitResult_penalty->status()==0 && fitResult_penalty->covQual()==3 && boundary->getValV()==0 ) {
	      isPhysical = true;
	      break;
	    }
	  }
	  if (isPhysical) break;
	}
	if (isPhysical) break;
      }
    
    }

    subTime.Stop();
    fitTime->setVal(subTime.CpuTime());
    // fitTime->setVal(subTime.RealTime());

    co1->setVal(0);
    co4->setVal(0);
    co5->setVal(0);

    if (usedPenalty) {
      if (isPhysical) {
	cout<<"Physical result with coeff: "<<coeff1<<" "<<coeff4<<" "<<coeff5<<endl;
	co1->setVal(coeff1);
	co4->setVal(coeff4);
	co5->setVal(coeff5);
      }
      else cout<<"No physical result"<<endl;
    }

    double boundCheck = boundary->getValV();

    ++cnt[8];
    int iCnt = 0;
    if (usedPenalty) {
      iCnt += 1;
      if (fitResult_penalty->status()!=0 || fitResult_penalty->covQual()!=3) iCnt += 4;
    } else
      if (fitResult->status()!=0 || fitResult->covQual()!=3) iCnt += 4;
    if (boundCheck>0) iCnt += 2;
    ++cnt[iCnt];

    subResults->add(savePars);
    if (boundCheck>0) {
      if (fitResult->status()==0 && fitResult->covQual()==3) subNegConv->add(savePars);
      else subNegNotc->add(savePars);
    } else {
      if (fitResult->status()==0 && fitResult->covQual()==3) subPosConv->add(savePars);
      else subPosNotc->add(savePars);
    }

    // Save fit results in file
    if (save) {
      fout->cd();
      if (usedPenalty) fitResult_penalty->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
      else fitResult->Write(("simFitResult_"+shortString+ Form("subs%d",is)).c_str(),TObject::kWriteDelete);
    }
  }  

  double time90quant = 0;
  double quant = 0;
  double totEntries = subResults->sumEntries();
  for (time90quant = 0; quant<0.9; time90quant += 0.1)
    quant = subResults->sumEntries(Form("fitTime<%.2f",time90quant))/totEntries;
  cout<<"Average fit time: "<<subResults->mean(*fitTime)<<" sec (90% quantile: "<<time90quant<<" sec)"<<endl;

  cout<<"Fitted subsamples: "<<cnt[8]<<" of which good: "<<cnt[0]+cnt[1]<<" ("<<cnt[1]<<" with the use of the penalty term)"<<endl;
  cout<<"Bad fits: "<<cnt[3]<<" converging outside physical region, "<<cnt[5]+cnt[7]<<" not converged ("<<cnt[5]<<" in ph region)"<<endl;

  if (save) {
    RooWorkspace* wksp = new RooWorkspace(("ws_"+shortString).c_str(),"Workspace with RECO subsamples fit results");
    wksp->import(*subResults);
    wksp->import(*subPosConv);
    wksp->import(*subPosNotc);
    wksp->import(*subNegConv);
    wksp->import(*subNegNotc);
    if (nSamples==1) wksp->import(*nll);

    fout->cd();
    wksp->Write();
  }


  if (!plot || datalike > 1) {
    fout->Close();
    return;
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

  TCanvas* cNll = new TCanvas("nllerrorhandling","nllerrorhandling",1200,900) ;
  cNll->Divide(3,3) ;
  cNll->cd(1) ; gPad->SetLeftMargin(0.15) ; frameFl->GetYaxis()->SetTitleOffset(1.4)  ; frameFl->Draw() ;
  cNll->cd(2) ; gPad->SetLeftMargin(0.15) ; frameP1->GetYaxis()->SetTitleOffset(1.4)  ; frameP1->Draw() ;
  cNll->cd(3) ; gPad->SetLeftMargin(0.15) ; frameP2->GetYaxis()->SetTitleOffset(1.4)  ; frameP2->Draw() ;
  cNll->cd(4) ; gPad->SetLeftMargin(0.15) ; frameP3->GetYaxis()->SetTitleOffset(1.4)  ; frameP3->Draw() ;
  cNll->cd(5) ; gPad->SetLeftMargin(0.15) ; frameP4p->GetYaxis()->SetTitleOffset(1.4) ; frameP4p->Draw() ;
  cNll->cd(6) ; gPad->SetLeftMargin(0.15) ; frameP5p->GetYaxis()->SetTitleOffset(1.4) ; frameP5p->Draw() ;
  cNll->cd(7) ; gPad->SetLeftMargin(0.15) ; frameP6p->GetYaxis()->SetTitleOffset(1.4) ; frameP6p->Draw() ;
  cNll->cd(8) ; gPad->SetLeftMargin(0.15) ; frameP8p->GetYaxis()->SetTitleOffset(1.4) ; frameP8p->Draw() ;
  cNll->SaveAs(("test_nll_Constr_" + all_years + stat + Form("_b%i.pdf", q2Bin)).c_str());


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

    combData->plotOn(xframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year+"_subs0").c_str()), Name(("plData"+year).c_str()));
    combData->plotOn(yframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year+"_subs0").c_str()));
    combData->plotOn(zframe,MarkerColor(kRed+1),LineColor(kRed+1),Binning(40), Cut(("sample==sample::data"+year+"_subs0").c_str()));

    simPdf->plotOn(xframe,Slice(sample, ("data"+year+"_subs0").c_str()), ProjWData(RooArgSet(sample), *combData), LineWidth(1),Name(("plPDF"+year).c_str()));
    simPdf->plotOn(yframe,Slice(sample, ("data"+year+"_subs0").c_str()), ProjWData(RooArgSet(sample), *combData), LineWidth(1));
    simPdf->plotOn(zframe,Slice(sample, ("data"+year+"_subs0").c_str()), ProjWData(RooArgSet(sample), *combData), LineWidth(1));

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
  fout->Close(); 

}


void simfit_recoMC_fullAngularBin1(int q2Bin, int parity, bool plot, bool save, int datalike, std::vector<int> years, std::map<int,float> scale_to_data, double power)
{
  if ( parity==-1 )
    for (parity=0; parity<2; ++parity)
      simfit_recoMC_fullAngularBin(q2Bin, parity, plot, save, datalike,  years, scale_to_data, power);
  else
    simfit_recoMC_fullAngularBin(q2Bin, parity, plot, save, datalike, years, scale_to_data, power);
}

int main(int argc, char** argv)
{
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively
  // parity format: [0] even efficiency
  //                [1] odd efficiency
  //                [-1] for each parity recursively

  int q2Bin   = -1;
  int parity  = -1; 

  if ( argc >= 2 ) q2Bin   = atoi(argv[1]);
  if ( argc >= 3 ) parity  = atoi(argv[2]);

  double power = 0.5;

  if ( argc >= 4 ) power = atof(argv[3]);

  bool plot = true;
  bool save = true;
  int  datalike = 0;

  if ( argc >= 5 && atoi(argv[4]) == 0 ) plot     = false;
  if ( argc >= 6 && atoi(argv[5]) == 0 ) save     = false;
  if ( argc >= 7 && atoi(argv[6]) > 0  ) datalike = atoi(argv[6]);

  std::vector<int> years;
  if ( argc >= 8 && atoi(argv[7]) != 0 ) years.push_back(atoi(argv[7]));
  else {
    cout << " no specific years selected, using default: 2016";
    years.push_back(2016);
  }
  if ( argc >= 9  && atoi(argv[8]) != 0 ) years.push_back(atoi(argv[8]));
  if ( argc >= 10 && atoi(argv[9]) != 0 ) years.push_back(atoi(argv[9]));

  if ( q2Bin   < -1 || q2Bin   >= nBins ) return 1;
  if ( parity  < -1 || parity  > 1      ) return 1;

  if ( q2Bin==-1 )   cout << "Running all the q2 bins" << endl;
  if ( parity==-1 )  cout << "Running both the parity datasets" << endl;
  if ( datalike )    cout << "Considering data-like statistics" << endl;

  std::map<int,float> scale_to_data;
  // https://docs.google.com/spreadsheets/d/1gG-qowySO9WJpMmr_bAWmOAu05J8zr95yJXGIYCY9-A/edit?usp=sharing
  scale_to_data.insert(std::make_pair(2016, 0.006*2 /2.5  )); // *2 since we are using only odd/even events, second factor is "data-driven"
  scale_to_data.insert(std::make_pair(2017, 0.005*2 /2.05 ));
  scale_to_data.insert(std::make_pair(2018, 0.007*2 /1.9  ));

  if ( q2Bin==-1 )
    for (q2Bin=0; q2Bin<nBins; ++q2Bin)
      simfit_recoMC_fullAngularBin1(q2Bin, parity, plot, save, datalike, years, scale_to_data, power);
  else
    simfit_recoMC_fullAngularBin1(q2Bin, parity, plot, save, datalike, years, scale_to_data, power);

  return 0;

}
