using namespace RooFit;
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

double PDGB0Mass = 5.27958;
double PDGJpsiMass = 3.096916;
double PDGPsiPrimeMass = 3.686109;
double PDGKstMass = 0.896;

TCanvas* c [nBins];

void createDataset(int year, int q2Bin, int data = 0, bool plot = false)
{
  // year format: [6] for 2016
  //              [7] for 2017
  //              [8] for 2018
  // q2-bin format: [0-8] for one bin
  //                [-1] for each bin recursively

  if ( q2Bin<-1 || q2Bin>=nBins ) return;

  if ( year<6 || year>8 ) return;

  bool isJpsi = false;
  bool isPsi  = false;
  bool isLMNR = false;

  cout << "q2Bin = " << q2Bin << endl;

  if (q2Bin==4)      isJpsi = true;
  else if (q2Bin==6) isPsi = true;
  else isLMNR = true;

  string dataString = data ? "DATA" : "MC";

  // define angular variables and variable for PU-reweighting
  RooRealVar ctK ("ctK","cos(#theta_{K})",-1,1);
  RooRealVar ctL ("ctL","cos(#theta_{L})",-1,1);
  RooRealVar phi ("phi","#phi",-TMath::Pi(),TMath::Pi());
  RooArgSet vars (ctK, ctL, phi);
  RooRealVar wei ("weight","weight",1);
  RooRealVar mass("mass","mass", 0,10);
  // random variable [0,1] to keep or reject the event when performing data-like stat. studies
  RooRealVar rand("rand", "rand", 0,1);
  TRandom rand_gen(1029);
  RooArgSet reco_vars (ctK, ctL, phi, mass, rand);
  if (data==0) reco_vars.add(wei);

  // flags to mark which q2 bins should be filled
  bool runBin [nBins];
  string shortString [nBins];
  string longString  [nBins];
  for (int i=0; i<nBins; ++i) {
    runBin [i] = false;
    if ( q2Bin!=-1 && q2Bin!=i ) continue;
    runBin [i] = true;
    shortString [i] = Form("b%i",i);
    longString  [i] = Form("q2 bin %i",i);
  }

  // Load ntuples
  TChain* t_num = new TChain();
  string year_str = Form("201%i", year);
  if (data==0 && isLMNR)
    t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iMC_LMNR.root/ntuple", year, year));
  else if (data==0 && isJpsi)
    t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iMC_JPSI.root/ntuple", year, year));
  else if (data==0 && isPsi)
    t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iMC_PSI.root/ntuple", year, year));
  else
    t_num->Add(Form("/eos/cms/store/user/fiorendi/p5prime/201%i/skims/newphi/201%iData_All_finalSelection.root/ntuple", year, year));
  int numEntries = t_num->GetEntries();
  std::cout << numEntries << std::endl;
  int counter;

  // Import branches from ntuples:
  // angular variables
  double recoCosThetaK, recoCosThetaL, recoPhi;
  t_num->SetBranchAddress( "cos_theta_k"     , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l"     , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu"    , &recoPhi       );

  // dimuon mass variables
  double recoDimuMass;
  t_num->SetBranchAddress( "mumuMass", &recoDimuMass );

  // B0 mass variable
  double recoB0Mass;
  t_num->SetBranchAddress( "tagged_mass", &recoB0Mass );

  bool passB0Psi_lmnr, passB0Psi_jpsi, passB0Psi_psip;
  t_num->SetBranchAddress( "passB0Psi_lmnr", &passB0Psi_lmnr );
  t_num->SetBranchAddress( "passB0Psi_jpsi", &passB0Psi_jpsi );
  t_num->SetBranchAddress( "passB0Psi_psip", &passB0Psi_psip );

  // B0-kinematic variables
  // double recoB0pT, recoB0eta;
  // t_num->SetBranchAddress( "bPt"    , &recoB0pT  );
  // t_num->SetBranchAddress( "bEta"   , &recoB0eta );

  double tagB0;
  t_num->SetBranchAddress( "tagB0"    , &tagB0     );

  // event number for even/odd splitting
  Long64_t eventN;
  t_num->SetBranchAddress( "eventN", &eventN     );

  Long64_t runN;
  t_num->SetBranchAddress("runN", &runN);

  // cut to remove B+->Psi(2S)K->Jpsi pi pi K
  // will be a boolean in ntuples in the future
  // keep it here for now since not finalized
  // as from https://github.com/CMSKStarMuMu/RooSideBand/blob/master/testSidebandFit.cc#L2592-L2661

  double wt_mass, wt_kstarmass, kaonPt, pionPt, mmpiMass, mmkMass;
  t_num->SetBranchAddress( "wt_mass",      &wt_mass      );
  t_num->SetBranchAddress( "wt_kstarmass", &wt_kstarmass );
  t_num->SetBranchAddress( "kaonPt",       &kaonPt       );
  t_num->SetBranchAddress( "pionPt",       &pionPt       );
  t_num->SetBranchAddress( "mmpiMass",     &mmpiMass     );
  t_num->SetBranchAddress( "mmkMass",      &mmkMass      );

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

  int xBin;

  // event pileup weight
  if (data==0) {
    // flavour tagging variables
    double genSignal;
    t_num->SetBranchAddress( "genSignal", &genSignal );
    
    std::cout << "is MC"  << std::endl;
    float PUweight = 1;
    t_num->SetBranchAddress( "weight", &PUweight );

    // Define datasets for five efficiency terms
    RooDataSet* data_ctRECO_ev [nBins];
    RooDataSet* data_ctRECO_od [nBins];
    RooDataSet* data_wtRECO_ev [nBins];
    RooDataSet* data_wtRECO_od [nBins];

    for (int i=0; i<nBins; ++i) if (runBin[i]){
        data_ctRECO_ev [i] = new RooDataSet( ("data_ctRECO_ev_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (even)",
    					     reco_vars, "weight" );
        data_ctRECO_od [i] = new RooDataSet( ("data_ctRECO_od_"+shortString[i]).c_str(), "Correctly-tagged reconstructed candidates after selections (odd)",
					     reco_vars, "weight" );
        data_wtRECO_ev [i] = new RooDataSet( ("data_wtRECO_ev_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (even)",
					     reco_vars, "weight" );
        data_wtRECO_od [i] = new RooDataSet( ("data_wtRECO_od_"+shortString[i]).c_str(), "Wrongly-tagged reconstructed candidates after selections (odd)",
					     reco_vars, "weight" );
    }

    // Prepare numerator dataset
    cout<<"Starting numerator dataset filling..."<<endl;
    counter=0;
    for (int iCand=0; iCand<numEntries; ++iCand) {
      t_num->GetEntry(iCand);
     
      if( (year == 6) && (runN >= 272007) && (runN <= 277990) ){continue;} // remove problematic year of 2016

      // anti-radiation cut
      if (isLMNR && passB0Psi_lmnr == 0) continue;
      else if (isJpsi && passB0Psi_jpsi == 0) continue;
      else if (isPsi  && passB0Psi_psip == 0)  continue;

      // find q2 bin of current candidate
      xBin=-1;
      for (int i=0; i<nBins; ++i)
        if ( runBin[i] )
    	if ( ( pow(recoDimuMass,2) < binBorders[i+1] ) &&
    	     ( pow(recoDimuMass,2) > binBorders[i]   ) ) {
   	  xBin = i;
  	  break;
  	}
      if (xBin<0) continue;

      // apply cut for bin 4 
      bool XCut= (( (PDGB0Mass - wt_mass) - y0Cut ) / (y1Cut-y0Cut)) < (((wt_kstarmass-PDGKstMass)-x0Cut) / (x1Cut-x0Cut)) && \
                    kaonPt > pionPt && \
                    (wt_kstarmass-PDGKstMass)>0 && \
                    (mmpiMass > CutX1 && mmpiMass < CutX2) && \
                    (mmkMass >  CutY1 && mmkMass  < CutY2) && \
                    ((mmkMass - y_0Cut) / (y_1Cut - y_0Cut)) > ((mmpiMass-x_0Cut)/(x_1Cut-x_0Cut));

      if (XCut && xBin == 4) continue;

      // status display
      if ( iCand > 1.0*counter*numEntries/100 ) {
        cout<<counter<<"%"<<endl;
        counter += 10;
      }
      //generate random variable [0,1]
      // fill dataset
      ctK.setVal(recoCosThetaK);
      ctL.setVal(recoCosThetaL);
      phi.setVal(recoPhi);
      mass.setVal(recoB0Mass);
      rand.setVal(rand_gen.Uniform(1));
  
      if (genSignal != tagB0+1) { // correctly tagged events
        if (eventN%2==0) data_ctRECO_ev[xBin]->add(reco_vars,PUweight);
        else data_ctRECO_od[xBin]->add(reco_vars,PUweight);
      } else { // wrongly tagged events
      if (eventN%2==0) data_wtRECO_ev[xBin]->add(reco_vars,PUweight);
      else data_wtRECO_od[xBin]->add(reco_vars,PUweight);
      }
    }
    cout<<"Dataset prepared"<<endl;

    // Save datasets in workspaces
    RooWorkspace *ws_ev [nBins];
    RooWorkspace *ws_od [nBins];

    for (int i=0; i<nBins; ++i) if (runBin[i]) {
        // Skip the creation of a file when the correct-tag efficiency cannot be computed (empty numerators)
        // which usually means that either you are using a resonant MC, which does not fill signal q2 bins,
        // or using a bin too fine, or out of range
        // If correct-tag numerator is filled and wrong-tag is not, a warning is returned
        if ( data_ctRECO_ev[i]->numEntries()==0 || data_ctRECO_od[i]->numEntries()==0 ) {
          cout<<"Error: ctRECO is empty in q2 bin "<<i<<endl;
	  continue;
        }
        if ( data_wtRECO_ev[i]->numEntries()==0 || data_wtRECO_od[i]->numEntries()==0 ) {
	  cout<<"Warning: wtRECO is empty in q2 bin "<<i<<endl;
        }
        ws_ev[i] = new RooWorkspace(("ws_"+shortString[i]+"p0").c_str(),"Workspace with single-bin even datasets");
        ws_od[i] = new RooWorkspace(("ws_"+shortString[i]+"p1").c_str(),"Workspace with single-bin odd datasets");
        ws_ev[i]->import( *data_ctRECO_ev[i] );
        ws_od[i]->import( *data_ctRECO_od[i] );
        ws_ev[i]->import( *data_wtRECO_ev[i] );
        ws_od[i]->import( *data_wtRECO_od[i] );
        TFile* fout = new TFile( ( "reco"+dataString+"Dataset_"+shortString[i]+ "_" + year_str + ".root" ).c_str(), "RECREATE"); 
        ws_ev[i]->Write();
        ws_od[i]->Write();
        fout->Close();

     }

    // Plot 1D distributions of datasets
    if (plot) {

      // to keep all distributions visible in the same plot, the ones with higher stats (tipically denominators) need to be rescaled
      double rescFac1 = 1.0/12;
      double rescFac2 = 1.0;
      double rescFac3 = 1.0/25;

      TLegend* leg = new TLegend(0.4,0.8,0.9,0.9);

      RooPlot* xframe_ev [nBins];
      RooPlot* yframe_ev [nBins];
      RooPlot* zframe_ev [nBins];
      RooPlot* xframe_od [nBins];
      RooPlot* yframe_od [nBins];
      RooPlot* zframe_od [nBins];

      bool legFilled = false;

      RooDataSet* data_RECO_ev [nBins];
      RooDataSet* data_RECO_od [nBins];

      for (int i=0; i<nBins; ++i) if (runBin[i]) {
   	// Create dataset containing both correct-tag and wrong-tag events
	data_RECO_ev [i] = new RooDataSet( ("data_RECO_ev_"+shortString[i]).c_str(), "Reconstructed candidates after selections (even)", data_ctRECO_ev[i], vars );
	data_RECO_od [i] = new RooDataSet( ("data_RECO_od_"+shortString[i]).c_str(), "Reconstructed candidates after selections (odd)" , data_ctRECO_od[i], vars );
	data_RECO_ev[i]->append(*(data_wtRECO_ev[i]));
	data_RECO_od[i]->append(*(data_wtRECO_od[i]));

	// create frames (one for each bin/parity/variable, but all the six efficiency terms are plotted together)
	c [i] = new TCanvas(("c_"+shortString[i]).c_str(),("Num and Den 1D projections - "+longString[i]).c_str(),2000,1400);
	xframe_ev [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (even)").c_str()));
	yframe_ev [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (even)").c_str()));
	zframe_ev [i] = phi.frame(Title((longString[i]+" #phi distributions (even)").c_str()));
	xframe_od [i] = ctK.frame(Title((longString[i]+" cos(#theta_{K}) distributions (odd)").c_str()));
	yframe_od [i] = ctL.frame(Title((longString[i]+" cos(#theta_{L}) distributions (odd)").c_str()));
	zframe_od [i] = phi.frame(Title((longString[i]+" #phi distributions (odd)").c_str()));

	// plot datasets on frames
	if (!legFilled) { // the first time assign names to tag them in the legend
	  data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40),Name("plCTrecoNum"));
	  data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40),Name("plWTrecoNum"));
	  data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40),Name("plRecoNum"));
	} else {
	  data_ctRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	  data_wtRECO_ev[i]->plotOn(xframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	  data_RECO_ev  [i]->plotOn(xframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	}
	data_ctRECO_od[i]->plotOn(xframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_od[i]->plotOn(xframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_od  [i]->plotOn(xframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_ctRECO_ev[i]->plotOn(yframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_ev[i]->plotOn(yframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_ev  [i]->plotOn(yframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_ctRECO_od[i]->plotOn(yframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_od[i]->plotOn(yframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_od  [i]->plotOn(yframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_ctRECO_ev[i]->plotOn(zframe_ev[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_ev[i]->plotOn(zframe_ev[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_ev  [i]->plotOn(zframe_ev[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));
	data_ctRECO_od[i]->plotOn(zframe_od[i],MarkerColor(kMagenta) ,LineColor(kMagenta) ,Binning(40));
	data_wtRECO_od[i]->plotOn(zframe_od[i],MarkerColor(kViolet-3),LineColor(kViolet-3),Binning(40));
	data_RECO_od  [i]->plotOn(zframe_od[i],MarkerColor(kBlack)   ,LineColor(kBlack)   ,Binning(40));

	xframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
	yframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
	zframe_ev[i]->GetYaxis()->SetTitleOffset(1.6);
	xframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
	yframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
	zframe_od[i]->GetYaxis()->SetTitleOffset(1.6);
	xframe_ev[i]->SetMaximum(xframe_ev[i]->GetMaximum()*rescFac1);
	yframe_ev[i]->SetMaximum(yframe_ev[i]->GetMaximum()*rescFac1);
	zframe_ev[i]->SetMaximum(zframe_ev[i]->GetMaximum()*rescFac1);
	xframe_od[i]->SetMaximum(xframe_od[i]->GetMaximum()*rescFac1);
	yframe_od[i]->SetMaximum(yframe_od[i]->GetMaximum()*rescFac1);
	zframe_od[i]->SetMaximum(zframe_od[i]->GetMaximum()*rescFac1);

	// Fill legend (only one time)
	if (!legFilled) {
	  string strRescFac1 = (rescFac1<1?Form(" (*%1.2f)",rescFac1):"");
	  string strRescFac2 = (rescFac2<1?Form(" (*%1.2f)",rescFac2):"");
	  string strRescFac3 = (rescFac3<1?Form(" (*%1.2f)",rescFac3):"");
	  leg->AddEntry(xframe_ev[i]->findObject("plRecoNum"  ),"Post-selection distribution","lep");
	  leg->AddEntry(xframe_ev[i]->findObject("plCTrecoNum"),"Post-selection correct-tag distribution","lep");
	  leg->AddEntry(xframe_ev[i]->findObject("plWTrecoNum"),"Post-selection wrong-tag distribution","lep");
	  legFilled = true;
	}

	// Plot even distributions in the top row and odd ones in the bottom row
	c[i]->Divide(3,2);
	c[i]->cd(1);
	gPad->SetLeftMargin(0.17); 
	xframe_ev[i]->Draw();
	leg->Draw("same");
	c[i]->cd(2);
	gPad->SetLeftMargin(0.17); 
	yframe_ev[i]->Draw();
	leg->Draw("same");
	c[i]->cd(3);
	gPad->SetLeftMargin(0.17); 
	zframe_ev[i]->Draw();
	leg->Draw("same");
	c[i]->cd(4);
	gPad->SetLeftMargin(0.17); 
	xframe_od[i]->Draw();
	leg->Draw("same");
	c[i]->cd(5);
	gPad->SetLeftMargin(0.17); 
	yframe_od[i]->Draw();
	leg->Draw("same");
	c[i]->cd(6);
	gPad->SetLeftMargin(0.17); 
	zframe_od[i]->Draw();
	leg->Draw("same");

	c[i]->SaveAs( ("plotDist_d/dist_RECO_"+shortString[i]+".pdf").c_str() );
      }
    }
  }

  else{
    RooDataSet* data [nBins];
    for (int i=0; i<nBins; ++i) {
      if (runBin[i] ){
        data [i] = new RooDataSet( ("data_"+shortString[i]).c_str(), "Reconstructed candidates after selections",reco_vars);
      }
    }

    // Prepare numerator dataset
    cout<<"Starting numerator dataset filling..."<<endl;
    counter=0;
    for (int iCand=0; iCand<numEntries; ++iCand) {
      t_num->GetEntry(iCand);

      if( (year == 6) && (runN >= 272007) && (runN <= 277990) ){continue;} // remove problematic year of 2016

      // anti-radiation cut
      if (isLMNR && passB0Psi_lmnr == 0) continue;
      else if (isJpsi && passB0Psi_jpsi == 0) continue;
      else if (isPsi  && passB0Psi_psip == 0)  continue;

      // find q2 bin of current candidate
      xBin=-1;
      for (int i=0; i<nBins; ++i)
        if ( runBin[i] )
  	if ( ( pow(recoDimuMass,2) < binBorders[i+1] ) &&
  	     ( pow(recoDimuMass,2) > binBorders[i]   ) ) {
  	  xBin = i;
  	  break;
  	}
      if (xBin<0) continue;

      // apply cut for bin 4 
      bool XCut= (( (PDGB0Mass - wt_mass) - y0Cut ) / (y1Cut-y0Cut)) < (((wt_kstarmass-PDGKstMass)-x0Cut) / (x1Cut-x0Cut)) && \
                    kaonPt > pionPt && \
                    (wt_kstarmass-PDGKstMass)>0 && \
                    (mmpiMass > CutX1 && mmpiMass < CutX2) && \
                    (mmkMass >  CutY1 && mmkMass  < CutY2) && \
                    ((mmkMass - y_0Cut) / (y_1Cut - y_0Cut)) > ((mmpiMass-x_0Cut)/(x_1Cut-x_0Cut));

      if (XCut && xBin == 4) continue;

      // status display
      if ( iCand > 1.0*counter*numEntries/100 ) {
        cout<<counter<<"%"<<endl;
        counter += 10;
      }
      // fill dataset
      ctK.setVal(recoCosThetaK);
      ctL.setVal(recoCosThetaL);
      phi.setVal(recoPhi);
      mass.setVal(recoB0Mass);

      data[xBin]->add(reco_vars);
    }
    cout<<"Dataset prepared"<<endl;

    // Save datasets in workspaces
    RooWorkspace *ws [nBins];
    for (int i=0; i<nBins; ++i) if (runBin[i]) {
      if ( data[i]->numEntries()==0  ) {
        cout<<"Error: RECO data is empty in q2 bin "<<i<<endl;
  	continue;

      }
      ws[i] = new RooWorkspace(("ws_"+shortString[i]+"p0").c_str(),"Workspace with single-bin data datasets");
      ws[i]->import( *data[i] );
      TFile* fout = new TFile( ( "reco"+dataString+"Dataset_"+shortString[i]+ "_" + year_str + ".root" ).c_str(), "RECREATE"); 
      ws[i]->Write();
      fout->Close();
    }

  }


}

