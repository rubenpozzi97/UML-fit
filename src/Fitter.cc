#include "Riostream.h" 

#include "Fitter.h" 

ClassImp(Fitter) 

Fitter::Fitter(const char *_name, const char *_title,
	       RooArgList _angPars,
	       RooDataSet* _combData,
	       RooAbsPdf* _simPdf,
	       RooAbsPdf* _simPdf_penalty,
	       BoundCheck* _boundary,
	       BoundDist* _bound_dist,
	       Penalty* _penTerm,
               RooArgSet* _constrVars
	       ) :
  name(_name),
  title(_title),
  angPars(_angPars),
  combData(_combData),
  simPdf(_simPdf),
  simPdf_penalty(_simPdf_penalty),
  boundary(_boundary),
  bound_dist(_bound_dist),
  penTerm(_penTerm),
  constrVars(_constrVars)
{

  SetDefConf();
  usedPenalty = false;

}

Fitter::Fitter(const char *_name, const char *_title,
	       RooArgList _angPars,
	       RooDataSet* _combData,
	       RooAbsPdf* _simPdf,
	       RooAbsPdf* _simPdf_penalty,
	       BoundCheck* _boundary,
	       BoundDist* _bound_dist,
      	       RooArgSet* _constrVars 	       
	       ) :
  name(_name),
  title(_title),
  angPars(_angPars),
  combData(_combData),
  simPdf(_simPdf),
  simPdf_penalty(_simPdf_penalty),
  boundary(_boundary),
  bound_dist(_bound_dist),
  constrVars(_constrVars)
{

  penTerm = 0;

  SetDefConf();
  usedPenalty = false;

}

Fitter::Fitter(const Fitter& other, const char* name) :  
  name(other.name),
  title(other.title),
  angPars(other.angPars),
  combData(other.combData),
  simPdf(other.simPdf),
  simPdf_penalty(other.simPdf_penalty),
  boundary(other.boundary),
  bound_dist(other.bound_dist),
  penTerm(other.penTerm),
  constrVars(other.constrVars)
{


  fac1 = other.fac1;
  fac4 = other.fac4;
  base1 = other.base1;
  base4 = other.base4;
  max1 = other.max1;
  max4 = other.max4;
  base1_corr = other.base1_corr;
  base4_corr = other.base4_corr;
  min_base = other.min_base;

  maxCoeff = other.maxCoeff;
  coeff1 = other.coeff1;
  coeff4 = other.coeff4;
  coeff5 = other.coeff5;

  usedPenalty = other.usedPenalty;

}

void Fitter::SetDefConf()
{
  fac1 = 0.1;
  fac4 = 0.01;
  base1 = 0.05;
  base4 = 0.05;
  max1 = 0.;
  max4 = 0.;

  min_base = 1.05;

  base1_corr = base1*sqrt(combData->numEntries());
  base4_corr = base4*sqrt(combData->numEntries());
  if (base1_corr<min_base) base1_corr = min_base;
  if (base4_corr<min_base) base4_corr = min_base;

  maxCoeff = 1e8;
  coeff1 = 0;
  coeff4 = 0;
  coeff5 = 0;

  boundDist = -1.;
  boundDistTime = 0.;

  minParError = 0.01;
  widthScale = 0.1;

  vFitResult = std::vector<Double_t>(angPars.getSize(),0);
  vFitErrLow  = std::vector<Double_t>(angPars.getSize(),0);
  vFitErrHigh = std::vector<Double_t>(angPars.getSize(),0);
  vImprovResult = std::vector<Double_t>(angPars.getSize(),0);

  vResult = std::vector<Double_t>(angPars.getSize(),0);
  vConfInterLow  = std::vector<Double_t>(angPars.getSize(),0);
  vConfInterHigh = std::vector<Double_t>(angPars.getSize(),0);

  plotString = name;
  plotDir = ".";

}


Int_t Fitter::fit()
{

    coeff1 = 0;
    coeff4 = 0;
    coeff5 = 0;

    // set up free fit
    nll = simPdf->createNLL(*combData,
                            RooFit::Extended(kFALSE),
                            RooFit::NumCPU(1)
                            );
    if ( constrVars != nullptr && constrVars->getSize() > 0  ){
      nll = simPdf->createNLL(*combData,
                              RooFit::Extended(kFALSE),
                              RooFit::NumCPU(1),
                              RooFit::Constrain(*constrVars)
                              );
    }                              
         
    RooMinimizer m(*nll) ;
    m.optimizeConst (kTRUE); // do not recalculate constant terms
    m.setOffsetting(kTRUE);  //  Enable internal likelihood offsetting for enhanced numeric precision.
    // m.setVerbose(kTRUE);
    m.setPrintLevel(-1);
    m.setPrintEvalErrors(-1);
    //  Minuit2.setEps(1e-16) ;
    m.setMinimizerType("Minuit2");

    // free fit
    m.setStrategy(0);
    // m.setEvalErrorWall(false);
    m.migrad() ;
    m.hesse() ;
    // std::cout << std::endl;
    // std::cout << "######################### now strategy 2 #########################"<< std::endl;
    m.setStrategy(2);
    m.migrad() ;
    m.hesse() ;
    m.minos() ;
    
    result_free = m.save("result") ;

    result_penalty = 0;
    usedPenalty = false;

    // if free fit is good return its result
    if ( result_free->status()==0 && result_free->covQual()==3 && boundary->getValV() == 0 ) {
      computeBoundaryDistance();
      fillResultContainers();
      return 0;
    }

    usedPenalty = true;

    // optional: if a partial boundary is satisfied
    // do not apply the corresponding penalty term
    // bool inCTL4  = true;
    // bool inCTL15 = true;
    // if ( !boundary->isInCTL4()  ) inCTL4  = false;
    // if ( !boundary->isInCTL15() ) inCTL15 = false;

    for (int totCoeff=0; fac1*pow(base1_corr,totCoeff)<=maxCoeff; ++totCoeff) {

      for (int iCoeff1=totCoeff; iCoeff1>=0; --iCoeff1) {

	// set penalty coefficients
	coeff1 = fac1 * pow(base1_corr,iCoeff1);
	if (max1>0 && coeff1>max1) continue;

	coeff4 = fac4 * pow(base4_corr,totCoeff-iCoeff1);
	if (max4>0 && coeff4>max4) continue;

	coeff5 = pow(coeff1,1.5) / 316.2;

	// optional: if a partial boundary is satisfied
	// do not apply the corresponding penalty term
	// if ( inCTL15 ) {
	//   if ( iCoeff1>0 ) continue;
	//   coeff1 = 0;
	//   coeff5 = 0;
	// }
	// if ( inCTL4 ) {
	//   if ( totCoeff-iCoeff1>0 ) continue;
	//   coeff4 = 0;
	// }

	penTerm->setCoefficient(1,coeff1);
	penTerm->setCoefficient(4,coeff4);
	penTerm->setCoefficient(5,coeff5);

	// set up the penalised fit
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
	m_penalty.setPrintEvalErrors(-1);
	m_penalty.setStrategy(2);
    
	// penalised fit
	m_penalty.migrad() ;
	m_penalty.hesse() ;
	result_penalty = m_penalty.save("result");
	    
	// cout<<penTerm->getCoefficient(1)<<"\t"<<penTerm->getCoefficient(5)<<"\t"<<P5p->getValV()<<endl;
	// result_penalty->Print("v");

	// if a good fit is found return good status
	if ( result_penalty->status()==0 && result_penalty->covQual()==3 ) {
	  if ( boundary->getValV()==0 ) {
	    // cout<<"P "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
	    computeBoundaryDistance();
	    fillResultContainers();
	    return 0;
	  } // else cout<<"O "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
	} // else cout<<"N "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;

      }

    }
    
    // if no good fit is found return bad status
    return 1;

}

Int_t Fitter::improveAng(int seed, int nGen)
{

  Double_t preBoundDist = boundDist;
  // std::cout<<"Starting with dist "<<preBoundDist<<std::endl;
  
  // Improve global fit result
  TRandom3 randGen (seed);

  double NLL_before = nll->getValV();
  double improvNLL = NLL_before;
  double testNLL = 0;
  // std::cout<<NLL_before;

  std::vector<double> vTestPar(angPars.getSize());
  std::vector<double> vImprovPar(angPars.getSize());
  for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
    vImprovPar[iPar] = ((RooRealVar*)angPars.at(iPar))->getValV();
    // std::cout<<"  "<<vImprovPar[iPar];
  } 
  // std::cout<<std::endl;

  int iImprove = 0;

  do {

    for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
      RooRealVar* par = (RooRealVar*)angPars.at(iPar);
      do vTestPar[iPar] = randGen.Gaus(vImprovPar[iPar],TMath::Max(boundDist,0.002));
      while (vTestPar[iPar]>par->getMax() || vTestPar[iPar]<par->getMin());
      par->setVal(vTestPar[iPar]);
    }

    if (boundary->getValV()>0) continue;

    testNLL = nll->getValV();
    if (improvNLL>testNLL) {
      improvNLL = testNLL;
      // std::cout<<NLL_before-improvNLL;
      for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
	vImprovPar[iPar] = vTestPar[iPar];
	// std::cout<<"  "<<vImprovPar[iPar];
      } 
      // std::cout<<std::endl;
    }

    ++iImprove;

  } while (iImprove<nGen);

  for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
    ((RooRealVar*)angPars.at(iPar))->setVal(vImprovPar[iPar]);
  
  computeBoundaryDistance();
  std::cout<<"Improved fit result: deltaNLL = "<<NLL_before-improvNLL<<" bound dist: "<<preBoundDist<<" -> "<<boundDist<<std::endl;

  fillResultContainers(true);

  return 0;

}


Int_t Fitter::improve(int seed, int nGen)
{

  if (!result_penalty) {
    std::cout<<"Error: the penalised fit need to be run before improving the result"<<std::endl;
    return 1;
  }

  double weight = 0.05;

  Double_t preBoundDist = boundDist;
  // std::cout<<"Starting full with dist "<<preBoundDist<<std::endl;
  
  // Improve global fit result
  TRandom3 randGen (seed);

  double NLL_before = nll->getValV();
  double improvNLL = NLL_before;
  double testNLL = 0;
  // std::cout<<NLL_before;

  std::vector<double> vTestPar(angPars.getSize());
  std::vector<double> vImprovPar(angPars.getSize());
  for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
    vImprovPar[iPar] = ((RooRealVar*)angPars.at(iPar))->getValV();
    // std::cout<<"  "<<vImprovPar[iPar];
  } 
  // std::cout<<std::endl;

  auto nllVars = nll->getVariables();
  auto fitPars = result_penalty->floatParsFinal();
  // printf("S ");
  // for (int iPar = 0; iPar < fitPars.getSize(); ++iPar)
  //   printf(" %f",((RooRealVar*)fitPars.at(iPar))->getValV());
  // printf("\n");

  int iImprove = 0;

  do {

    for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
      RooRealVar* par = (RooRealVar*)angPars.at(iPar);
      do vTestPar[iPar] = randGen.Gaus(vImprovPar[iPar],TMath::Max(boundDist,0.002));
      while (vTestPar[iPar]>par->getMax() || vTestPar[iPar]<par->getMin());
      par->setVal(vTestPar[iPar]);
    }

    if (boundary->getValV()>0) continue;

    auto floatPars = result_penalty->randomizePars();
    for (int iPar = 0; iPar < floatPars.getSize(); ++iPar) {
      RooRealVar* par = (RooRealVar*)floatPars.at(iPar);
      if (angPars.contains(*par)) continue; // angular parameters already set
      auto nllPar = (RooRealVar*)nllVars->find(*par);
      double fitPar = ((RooRealVar*)fitPars.find(*par))->getValV();
      if (nllPar) nllPar->setVal(weight*par->getValV()+(1.-weight)*fitPar);
    }

    testNLL = nll->getValV();
    if (improvNLL>testNLL) {
      improvNLL = testNLL;
      // std::cout<<NLL_before-improvNLL;
      for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {
	vImprovPar[iPar] = vTestPar[iPar];
	// std::cout<<"  "<<vImprovPar[iPar];
      } 
      // std::cout<<std::endl;
      // printf(" I ");
      for (int iPar = 0; iPar < floatPars.getSize(); ++iPar) {
	RooRealVar* par = (RooRealVar*)floatPars.at(iPar);
	improvVars[par->GetName()] = par->getValV();
	// if (angPars.contains(*par)) continue;
	// printf(" %f",improvVars[par->GetName()]);
      }
      // printf("\n");

    }

    ++iImprove;

  } while (iImprove<nGen);

  for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
    ((RooRealVar*)angPars.at(iPar))->setVal(vImprovPar[iPar]);
  
  computeBoundaryDistance();
  std::cout<<"Fully improved fit result: deltaNLL = "<<NLL_before-improvNLL<<" bound dist: "<<preBoundDist<<" -> "<<boundDist<<std::endl;

  fillResultContainers(true);

  return 0;

}


Int_t Fitter::Minos(int seed, int nGenMINOS, int nImprovMINOS, bool plotMINOS)
{

  // NLL of the best-fit result
  double NLL_min = nll->getValV();

  TRandom3 randGenMinos (seed);

  std::vector<double> vLastHit(angPars.getSize());

  // Loop over the parameters
  for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)angPars.at(iPar);

    // get and print the best-fit result
    double p_best = vResult[iPar];
    std::cout<<par->GetName()<<" best: "<<p_best<<std::endl;

      // vectors for TGraph plots
    std::vector<double> *vPval = 0;
    std::vector<double> *vdNLL = 0;
    if (plotMINOS) {
      vPval = new std::vector<double> (0);
      vdNLL = new std::vector<double> (0);
    }

    singleMinosInstance(NLL_min, iPar, false, nGenMINOS, &vLastHit, &randGenMinos, vPval, vdNLL);
    if (nImprovMINOS>0) improveMINOS(NLL_min, iPar, false, nImprovMINOS, &vLastHit, seed);

    singleMinosInstance(NLL_min, iPar, true, nGenMINOS, &vLastHit, &randGenMinos, vPval, vdNLL);
    if (nImprovMINOS>0) improveMINOS(NLL_min, iPar, true, nImprovMINOS, &vLastHit, seed);

    if (plotMINOS) doPlotMINOS(iPar, vPval, vdNLL);

  }

  return 0;

}


Int_t Fitter::singleMinosInstance(double NLL_min, int iPar, bool isErrHigh, int nGenMINOS, std::vector<double> *vLastHit, TRandom3* randGenMinos, std::vector<double> *vPval, std::vector<double> *vdNLL)
{

  double probedNLL;

  RooRealVar* par = (RooRealVar*)angPars.at(iPar);

  double p_best = vResult[iPar];

  // if vectors for plot exist, fill them
  if (!isErrHigh && vPval && vdNLL) {
    vPval->push_back(p_best);
    vdNLL->push_back(0);
  }

  for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1)
    vLastHit->at(iPar1) = vResult[iPar1];

  TH1D* parRandomPool = 0;
  int nHistBins = 0;
  if (isErrHigh) {
    nHistBins = (int)((par->getMax()-p_best)/0.0005);
    // nHistBins = (int)((par->getMax()-p_best+0.0005)/0.001);
    if (nHistBins<1) nHistBins=1;
    parRandomPool = new TH1D(Form("hRandPoolH%i_%s",iPar,name.Data()),
			     Form("hRandPoolH%i_%s",iPar,name.Data()),
			     nHistBins,par->getMax()-0.0005*nHistBins,par->getMax());
    // nHistBins,par->getMax()+0.0005-0.001*nHistBins,par->getMax()+0.0005);
  } else {
    nHistBins = (int)((p_best-par->getMin())/0.0005);
    // nHistBins = (int)((p_best-par->getMin()+0.0005)/0.001);
    if (nHistBins<1) nHistBins=1;
    parRandomPool = new TH1D(Form("hRandPoolL%i_%s",iPar,name.Data()),
			     Form("hRandPoolL%i_%s",iPar,name.Data()),
			     nHistBins,par->getMin(),par->getMin()+0.0005*nHistBins);
    // nHistBins,par->getMin()-0.0005,par->getMin()-0.0005+0.001*nHistBins);
  }
  double sigma = fabs(isErrHigh?vFitErrHigh[iPar]:vFitErrLow[iPar]);
  if (sigma<minParError) sigma = minParError;
  for (int iBin=1; iBin<=nHistBins; ++iBin) {
    double x = (parRandomPool->GetBinCenter(iBin)-p_best) / sigma;
    parRandomPool->SetBinContent(iBin,fabs(x)*exp(-0.5*x*x));
  }

  double p_in = p_best;
  double p_test = 0;

  int iPnt=0;

  do {

    do p_test = parRandomPool->GetRandom();
    while (p_test>par->getMax() || p_test<par->getMin());
    par->setVal(p_test);
    for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1) {
      if (iPar1==iPar) continue;
      RooRealVar* par1 = (RooRealVar*)angPars.at(iPar1);
      double par1val = 0;
      do par1val = randGenMinos->Gaus(vLastHit->at(iPar1),widthScale*TMath::Max(0.5*(vFitErrHigh[iPar1]-vFitErrLow[iPar1]),minParError));
      while (par1val>par1->getMax() || par1val<par1->getMin());
      par1->setVal(par1val);
    }
    // check if the point is physical
    if (boundary->getValV()>0) continue;
    // get and test the local likelihood
    probedNLL = nll->getValV();
    if (probedNLL<=NLL_min+0.5) {
      p_in = p_test;
      if (isErrHigh) { if ( p_in > par->getMax()-parRandomPool->GetBinWidth(1) ) break; }
      else if ( p_in < par->getMin()+parRandomPool->GetBinWidth(1) ) break;
      for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1) {
	RooRealVar* par1 = (RooRealVar*)angPars.at(iPar1);
	vLastHit->at(iPar1) = par1->getValV();
      }
      if (isErrHigh)
	for (int iBin=1; iBin<=parRandomPool->FindBin(p_test); ++iBin)
	  parRandomPool->SetBinContent(iBin,0);
      else 
	for (int iBin=parRandomPool->FindBin(p_test); iBin<=nHistBins; ++iBin)
	  parRandomPool->SetBinContent(iBin,0);
    } else {
      parRandomPool->Fill(p_test,0.02/(probedNLL-NLL_min-0.5));
    }

    if (vPval && vdNLL && probedNLL-NLL_min<4.5) {
      // fill the vectors for plot
      vPval->push_back(p_test);
      vdNLL->push_back(probedNLL-NLL_min);
    }

    ++iPnt;
    // apply conditions
  } while ( iPnt < nGenMINOS );

  if (isErrHigh) {
    vConfInterHigh[iPar] = p_in;
    std::cout<<par->GetName()<<" high: "<<p_in<<std::endl;
  } else {
    vConfInterLow[iPar] = p_in;
    std::cout<<par->GetName()<<" low:  "<<p_in<<std::endl;
  }

  return 0;

}


Int_t Fitter::improveMINOS(double NLL_min, int iPar, bool isErrHigh, int nGen, std::vector<double> *vLastHit, int seed)
{

  double weight = 0.01;
  double genWidth = 0.002;

  // Improve global fit result
  TRandom3 randGen (seed);

  std::vector<double> vTestPar(angPars.getSize());
  for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1)
    ((RooRealVar*)angPars.at(iPar1))->setVal(vLastHit->at(iPar1));

  double testNLL = 0;

  // std::cout<<nll->getValV();
  // for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1)
  //   std::cout<<"  "<<vLastHit->at(iPar1);
  auto nllVars = nll->getVariables();
  auto fitPars = result()->floatParsFinal();
  // printf(" S ");
  // for (int iPar1 = 0; iPar1 < fitPars.getSize(); ++iPar1)
  //   printf(" %f",((RooRealVar*)fitPars.at(iPar1))->getValV());
  // printf("\n");

  double intialVal = vLastHit->at(iPar);
  double roundVal = 1000*vLastHit->at(iPar)+0.5;
  if (roundVal<0) roundVal -= 1.;
  roundVal = ((int)roundVal)/1000.;	
  if (isErrHigh) { if ( roundVal+0.0005 > ((RooRealVar*)angPars.at(iPar))->getMax() ) return 0; }
  else if ( roundVal-0.0005 < ((RooRealVar*)angPars.at(iPar))->getMin() ) return 0;

  int iImprove = 0;

  do {

    for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1) {
      RooRealVar* par = (RooRealVar*)angPars.at(iPar1);
      do {
	if (iPar1==iPar) 
	  if (isErrHigh) vTestPar[iPar1] = randGen.Gaus(roundVal+0.0015,0.001);
	  else vTestPar[iPar1] = randGen.Gaus(roundVal-0.0015,0.001);
	else vTestPar[iPar1] = randGen.Gaus(vLastHit->at(iPar1),genWidth);
      } while (vTestPar[iPar1]>par->getMax() || 
	       vTestPar[iPar1]<par->getMin() || 
	       ( iPar1==iPar && ( ( isErrHigh  && vTestPar[iPar1]<roundVal+0.0005 ) ||
				  ( !isErrHigh && vTestPar[iPar1]>roundVal-0.0005 ) ) ) );
      par->setVal(vTestPar[iPar1]);
    }

    if (boundary->getValV()>0) continue;

    auto floatPars = result()->randomizePars();
    for (int iPar1 = 0; iPar1 < floatPars.getSize(); ++iPar1) {
      RooRealVar* par = (RooRealVar*)floatPars.at(iPar1);
      if (angPars.contains(*par)) continue; // angular parameters already set
      auto nllPar = (RooRealVar*)nllVars->find(*par);
      double fitPar = ((RooRealVar*)fitPars.find(*par))->getValV();
      if (nllPar) nllPar->setVal(weight*par->getValV()+(1.-weight)*fitPar);
    }

    testNLL = nll->getValV();
    if (testNLL-NLL_min<0.5) {
      // std::cout<<testNLL-NLL_min;
      for (int iPar1 = 0; iPar1 < angPars.getSize(); ++iPar1) {
	vLastHit->at(iPar1) = vTestPar[iPar1];
	// std::cout<<"  "<<vLastHit->at(iPar1);
      }
      // std::cout<<std::endl;
      // printf(" I ");
      // for (int iPar1 = 0; iPar1 < floatPars.getSize(); ++iPar1) {
      // 	RooRealVar* par = (RooRealVar*)floatPars.at(iPar1);
      // 	if (angPars.contains(*par)) continue;
      // 	printf(" %f",par->getValV());
      // }
      // printf("\n");

      roundVal = 1000*vLastHit->at(iPar)+0.5;
      if (roundVal<0) roundVal -= 1.;
      roundVal = ((int)roundVal)/1000.;	
      if (isErrHigh) { if ( roundVal+0.0005 > ((RooRealVar*)angPars.at(iPar))->getMax() ) break; }
      else if ( roundVal-0.0005 < ((RooRealVar*)angPars.at(iPar))->getMin() ) break;

    }

    ++iImprove;

  } while (iImprove<nGen);

  if ( vLastHit->at(iPar) != intialVal ) {
    if (isErrHigh) {
      vConfInterHigh[iPar] = vLastHit->at(iPar);
      std::cout<<((RooRealVar*)angPars.at(iPar))->GetName()<<" high impr: "<<vLastHit->at(iPar)<<std::endl;
    } else {
      vConfInterLow[iPar] = vLastHit->at(iPar);
      std::cout<<((RooRealVar*)angPars.at(iPar))->GetName()<<" low impr:  "<<vLastHit->at(iPar)<<std::endl;
    }
  }

  return 0;

}


void Fitter::doPlotMINOS(int iPar, std::vector<double> *vPval, std::vector<double> *vdNLL)
{

  TString parName  = ((RooRealVar*)angPars.at(iPar))->GetName();
  TString parTitle = ((RooRealVar*)angPars.at(iPar))->GetTitle();

  // produce deltaNLL vs parameter graph, with the probed points
  TCanvas* canNLL = new TCanvas(Form("canNLL_%s",parName.Data()),"canNLL",1000,1000);
  TGraph* grNLL = new TGraph(vPval->size(),&vPval->at(0),&vdNLL->at(0));
  grNLL->SetName(Form("grNLL_%s",parName.Data()));
  grNLL->SetTitle(Form("deltaNLL scan for %s",parTitle.Data()));
  grNLL->GetXaxis()->SetTitle(parName.Data());
  grNLL->GetYaxis()->SetTitle("deltaNLL");
  grNLL->SetMarkerStyle(7);
  // grNLL->SetMarkerStyle(20);
  // grNLL->SetMarkerSize(2);
  grNLL->SetMarkerColor(9);
  canNLL->cd();
  grNLL->Draw("AP");
	  
  TLine errLow (vConfInterLow[iPar],grNLL->GetYaxis()->GetXmin(),vConfInterLow[iPar],grNLL->GetYaxis()->GetXmax());
  TLine errHigh(vConfInterHigh[iPar],grNLL->GetYaxis()->GetXmin(),vConfInterHigh[iPar],grNLL->GetYaxis()->GetXmax());
  TLine DeltaNLL0p5 (grNLL->GetXaxis()->GetXmin(),0.5,grNLL->GetXaxis()->GetXmax(),0.5);
  errLow .SetLineColor(46);
  errHigh.SetLineColor(46);
  DeltaNLL0p5.SetLineColor(13);
  DeltaNLL0p5.SetLineStyle(9);
  errLow .Draw();
  errHigh.Draw();
  DeltaNLL0p5.Draw();

  canNLL->SaveAs(Form("%s/profiledNLL-%s_%s.pdf",plotDir.Data(),parName.Data(),plotString.Data()));

}


void Fitter::fillResultContainers(bool fromImprov)
{

  // fill results' containers

  for (int iPar = 0; iPar < angPars.getSize(); ++iPar) {

    RooRealVar* par = (RooRealVar*)angPars.at(iPar);
    vResult[iPar] = par->getValV();

    if (!fromImprov) {
      vFitResult[iPar] = vResult[iPar];
      vFitErrLow[iPar] = par->getErrorLo();
      vFitErrHigh[iPar] = par->getErrorHi();
    } else vImprovResult[iPar] = vResult[iPar];

  }

}


Double_t Fitter::computeBoundaryDistance()
{

  TStopwatch distTime;
  distTime.Start(true);

  // Compute distance from boundary
  boundDist = bound_dist->getValV();

  distTime.Stop();
  // if the result was not cached (producing time=~0) update the time
  if (distTime.CpuTime()>0.01)
    boundDistTime = distTime.CpuTime();

  std::cout<<"Distance from boundary: "<<boundDist<<" (computed in "<<boundDistTime<<" s)"<<std::endl;
  
  return boundDist;

}
