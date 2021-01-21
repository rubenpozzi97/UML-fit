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
	       Penalty* _penTerm) :
  name(_name),
  title(_title),
  angPars(_angPars),
  combData(_combData),
  simPdf(_simPdf),
  simPdf_penalty(_simPdf_penalty),
  boundary(_boundary),
  bound_dist(_bound_dist),
  penTerm(_penTerm)
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
	       BoundDist* _bound_dist) :
  name(_name),
  title(_title),
  angPars(_angPars),
  combData(_combData),
  simPdf(_simPdf),
  simPdf_penalty(_simPdf_penalty),
  boundary(_boundary),
  bound_dist(_bound_dist)
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
  penTerm(other.penTerm)
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
      boundDist = bound_dist->getValV();
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
	    boundDist = bound_dist->getValV();
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
  
  // Improve global fit result
  TRandom3 randGen (1);

  std::vector<double> vTestPar(angPars.getSize());
  std::vector<double> vImprovPar(angPars.getSize());
  for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
    vImprovPar[iPar] = ((RooRealVar*)angPars.at(iPar))->getValV();

  double NLL_before = nll->getValV();
  double improvNLL = NLL_before;
  double testNLL = 0;
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
      for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
	vImprovPar[iPar] = vTestPar[iPar];
    }

    ++iImprove;

  } while (iImprove<nGen);

  for (int iPar = 0; iPar < angPars.getSize(); ++iPar)
    ((RooRealVar*)angPars.at(iPar))->setVal(vImprovPar[iPar]);
  
  boundDist = bound_dist->getValV();
  std::cout<<"Improved fit result: deltaNLL = "<<NLL_before-improvNLL<<" bound dist: "<<preBoundDist<<" -> "<<boundDist<<std::endl;

  return 0;

}
