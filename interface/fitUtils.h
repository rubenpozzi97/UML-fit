#include <RooFitResult.h>
#include <RooDataSet.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooMinimizer.h>

#include "BoundCheck.h"
#include "Penalty.h"


RooFitResult* fit (RooDataSet* combData,
		   RooAbsPdf* simPdf,
		   RooAbsPdf* simPdf_penalty,
		   RooAbsReal* & nll,
		   RooAbsReal* & nll_penalty,
		   BoundCheck* boundary,
		   Penalty* penTerm,
		   double fac1,
		   double fac4,
		   double base1,
		   double base4,
		   double max1,
		   double max4)
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
    
    RooFitResult* fitResult = m.save("result") ;

    RooFitResult* fitResult_penalty = 0;
    usedPenalty = false;

    // if free fit is good return its result
    if ( fitResult->status()==0 && fitResult->covQual()==3 && boundary->getValV() == 0 ) return fitResult;

    usedPenalty = true;

    // optional: if a partial boundary is satisfied
    // do not apply the corresponding penalty term
    // bool inCTL4  = true;
    // bool inCTL15 = true;
    // if ( !boundary->isInCTL4()  ) inCTL4  = false;
    // if ( !boundary->isInCTL15() ) inCTL15 = false;

    for (int totCoeff=0; fac1*pow(base1,totCoeff)<=maxCoeff; ++totCoeff) {

      for (int iCoeff1=totCoeff; iCoeff1>=0; --iCoeff1) {

	// set penalty coefficients
	coeff1 = fac1 * pow(base1,iCoeff1);
	if (max1>0 && coeff1>max1) continue;

	coeff4 = fac4 * pow(base4,totCoeff-iCoeff1);
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
	fitResult_penalty = m_penalty.save("result");
	    
	// cout<<penTerm->getCoefficient(1)<<"\t"<<penTerm->getCoefficient(5)<<"\t"<<P5p->getValV()<<endl;
	// fitResult_penalty->Print("v");

	// if a good fit is found return its result
	if ( fitResult_penalty->status()==0 && fitResult_penalty->covQual()==3 ) {
	  if ( boundary->getValV()==0 ) {
	    // cout<<"P "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
	    return fitResult_penalty;
	  } // else cout<<"O "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;
	} // else cout<<"N "<<coeff1<<"\t"<<coeff4<<"\t"<<coeff5<<endl;

      }

    }
    
    // if no good fit is found return a null pointer
    return 0;

}
