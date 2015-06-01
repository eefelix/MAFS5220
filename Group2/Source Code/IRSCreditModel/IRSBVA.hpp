#ifndef irsbva_hpp
#define irsbva_hpp

#include <ql\types.hpp>
#include <ql\time\date.hpp>

using namespace QuantLib;

/**
* Class Name: IRSBVA.
* Usage:
*	1. Bilateral-VA calculator for interest rate swap
*	2. Use Monte Carlo simulation
*/
class IRSBVA{

public:
	enum Type { Receiver = -1, Payer = 1 };
	IRSBVA(){};
	~IRSBVA(){};

	/**
		Utility functions to setup process parameters
		*/

	IRSBVA& WithShortRateCIRParam(Rate r0, Real ar, Real br, Real sigmar){
		this->r0 = r0; this->ar = ar; this->br = br; this->sigmar = sigmar;
		return *this;
	}

	IRSBVA& WithBIntensityCIRParam(Real lamda0b, Real ab, Real bb, Real sigmab){
		this->lamda0b = lamda0b; this->ab = ab; this->bb = bb; this->sigmab = sigmab;
		return *this;
	}

	IRSBVA& WithCIntensityCIRParam(Real lamda0c, Real ac, Real bc, Real sigmac){
		this->lamda0c = lamda0c; this->ac = ac; this->bc = bc; this->sigmac = sigmac;
		return *this;
	}

	IRSBVA& WithCIRCorrelation(Real rhobc, Real rhobr, Real rhocr) {
		this->rhobc = rhobc; this->rhobr = rhobr; this->rhocr = rhocr;
		return *this; 
	}

	IRSBVA& WithRecoveryRate(Real recoveryRateb, Real recoveryRatec) {
		this->recoveryRateb = recoveryRateb; 
		this->recoveryRatec = recoveryRatec;
		return *this; 
	}
	
	IRSBVA& WithIRSSwap(Real nominal, Real spread, Real fixedRate, Date startDate, Frequency fixedFrequency,
		Frequency floatingFrequency, Period tenor, Type type){
		this->nominal = nominal; this->spread = spread; this->fixedRate = fixedRate; this->startDate = startDate;
		this->fixedFrequency = fixedFrequency; this->floatingFrequency = floatingFrequency; this->tenor = tenor; this->type = type;
		return *this;
	}

	IRSBVA& WithReferenceDate(Date referenceDate) { this->referenceDate = referenceDate; return *this; }

	void initialization();
	Real A(Date t1, Date t2);
	Real B(Date t1, Date t2);
	Real P(Date t1, Date t2, Real rt);

	/**
	* Function Name: bvaPath
	* Usage: calculate the net present value of vanilla swap
	*/
	Real npv(Date tao);
	/**
	* Function Name: bvaPath
	* Usage: Single simulation path
	*/
	Real bvaPath(Date referenceDate, Date endDate);

	/**
	* Function Name: bvaCalculation
	* Usage: Key entry for bva simulations
	*/
	Real bvaCalculation(int pathNum);

private:
	/**
		Most of the parameters are similar with the setting of UCVA WWR simulation
		Please kindly refer to CVAWithWWR.hpp
	*/
	Real ar;
	Real br;
	Real sigmar;
	Real r0;

	Real ab;
	Real bb;
	Real sigmab;
	Real lamda0b;

	Real ac;
	Real bc;
	Real sigmac;
	Real lamda0c;

	Real rhobc;
	Real rhobr;
	Real rhocr;

	Real recoveryRateb;
	Real recoveryRatec;

	Real nominal = 1.0;
	Real spread = 0.0;
	Rate fixedRate;
	Rate rtao;
	Date referenceDate;
	Date startDate;
	Date endDate;
	Frequency fixedFrequency;
	Frequency floatingFrequency;
	Period tenor;
	Type type;
	std::vector<Date> floatingLeg;
	std::vector<Date> fixedLeg;

	const DayCounter dc = Actual360();
};

#endif