
#ifndef cdscalibrator_hpp
#define cdscalibrator_hpp

#include "AT1Pmodel.hpp"
#include <ql/types.hpp>
#include <ql/time/date.hpp>
#include <vector>
using namespace QuantLib;


/**
* Class Name: CDSCalibrateProblemFunction.
* Usage:
*	1. Calculate CDS value
*	2. Derived from CostFunction class to do the optimization.
*/
class CDSCalibrateProblemFunction : public CostFunction{

public:
	CDSCalibrateProblemFunction()
	{
		at1p = new AT1Pmodel();
		dc = Actual360();
		calendar = TARGET();
	}
	~CDSCalibrateProblemFunction(){ delete at1p; }

	void setRecoverRate(Real RecoveryRate) { this->recoveryRate = RecoveryRate; }
	void setTermStructure(std::vector<Date> dateVector) {
		CDSTermStructure.clear();
		CDSTermStructure.insert(CDSTermStructure.end(), dateVector.begin(), dateVector.end());
	}
	void setCdsQuotes(std::vector<Real> quoteVector) {
		termNum = quoteVector.size();
		CDSMarketQuote.clear();
		CDSMarketQuote.insert(CDSMarketQuote.end(), quoteVector.begin(), quoteVector.end());
	}
	void setImpliedSurvivalProb(std::vector<Probability> impliedProbs) {
		impliedProbability.clear();
		impliedProbability.insert(impliedProbability.end(), impliedProbs.begin(), impliedProbs.end());
	}

	/**
	* Function Name: value, values
	* Usage: Key functions for optimization class
	*/
	Real value(const Array& x) const;
	Disposable<Array> values(const Array& x) const;



private:

	AT1Pmodel* at1p;
	DayCounter dc;
	Calendar calendar;
	Rate recoveryRate;

	int termNum;
	std::vector<Date> CDSTermStructure;
	std::vector<Spread> CDSMarketQuote;
	std::vector<Probability> impliedProbability;

	//Real differentialization(Date testDate) const;
	//Real IntergratedCallFunction(Spread premium, Date currDate) const;
	//Real CDSValueFormula(Spread premium, Date startDate, Date endDate) const;
};




/**
* Class Name: CDSCalibrator.
* Usage: 
*	1. Use CDS market quote to calibrate implied piecewise constant volatility. 
*	2. Construct the survival probability curve.
*/
class CDSCalibrator {
public:
	CDSCalibrator(Real RecoveryRate, Real interestRate = 0.05){
		at1p = new AT1Pmodel(); 
		this->recoveryRate = RecoveryRate;
		this->interestRate = interestRate;
		optimizeFunc.setRecoverRate(RecoveryRate);
	}
	~CDSCalibrator() {
		delete at1p; 
	}

	void setStartDate(Date startDate) { this->startDate = startDate; }
	void setTermStructure(std::vector<Period> tenorVector) {
		Tenors.clear();
		Tenors.insert(Tenors.end(), tenorVector.begin(), tenorVector.end());
	}
	void setQuotes(std::vector<Real> quoteVector) {
		Quotes.clear();
		Quotes.insert(Quotes.end(), quoteVector.begin(), quoteVector.end());
	}
	/**
	* Function Name: extractDefaultProbabilityCurve
	* Usage: extract model independent survival probability using spreadCdsHelper
	*/
	void extractDefaultProbabilityCurve();

	/**
	* Function Name: StartCalibrate.
	* Return: Null
	* Usage:
	*	1. Use Optimizer to calibrate the optimal volatilities.
	*/
	void StartCalibrate();

	/* -------- Utility Function -------- */

	Probability FinalizedSurvivalProbability(Date testDate);
	Volatility FinalizedVolatility(Date testDate);


private:

	Date startDate;
	Real recoveryRate;
	Real interestRate;
	std::vector<Real> Quotes;
	std::vector<Period> Tenors;
	std::vector<Probability> impliedProbs;
	std::vector<Volatility> segmentedVol;

	AT1Pmodel* at1p;
	CDSCalibrateProblemFunction optimizeFunc;
	
};



#endif