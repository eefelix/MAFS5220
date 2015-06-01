#ifndef irsucvawithwwr_hpp
#define irsucvawithwwr_hpp

#include <ql\types.hpp>
#include <ql\time\date.hpp>
#include<vector>

using namespace QuantLib;

/**
* Class Name: CVAwithWWR.
* Usage:
*	1. Unilateral CVA calculator for interest rate swap
*	2. Considering the wrong way risk due to the correlation between hazard rate and short rate dynamics
*/
class CVAwithWWR{

public:
	//Receiver will receive fixed rate, while Payer will pay fixed rate
	enum Type { Receiver = -1, Payer = 1 };  
	CVAwithWWR(){};
	~CVAwithWWR(){};

	/**
		Utility functions to setup process parameters
	*/
	CVAwithWWR& WithShortRateCIRParam(Rate r0, Real theta, Real k, Real sigma1){
		this->r0 = r0; this->theta = theta; this->k = k; this->sigma1 = sigma1;
		return *this;
	}
	CVAwithWWR& WithIntensityCIRParam(Real lamda0, Real a, Real b, Real sigma2){
		this->lamda0 = lamda0; this->a = a; this->b = b; this->sigma2 = sigma2;
		return *this;
	}
	CVAwithWWR& WithCIRCorrelation(Real rho) { this->rho = rho; return *this; }
	CVAwithWWR& WithRecoveryRate(Real recoveryRate) { this->recoveryRate = recoveryRate; return *this; }

	CVAwithWWR& WithIRSSwap(Real nominal, Real spread, Real fixedRate, Date startDate, Frequency fixedFrequency,
	Frequency floatingFrequency, Period tenor, Type type){
		this->nominal = nominal; this->spread = spread; this->fixedRate = fixedRate; this->startDate = startDate;
		this->fixedFrequency = fixedFrequency; this->floatingFrequency = floatingFrequency; this->tenor = tenor; this->type = type;
		return *this;
	}

	CVAwithWWR& WithReferenceDate(Date referenceDate) { this->referenceDate = referenceDate; return *this; }

	void initialization();
	Real A(Date t1, Date t2);
	Real B(Date t1, Date t2);
	Real P(Date t1, Date t2, Real rt);
	/**
	* Function Name: bvaPath
	* Usage: calculate the net present value of vanilla swap
	*/
	Real npv(Date tao);
	Real ucvaPath(Date referenceDate, Date endDate);

	/**
	* Function Name: cvaWithWWRCalculation
	* Usage: Key entry for CVAWithWWR simulations
	*/
	Real cvaWithWWRCalculation(int pathNum);
	

private:
	// dr(t)=k*(theta-r(t))*dt+sigma1*sqrt(r(t))*dw(t)
	Real k;
	Real theta;
	Real sigma1;
	Real r0; //initial short rate

	// dlamda(t)=a*(b-lamda(t))*dt+sigma2*sqrt(lamda(t))*dw(t)
	Real a;
	Real b;
	Real sigma2;
	Real lamda0; //initial hazard rate

	Real rho; //correlation coefficient between short rate rt and hazard rate lamdat

	Real recoveryRate; //the recovery rate under counterparty default 
	Real nominal = 1.0;
	Real spread = 0.0; //the spread of floatingleg
	Rate fixedRate;
	Rate rtao; // short rate at default time tao
	Date referenceDate; //the date to calculate cva,dva.etc
	Date startDate; //the start date of swap contract
	Date endDate; //the last date of swap contract, which is also the last date to exchange cashflow 
	Frequency fixedFrequency;
	Frequency floatingFrequency;
	Period tenor;
	Type type;

	std::vector<Date> floatingLeg; //the vector that stores the cashflow exchange date of floatingLeg
	std::vector<Date> fixedLeg; //the vector that stores the cashflow exchange date of fixedLeg
	
	const DayCounter dc = Actual360();
};

#endif