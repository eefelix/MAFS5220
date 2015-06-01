
#ifndef at1pmodel_hpp
#define at1pmodel_hpp

#include <ql/types.hpp>
#include <ql/time/date.hpp>
#include <vector>
#include "DefaultModel.hpp"
using namespace QuantLib;


/**
* Class Name: AT1Pmodel.
* Usage:
*	1. Modeling the probability curve.
*/
class AT1Pmodel : public DefaultModel {
public:
	AT1Pmodel(Real B = 0.0, Real Ratio = 0.4)
		:BarrierExponentialVolParameter(B), BarrierFirmvalueRatio(Ratio)
	{ dc = Actual360();	}
	~AT1Pmodel(){}

	/**
	* Function Name: SurvivalProbability
	* Usage: get survival probability under AT1P model
	*/
	Probability SurvivalProbability(Date date);
	/**
	* Function Name: Initialization
	* Usage: setup model parameters for AT1P model
	*/
	void Initialization(std::vector<Date> dateVector, std::vector<Volatility> volVector);

private:
	Real BarrierExponentialVolParameter;
	Real BarrierFirmvalueRatio;
	Date startDate;
	std::vector<Real> segmentedDate;
	std::vector<Volatility> segmentedVol;
	DayCounter dc;

	Probability at1pFormula(Volatility cumulativeVariance);
};


#endif