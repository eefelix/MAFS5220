
#ifndef defaultcurve_hpp
#define defaultcurve_hpp

#include <ql\types.hpp>
#include <ql\time\date.hpp>
#include <vector>
using namespace QuantLib;

/**
* Class Name: DefaultCurve.
* Usage:
*	1. Base class for different types of default profiles
	2. user-defined two enum type of @{ConstantIntensityGiven} and @{DefaultProbabilityGiven}
	3. Mainly used in UCVA calculator
*/
class DefaultCurve :public Observable{
public:
	enum DefaultCurveType {
		ConstantIntensityGiven,
		DefaultProbabilityGiven,
	};

	DefaultCurve(DefaultCurveType type, 
			std::vector<Date> creditTerm, 
			std::vector<Probability> survivalProb);
	DefaultCurve(DefaultCurveType type,
			Date startDate,
			Rate ConstantIntensity);
	~DefaultCurve(){}

	DefaultCurveType getType() { return this->type; }
	Date getStartDate() { return startDate; }
	std::vector<Date> getDefaultTerm() { return cpTermDate; }
	Probability getSurvivalProb(size_t index) { return cpSurvivalProb[index]; }
	void buildCurve(Period maturity, Period interval);

private:
	DefaultCurveType type;
	DayCounter dc;
	Date startDate;
	Date maturityDate;
	Rate Intensity;
	std::vector<Probability> cpSurvivalProb;
	std::vector<Real> cpTerm;
	std::vector<Date> cpTermDate;
};

#endif