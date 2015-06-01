#ifndef intensitymodel_hpp
#define intensitymodel_hpp

#include <ql\types.hpp>
#include <ql\time\date.hpp>
#include "DefaultModel.hpp"

/**
* Class Name: IntensityModel.
* Usage:
*	1. Modeling the probability curve.
*	2. Currently support for piecewise constant hazard rate
*/
class IntensityModel : public DefaultModel {
public:
	IntensityModel()
	{
		dc = Actual360();
	}
	~IntensityModel(){}
	Probability SurvivalProbability(Date date);
	void Initialization(std::vector<Date> dateVector, std::vector<Rate> intensityVector);

private:
	Date startDate;
	std::vector<Real> segmentedDate;
	std::vector<Rate> segmentedIntensity;
	DayCounter dc;
};

#endif