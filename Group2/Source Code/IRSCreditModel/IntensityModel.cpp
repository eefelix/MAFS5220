
#include "pch.h"
#include "IntensityModel.hpp"
#include <math.h>

using namespace QuantLib;

void IntensityModel::Initialization(std::vector<Date> dateVector, std::vector<Rate> intensityVector)
{
	startDate = dateVector.front();
	segmentedDate.clear();
	for (std::vector<Date>::const_iterator iterDate = dateVector.cbegin();
		iterDate != dateVector.cend(); iterDate++)
	{
		segmentedDate.push_back(dc.yearFraction(startDate, *iterDate));
	}
	segmentedIntensity.clear();
	segmentedIntensity.insert(segmentedIntensity.end(), intensityVector.begin(), intensityVector.end());
}

Probability IntensityModel::SurvivalProbability(Date date)
{
	Real testDate = dc.yearFraction(startDate, date);
	QL_REQUIRE(testDate >= segmentedDate.front() && testDate <= segmentedDate.back(),
		"Current time is not in the range, please re-enter the valid date.");
	Real cumulativeIntensity = 0.0;
	std::vector<Rate>::const_iterator iterIntensity = segmentedIntensity.cbegin();
	for (std::vector<Real>::const_iterator iterDate = ++segmentedDate.cbegin();
		iterDate != segmentedDate.end(); iterDate++) {
		if (testDate <= *iterDate) {
			cumulativeIntensity += *iterIntensity * (testDate - *(iterDate - 1)); 
			break;
		}
		else {
			cumulativeIntensity += *iterIntensity * (*iterDate - *(iterDate - 1));
			iterIntensity++;
		}
	}
	return exp(-cumulativeIntensity);
}