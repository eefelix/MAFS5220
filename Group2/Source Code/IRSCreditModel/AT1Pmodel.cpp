
#include "pch.h"
#include "AT1Pmodel.hpp"
#include <boost/math/distributions.hpp>

using namespace QuantLib;

void AT1Pmodel::Initialization(std::vector<Date> dateVector, std::vector<Volatility> volVector)
{
	startDate = dateVector.front();
	segmentedDate.clear();
	for (std::vector<Date>::const_iterator iterDate = dateVector.cbegin();
		iterDate != dateVector.cend(); iterDate++) 
	{
		segmentedDate.push_back(dc.yearFraction(startDate, *iterDate));
	}

	segmentedVol.clear();
	segmentedVol.insert(segmentedVol.end(), volVector.begin(), volVector.end());
}


Probability AT1Pmodel::at1pFormula(Volatility cumulativeVariance) {
	boost::math::normal_distribution<> d(0, 1);
	Real firstTerm = cdf(d, (log(1 / BarrierFirmvalueRatio) + (BarrierExponentialVolParameter - 0.5) * cumulativeVariance) / sqrt(cumulativeVariance));
	Real secondTerm = cdf(d, (log(BarrierFirmvalueRatio) + (BarrierExponentialVolParameter - 0.5) * cumulativeVariance) / sqrt(cumulativeVariance));
	Probability tmpProb = firstTerm - pow(BarrierFirmvalueRatio, 2 * BarrierExponentialVolParameter - 1) * secondTerm;
	return tmpProb;
}

Probability AT1Pmodel::SurvivalProbability(Date date) {
	Real testDate = dc.yearFraction(startDate, date);
	QL_REQUIRE(testDate >= segmentedDate.front() && testDate <= segmentedDate.back(),
		"Current time is not in the range, please re-enter the valid date.");
	Volatility cumulativeVariance = 0.0;
	std::vector<Volatility>::iterator iterVol = segmentedVol.begin();
	Real testDatePeriod = testDate;
	for (std::vector<Real>::iterator iterDatePeriod = ++segmentedDate.begin();
		iterDatePeriod != segmentedDate.end(); ++iterDatePeriod) {
		if (testDatePeriod <= *iterDatePeriod) {
			cumulativeVariance += pow((*iterVol), 2) * (testDatePeriod - *(iterDatePeriod - 1));
			break;
		}
		else {
			cumulativeVariance += pow((*iterVol), 2) * (*iterDatePeriod - *(iterDatePeriod - 1));
			iterVol++;
		}
	}
	Probability survivalProb = at1pFormula(cumulativeVariance);
	return survivalProb;
}
