#include "pch.h"
#include "DefaultCurve.hpp"

using namespace QuantLib;

DefaultCurve::DefaultCurve(DefaultCurveType type, 
						std::vector<Date> creditTerm, 
						std::vector<Probability> survivalProb)
{
	this->type = type;
	this->dc = Actual360();
	cpSurvivalProb.clear();
	cpSurvivalProb.insert(cpSurvivalProb.end(), survivalProb.begin(), survivalProb.end());
	cpTermDate.clear();
	cpTermDate.insert(cpTermDate.end(), creditTerm.begin(), creditTerm.end());
	startDate = creditTerm.front();
	cpTerm.clear();
	for (size_t i = 0; i < creditTerm.size(); i++){
		cpTerm.push_back(dc.yearFraction(startDate, creditTerm[i]));
	}
}

DefaultCurve::DefaultCurve(DefaultCurveType type,
						Date startDate,
						Rate ConstantIntensity)
{
	this->type = type;
	this->dc = Actual360();
	this->startDate = startDate;
	this->Intensity = ConstantIntensity;
}

void DefaultCurve::buildCurve(Period maturity, Period interval)
{
	maturityDate = startDate + maturity;
	if (type == DefaultProbabilityGiven){
		return;
	}
	else if (type == ConstantIntensityGiven) {
		//Initialize Default Model
		cpTerm.clear();
		cpTermDate.clear();
		Date tmpDate = startDate;
		do {
			cpTerm.push_back(dc.yearFraction(startDate, tmpDate));
			cpTermDate.push_back(tmpDate);
			tmpDate += interval;
		} while (tmpDate <= maturityDate);

		cpSurvivalProb.clear();
		for (std::vector<Real>::const_iterator termIter = cpTerm.cbegin();
			termIter != cpTerm.cend(); termIter++)
		{
			cpSurvivalProb.push_back(exp(-Intensity * (*termIter)));
		}
	}
	else {
		//do nothing
	}
}