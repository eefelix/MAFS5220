#include <pch.h>
#include "mcvanillaswapexposuremodel.hpp"

//! \file mcvanillaswapexposuremodel.cpp
//!
using namespace QuantLib;

Real MCVanillaSwapExposureModel::exposure(const MultiPath& path, const Time issuerDefaultTime,
	const Time investorDefaultTime, const Handle<YieldTermStructure> disTS) const {
	// find default time
	Time defaultTime = std::min(issuerDefaultTime, investorDefaultTime);

	boost::shared_ptr<const VanillaSwap> inst = boost::dynamic_pointer_cast<const VanillaSwap>(instrument_);
	Date defaultDate = disTS->referenceDate();
	std::vector<Date> fixedDates = arguments_.fixedPayDates;
	std::vector<Date> floatingDates = arguments_.floatingPayDates;
	// transfer payment Dates to Times
	std::vector<Real> fixedTimes;
	std::vector<Real> floatingTimes;
	for (int i = 0; i != fixedDates.size(); ++i){
		fixedTimes.push_back(inst->fixedDayCount().yearFraction(defaultDate, fixedDates[i]));
	}
	for (int i = 0; i != floatingDates.size(); ++i){
		floatingTimes.push_back(inst->floatingDayCount().yearFraction(defaultDate, floatingDates[i]));
	}
	// fixed discount factor P(defaultTime, Ti) for fixed leg
	// floating discount factor P(defaultTime, Ti) for floating leg
	std::vector<Real> fixedDiscountFactor;
	std::vector<Real> floatingDiscountFactor;
	int fixedIndex = 0;
	int floatingIndex = 0;

	// calculate the corresponding discount factors of each leg to default time
	for (int i = 0; i != fixedTimes.size(); ++i) {
		if (fixedTimes[i] >= 0.0) {
			fixedDiscountFactor.push_back(disTS->discount(fixedTimes[i], true));
		}
		else {
			++fixedIndex;
		}
	}
	for (int i = 0; i != floatingTimes.size(); ++i) {
		if (floatingTimes[i] >= 0.0) {
			floatingDiscountFactor.push_back(disTS->discount(floatingTimes[i], true));
		}
		else {
			++floatingIndex;
		}
	}
	//calculate the NPV(default)
	Real NPVdefault = 0.0;
	NPVdefault += arguments_.nominal * (1. - floatingDiscountFactor.back());
	for (int m = 0; m != floatingDiscountFactor.size(); ++m){
		NPVdefault += arguments_.nominal * arguments_.floatingSpreads[m + floatingIndex] * floatingDiscountFactor[m];
	}
	for (int m = 0; m != fixedDiscountFactor.size(); ++m){
		NPVdefault -= arguments_.fixedCoupons[m + fixedIndex] * fixedDiscountFactor[m];
	}

	if (arguments_.type == QuantLib::VanillaSwap::Payer){
		NPVdefault = NPVdefault;
	}
	else{
		NPVdefault = -NPVdefault;
	}
	return NPVdefault;
}