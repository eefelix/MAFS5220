#include <pch.h>
#include "shortratetermstructure.hpp"
#include <ql/termstructures/yield/discountcurve.hpp>

//! \file shortratetermstructure.cpp
//!	
using namespace QuantLib;

ShortRateTermStructure::ShortRateTermStructure(
	const Calendar& calendar,
	const DayCounter& daycounter,
	const boost::shared_ptr<OneFactorAffineModel> &model) :
	calendar_(calendar), daycounter_(daycounter), model_(model)
{
}

ShortRateTermStructure::~ShortRateTermStructure()
{
}

boost::shared_ptr<YieldTermStructure>
ShortRateTermStructure::GetTermStructure(
Date referenceDate, Date maxDate, Real rate, Frequency freq/*= Quarterly*/) {
	std::vector<Date> dates;
	std::vector<Real> zcbondprices;
	Period tenor = Period();

	switch (freq)
	{
	case QuantLib::Annual:
		tenor = Period(12, Months);
		break;
	case QuantLib::Semiannual:
		tenor = Period(6, Months);
		break;
	case QuantLib::EveryFourthMonth:
		tenor = Period(4, Months);
		break;
	case QuantLib::Quarterly:
		tenor = Period(3, Months);
		break;
	case QuantLib::Bimonthly:
		tenor = Period(2, Months);
		break;
	case QuantLib::Monthly:
		tenor = Period(1, Months);
		break;
	case QuantLib::EveryFourthWeek:
		tenor = Period(28, Days);
		break;
	case QuantLib::Biweekly:
		tenor = Period(14, Days);
		break;
	case QuantLib::Weekly:
		tenor = Period(7, Days);
		break;
	case QuantLib::Daily:
		tenor = Period(1, Days);
		break;
	default:
		break;
	}

	int n = 0;
	while (referenceDate + n * tenor < maxDate) {
		dates.push_back(referenceDate + n * tenor);
		if (model_->params()[2] == QL_EPSILON){
			zcbondprices.push_back(std::exp(-rate*daycounter_.yearFraction(referenceDate, dates[n])));
		}
		else{
			zcbondprices.push_back(model_->discountBond(0,
				daycounter_.yearFraction(referenceDate, dates[n]), rate));
		}
		++n;
	}
	dates.push_back(maxDate);
	if (model_->params()[2] == QL_EPSILON){
		zcbondprices.push_back(std::exp(-rate*daycounter_.yearFraction(referenceDate, maxDate)));
	}
	else{
		zcbondprices.push_back(model_->discountBond(0, daycounter_.yearFraction(referenceDate, maxDate), rate));
	}

	return boost::make_shared<DiscountCurve>(dates, zcbondprices, daycounter_, calendar_);
}

boost::shared_ptr<YieldTermStructure>
ShortRateTermStructure::GetTermStructure(
Date referenceDate, Time Tmax, Real rate, Frequency freq/* = Quarterly*/) {
	Date maxDate = referenceDate+(int)Tmax*360;
	
	/*while (daycounter_.yearFraction(referenceDate, maxDate) < Tmax) {
		maxDate += 1;
	}*/
	return GetTermStructure(referenceDate, maxDate, rate, freq);
}