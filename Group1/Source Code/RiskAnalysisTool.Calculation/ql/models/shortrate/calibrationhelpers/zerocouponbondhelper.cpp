#include <pch.h>
#include <ql/pricingengines/bond/discountingbondengine.hpp>
#include "zerocouponbondhelper.hpp"

//! \file zerocouponbondhelper.cpp
//!
using namespace QuantLib;

ZerocouponbondHelper::ZerocouponbondHelper(
	Natural settlementDays,
	const Calendar& calendar,
	Real faceAmount,
	const Date& maturityDate,
	BusinessDayConvention paymentConvention,
	Real redemption,
	const Date& issueDate,
	const Handle<Quote>& volatility,
	const Handle<YieldTermStructure>& termStructure,
	CalibrationErrorType errorType) :
	zcbond_(new ZeroCouponBond(settlementDays,calendar,faceAmount,maturityDate,paymentConvention,redemption,issueDate)),
	CalibrationHelper(volatility, termStructure, errorType)
{
}

void ZerocouponbondHelper::addTimesTo(std::list<Time>& times) const {
	//this is to be used with tree-based methods 
	QL_FAIL("not implement!");
}

Real ZerocouponbondHelper::modelValue() const {
	calculate();
	zcbond_->setPricingEngine(engine_);
	return zcbond_->NPV();
}

Real ZerocouponbondHelper::blackPrice(Volatility sigma) const {
	//market data
	calculate();
	boost::shared_ptr<PricingEngine> marketengine(new DiscountingBondEngine(termStructure_));
	zcbond_->setPricingEngine(marketengine);
	Real value = zcbond_->NPV();
	zcbond_->setPricingEngine(engine_);
	return value;
}

void ZerocouponbondHelper::performCalculations() const {
	CalibrationHelper::performCalculations();
}