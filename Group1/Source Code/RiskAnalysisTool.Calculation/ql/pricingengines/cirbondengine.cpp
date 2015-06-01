#include <pch.h>
#include "cirbondengine.hpp"

//! \file cirbondengine.cpp

using namespace QuantLib;

CIRBondEngine::CIRBondEngine(const boost::shared_ptr<CoxIngersollRoss>& model,
	const Handle<YieldTermStructure>& termStructure)
	:GenericModelEngine<CoxIngersollRoss, ZeroCouponBond::arguments, ZeroCouponBond::results>(model),
	termStructure_(termStructure)
{
}

void CIRBondEngine::calculate() const {
	Date maturityDate = arguments_.cashflows.back()->date();
	Date referenceDate = termStructure_->referenceDate();
	Real rate = model_->params()[3];
	results_.value = model_->discountBond(0, termStructure_->dayCounter().yearFraction(referenceDate, maturityDate), rate);
}