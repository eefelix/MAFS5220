#include <pch.h>
#include <ql/pricingengines/swap/analyticequityswapengine.hpp>

//! \file analyticequityswapengine.cpp

using namespace QuantLib;

void AnalyticESEngine::calculate() const
{
	QuantLib::Real leg1npv = 0;
	QuantLib::Real leg2npv = 0;

	for (size_t i = 0; i < arguments_.legs[0].size(); ++i) {
		leg1npv += arguments_.legs[0].at(i)->amount()*arguments_.discountCurve->discount(arguments_.legs[0].at(i)->date(), true);
		leg2npv += arguments_.legs[1].at(i)->amount()*arguments_.discountCurve->discount(arguments_.legs[1].at(i)->date(), true);
	}
	results_.legNPV.push_back(leg1npv);
	results_.legNPV.push_back(leg2npv);


	if (arguments_.type == EquitySwap::Type::Payer) {
		results_.value = arguments_.amount*(leg1npv - leg2npv);
	}
	else results_.value = arguments_.amount*(leg2npv - leg1npv);
}