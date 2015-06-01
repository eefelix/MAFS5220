#include <pch.h>
#include "mcequityswapexposuremodel.hpp"
#include <ql/pricingengines/swap/analyticequityswapengine.hpp>
#include <ql/instruments/swap/equityswap.hpp>

//! \file mcequityswapexposuremodel.cpp
//!
using namespace QuantLib;

Real MCEquitySwapExposureModel::exposure(const MultiPath& path, const Time issuerDefaultTime,
	const Time investorDefaultTime, const Handle<YieldTermStructure> disTS) const {

	Time defaultTime = std::min(issuerDefaultTime, investorDefaultTime);

	TimeGrid timeGrid = path[1].timeGrid();
	std::vector<Time>::const_iterator defaultPos = std::find(timeGrid.begin(), timeGrid.end(), defaultTime);
	double spotprice = path[1][defaultPos - timeGrid.begin()];

	double cumulativedividend = 0;
	double cumulativecoupon = 0;

	size_t defaultgrid = std::min<size_t>(timeGrid.closestIndex(defaultTime), timeGrid.size() - 2);
	for (std::vector<double>::reverse_iterator it = arguments_->yearfraction.rbegin();
		it != arguments_->yearfraction.rend(); ++it){
		if (*it < defaultTime){
			cumulativecoupon = (defaultTime - *it)
				*arguments_->fixedRate*arguments_->startPrice;
			for (size_t j = timeGrid.closestIndex(*it); j <= defaultgrid; ++j){
				cumulativedividend += path[1][j] * arguments_->dividend*timeGrid.dt(j);
			}
			if (it == (arguments_->yearfraction.rend() - 1)){
				cumulativecoupon += arguments_->cumulativeCoupoon;
				cumulativedividend += arguments_->cumulativeDividend;
			}
			break;
		}
	}
	// calculate default date
	Date defaultDate = disTS->referenceDate();

	EquitySwap ers(arguments_->type, arguments_->startDate, defaultDate,
		arguments_->startPrice, spotprice, 1,
		cumulativecoupon, cumulativedividend,
		arguments_->fixedRate, arguments_->dividend, arguments_->maturity, arguments_->sigma,
		*disTS);
	/* the risk free rate should be adapted to the time level tau*/
	ers.setPricingEngine(boost::shared_ptr<AnalyticESEngine>(new AnalyticESEngine()));
	// cva is the positive part of NPV
	return ers.NPV();

}