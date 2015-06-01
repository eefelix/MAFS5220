#include <pch.h>
#include "equityswap.hpp"
#include <ql/cashflows/simplecashflow.hpp>

//! \file equityswap.cpp

using namespace QuantLib;


EquitySwap::EquitySwap(
	Type type,
	QuantLib::Date startdate,
	QuantLib::Date referencedate,
	QuantLib::Real startprice,
	QuantLib::Real spotprice,
	QuantLib::Integer amount,
	QuantLib::Real cumulativecoupoon,
	QuantLib::Real cumulativedividend,
	QuantLib::Rate fixedrate,
	QuantLib::Rate dividend,
	QuantLib::Time maturity,
	QuantLib::Volatility sigma,
	boost::shared_ptr<QuantLib::YieldTermStructure> discountcurve,
	const QuantLib::Frequency& freq,
	const QuantLib::BusinessDayConvention& busDayConvention,
	const QuantLib::DateGeneration::Rule& rule,
	const QuantLib::Calendar& calendar) :
	QuantLib::Swap(2),
	type_(type), startPrice_(startprice), spotPrice_(spotprice), amount_(amount), referenceDate_(referencedate), startDate_(startdate),
	cumulativeCoupoon_(cumulativecoupoon), fixedRate_(fixedrate), dividend_(dividend),
	maturity_(maturity), cumulativeDividend_(cumulativedividend), discountCurve_(discountcurve),
	sigma_(sigma), frequency_(freq), BusDayConvention_(busDayConvention), rule_(rule), calendar_(calendar)
{
	legNPV_.reserve(2);
	initialize();
}

EquitySwap::EquitySwap(const EquitySwap& ers) :QuantLib::Swap(2){
	type_ = ers.type_;
	referenceDate_ = ers.referenceDate_;
	startDate_ = ers.startDate_;
	startPrice_ = ers.startPrice_;
	spotPrice_ = ers.spotPrice_;
	amount_ = ers.amount_;
	cumulativeDividend_ = ers.cumulativeDividend_;
	cumulativeCoupoon_ = ers.cumulativeCoupoon_;
	fixedRate_ = ers.fixedRate_;
	dividend_ = ers.dividend_;
	maturity_ = ers.maturity_;
	sigma_ = ers.sigma_;
	discountCurve_ = ers.discountCurve_;
	frequency_ = ers.frequency_;
	BusDayConvention_ = ers.BusDayConvention_;
	rule_ = ers.rule_;

	initialize();
}

void EquitySwap::initialize()
{
	int n = 0;
	int paytimes = 0;

	// the payment frequency is supposed to be
	// 1 year, 6 month, 3 month

	switch (frequency_)
	{
	case Annual:
		n = 12;
		paytimes = static_cast<int>(std::floor(maturity_));
		break;
	case Semiannual:
		n = 6;
		paytimes = static_cast<int>(std::floor(maturity_ * 2));
		break;
	case Quarterly:
		n = 3;
		paytimes = static_cast<int>(std::floor(maturity_ * 4));
		break;
	default:
		break;
	}

	// for the time being the interest rate is supposed constant
	double constforwardrate = discountCurve_->forwardRate(0, 0.00001, QuantLib::Continuous);

	// make the complete schedule for the equity swap, from the start date of the contract
	// to the end date of the contract.
	schedule_ = QuantLib::Schedule(
		startDate_, startDate_ + ((int)maturity_) * Years, n * Months,
		calendar_, ModifiedFollowing, ModifiedFollowing, DateGeneration::Forward, false);

	yearfraction_.push_back(0);
	std::vector<QuantLib::Real> expectedstockprice;

	// find the first payment day after the reference date
	// the date should not excess the last payment date
	std::vector<Date>::const_iterator it = schedule_.dates().begin();
	while (*it <= referenceDate_ && (it != schedule_.dates().end() - 1)){
		++it;
	}
	const std::vector<Date>::const_iterator validit = it;

	do{
		yearfraction_.push_back(discountCurve_->timeFromReference(*it));
		expectedstockprice.push_back(spotPrice_*exp(-dividend_*yearfraction_[it - validit + 1])
			/ discountCurve_->discount(yearfraction_[it - validit + 1], true));
		++it;
	} while (it != schedule_.dates().end());

	count_ = yearfraction_.size() - 1; // the number of payments to make

	// if there is only one payment
	// the cumulative coupon and dividend and the notional payments
	// should be taken into consideration in the same time

	if (count_ == 1 || (it == schedule_.dates().end() - 1)){
		legs_[1].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
			(new QuantLib::SimpleCashFlow(
			cumulativeCoupoon_ + fixedRate_*startPrice_*yearfraction_[1] + startPrice_, schedule_.endDate())));
		legs_[0].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
			(new QuantLib::SimpleCashFlow(
			cumulativeDividend_ + dividend_*(expectedstockprice[0] - spotPrice_) / (constforwardrate - dividend_)
			+ expectedstockprice[0], schedule_.endDate())));
	}

	// if there is more than one payments
	// the first payment should consider the cumulative coupon and dividend
	// the last payment should consider the exchange of notional

	else{
		legs_[1].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
			(new QuantLib::SimpleCashFlow(cumulativeCoupoon_ + fixedRate_*startPrice_*yearfraction_[1], *validit)));
		legs_[0].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
			(new QuantLib::SimpleCashFlow(cumulativeDividend_ +
			dividend_*(expectedstockprice[0] - spotPrice_) / (constforwardrate - dividend_), *validit)));
		for (std::vector<Date>::const_iterator iter = validit + 1; iter != schedule_.end() - 1; ++iter){
			size_t m = iter - validit;
			legs_[1].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
				(new QuantLib::SimpleCashFlow(
				startPrice_*fixedRate_*(yearfraction_[m + 1] - yearfraction_[m]), *iter)));
			legs_[0].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
				(new QuantLib::SimpleCashFlow(dividend_*(expectedstockprice[m] - expectedstockprice[m - 1])
				/ (constforwardrate - dividend_), *iter)));
		}
		// the final leg: exchange the notional with stocks
		legs_[1].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
			(new QuantLib::SimpleCashFlow(
			(yearfraction_[count_] - yearfraction_[count_ - 1])*startPrice_*fixedRate_ + startPrice_, schedule_.endDate())));
		legs_[0].push_back(boost::shared_ptr<QuantLib::SimpleCashFlow>
			(new QuantLib::SimpleCashFlow(
			dividend_*(expectedstockprice[count_ - 1] - expectedstockprice[count_ - 2]) / (constforwardrate - dividend_)
			+ expectedstockprice[count_ - 1], schedule_.endDate())));
	}
}

void EquitySwap::setupArguments(QuantLib::PricingEngine::arguments* args) const
{
	EquitySwap::arguments* arguments =
		dynamic_cast<EquitySwap::arguments*>(args);
	QL_REQUIRE(arguments != 0, "wrong argument type");
	arguments->count = count_;
	arguments->discountCurve = discountCurve_;
	arguments->spotPrice = spotPrice_;
	arguments->startPrice = startPrice_;
	arguments->cumulativeCoupoon = cumulativeCoupoon_;
	arguments->cumulativeDividend = cumulativeDividend_;
	arguments->amount = amount_;
	arguments->referenceDate = referenceDate_;
	arguments->startDate = startDate_;
	arguments->type = type_;
	arguments->fixedRate = fixedRate_;
	arguments->dividend = dividend_;
	arguments->sigma = sigma_;
	arguments->maturity = maturity_;
	arguments->freq = frequency_;
	arguments->yearfraction = yearfraction_;
	arguments->schedule = boost::make_shared<QuantLib::Schedule>(schedule_);
	// original swap arguments
	arguments->legs = legs_;
	arguments->payer = payer_;
}

void EquitySwap::fetchResults(const QuantLib::PricingEngine::results* r) const {
	Instrument::fetchResults(r);

	const EquitySwap::results* results = dynamic_cast<const EquitySwap::results*>(r);
	QL_REQUIRE(results != 0, "wrong result type");
	if (!results->legNPV.empty()) {
		QL_REQUIRE(results->legNPV.size() == legNPV_.size(),
			"wrong number of leg NPV returned");
		legNPV_ = results->legNPV;
	}
	else {
		std::fill(legNPV_.begin(), legNPV_.end(), QuantLib::Null<double>());
	}
}