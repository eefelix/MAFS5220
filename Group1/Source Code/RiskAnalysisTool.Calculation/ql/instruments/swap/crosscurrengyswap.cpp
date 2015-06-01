#include <pch.h>
#include "crosscurrencyswap.hpp"
#include <ql/cashflows/fixedratecoupon.hpp>

//! \file crosscurrencyswap.cpp
namespace QuantLib {

	CrossCurrencySwap::CrossCurrencySwap(
		Type type,
		Rate fxRate,
		Volatility fxVol,
		const Date& referenceDate,
		Real payNominal,
		Real receiveNominal,
		const Schedule& paySchedule,
		Rate payRate,
		const DayCounter& payDayCount,
		const Schedule& receiveSchedule,
		Rate receiveRate,
		const DayCounter& receiveDayCount,
		boost::shared_ptr<Currency> payCurrency,
		boost::shared_ptr<Currency> receiveCurrency,
		boost::optional<BusinessDayConvention> paymentConvention /*= boost::none*/)
		:Swap(2), type_(type), fxRate_(fxRate), fxVol_(fxVol), referenceDate_(referenceDate), payNominal_(payNominal), receiveNominal_(receiveNominal),
		paySchedule_(paySchedule), payRate_(payRate), payDayCount_(payDayCount),
		receiveSchedule_(receiveSchedule), receiveRate_(receiveRate), receiveDayCount_(receiveDayCount),
		payCurrency_(payCurrency), receiveCurrency_(receiveCurrency) {

		if (paymentConvention)
			paymentConvention_ = *paymentConvention;
		else
			paymentConvention_ = receiveSchedule_.businessDayConvention();

		legs_[0] = FixedRateLeg(paySchedule_)
			.withNotionals(payNominal_)
			.withCouponRates(payRate_, payDayCount_)
			.withPaymentAdjustment(paymentConvention_);

		legs_[1] = FixedRateLeg(receiveSchedule_)
			.withNotionals(receiveNominal_)
			.withCouponRates(receiveRate_, receiveDayCount_)
			.withPaymentAdjustment(paymentConvention_);

	}

	CrossCurrencySwap::~CrossCurrencySwap() {
	}

	void CrossCurrencySwap::setupArguments(PricingEngine::arguments* args) const {
		Swap::setupArguments(args);

		CrossCurrencySwap::arguments* arguments =
			dynamic_cast<CrossCurrencySwap::arguments*>(args);

		if (!arguments)  // it's a swap engine...
			return;
		arguments->type = type_;
		arguments->fxRate = fxRate_;
		arguments->fxVol = fxVol_;
		arguments->referenceDate = referenceDate_;
		arguments->payNominal = payNominal_;
		arguments->receiveNominal = receiveNominal_;
		arguments->payDayCount = payDayCount_;
		arguments->receiveDayCount = receiveDayCount_;

		const Leg& payCoupons = payLeg();

		arguments->payCoupons = std::vector<Real>(payCoupons.size());
		arguments->payDates.resize(payCoupons.size());
		for (Size i = 0; i < payCoupons.size(); ++i) {
			boost::shared_ptr<FixedRateCoupon> coupon =
				boost::dynamic_pointer_cast<FixedRateCoupon>(payCoupons[i]);

			arguments->payDates[i] = coupon->date();
			arguments->payCoupons[i] = coupon->amount();
		}
		arguments->payCoupons.back() += payNominal_;
		const Leg& receiveCoupons = receiveLeg();

		arguments->receiveCoupons = std::vector<Real>(receiveCoupons.size());
		arguments->receiveDates.resize(receiveCoupons.size());
		for (Size i = 0; i < receiveCoupons.size(); ++i) {
			boost::shared_ptr<FixedRateCoupon> coupon =
				boost::dynamic_pointer_cast<FixedRateCoupon>(receiveCoupons[i]);

			arguments->receiveDates[i] = coupon->date();
			arguments->receiveCoupons[i] = coupon->amount();
		}
		arguments->receiveCoupons.back() += receiveNominal_;
	}

	void CrossCurrencySwap::arguments::validate() const {
		Swap::arguments::validate();
		QL_REQUIRE(payNominal != Null<Real>(), "pay nominal null or not set");
		QL_REQUIRE(receiveNominal != Null<Real>(), "receive nominal null or not set");
		QL_REQUIRE(payDates.size() == payCoupons.size(),
			"number of fixed payment dates different from "
			"number of fixed coupon amounts");
		QL_REQUIRE(receiveDates.size() == receiveCoupons.size(),
			"number of receive payment dates different from "
			"number of receive coupon amounts");
	}

	void CrossCurrencySwap::fetchResults(const PricingEngine::results* r) const {
		Swap::fetchResults(r);
	}
}