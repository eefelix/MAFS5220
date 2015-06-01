#include <pch.h>
#include "mccrosscurrencyswapexposuremodel.hpp"
#include <ql\math\randomnumbers\seedgenerator.hpp>
#include <ql\math\randomnumbers\boxmullergaussianrng.hpp>

//! \file mccrosscurrencyswapexposuremodel.cpp
//!
using namespace QuantLib;

Real MCCrossCurrencySwapExposureModel::exposure(const MultiPath& path, const Time issuerDefaultTime,
	const Time investorDefaultTime, const Handle<YieldTermStructure> disTS) const {

	// find fx rate at default time
	Time defaultTime = std::min(issuerDefaultTime, investorDefaultTime);
	TimeGrid tg = path[0].timeGrid();
	Time dt = tg[1] - tg[0];
	double domesticRate = path[0][tg.index(defaultTime)];
	double foreignRate = path[1][tg.index(defaultTime)];
	double fxRate = arguments_.fxRate;
	double fxVol = arguments_.fxVol;
	fxRate *= std::exp(-0.5*fxVol*fxVol*(defaultTime));
	for (int i = 1; i <= defaultTime; ++i) {
		fxRate *= std::exp((path[0][static_cast<std::size_t>(tg[i - 1])] - path[1][static_cast<std::size_t>(tg[i - 1])])*dt);
	}

	long seed = QuantLib::SeedGenerator::instance().get();
	MersenneTwisterUniformRng rnd(seed);
	QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> normrnd(rnd);
	fxRate *= std::exp(fxVol*normrnd.next().value*std::sqrt(defaultTime));

	// calculate default date
	Date defaultDate = disTS->referenceDate();

	// calculate remaining pay and receive dates
	std::function<bool(const Date&)> dateCompare = [&defaultDate](const Date& d) { return d > defaultDate; };
	std::vector<Date>::iterator
		payDateItr = std::find_if(arguments_.payDates.begin(), arguments_.payDates.end(), dateCompare),
		receiveDateItr = std::find_if(arguments_.receiveDates.begin(), arguments_.receiveDates.end(), dateCompare);

	Real payNPV = 0., receiveNPV = 0.;
	if (arguments_.type == CrossCurrencySwap::Type::payDomestic) {
		while (payDateItr != arguments_.payDates.end()) {
			payNPV +=
				arguments_.payCoupons[payDateItr - arguments_.payDates.begin()] * disTS->discount(*payDateItr - defaultDate, true);
			++payDateItr;
		}

		while (receiveDateItr != arguments_.receiveDates.end()) {
			receiveNPV +=
				arguments_.receiveCoupons[receiveDateItr - arguments_.receiveDates.begin()]
				* fxRate*std::exp((domesticRate - foreignRate)*arguments_.receiveDayCount.yearFraction(defaultDate, *receiveDateItr))
				*disTS->discount(*receiveDateItr - defaultDate, true);
			++receiveDateItr;
		}
	} else {
		while (payDateItr != arguments_.payDates.end()) {
			payNPV +=
				arguments_.payCoupons[payDateItr - arguments_.payDates.begin()]
				* fxRate*std::exp((domesticRate - foreignRate)*arguments_.payDayCount.yearFraction(defaultDate, *payDateItr))
				*disTS->discount(*payDateItr - defaultDate, true);
			++payDateItr;
		}

		while (receiveDateItr != arguments_.receiveDates.end()) {
			receiveNPV +=
				arguments_.receiveCoupons[receiveDateItr - arguments_.receiveDates.begin()]
				* disTS->discount(*receiveDateItr - defaultDate, true);
			++receiveDateItr;
		}
	}

	return receiveNPV - payNPV;

}