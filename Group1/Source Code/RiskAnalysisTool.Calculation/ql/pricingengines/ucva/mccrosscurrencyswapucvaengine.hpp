/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Fang, XIONG

This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - http://quantlib.org/

QuantLib is free software: you can redistribute it and/or modify it
under the terms of the QuantLib license.  You should have received a
copy of the license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<http://quantlib.org/license.shtml>.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#pragma once

#include <Calculation/Calculation.h>
#include <boost/enable_shared_from_this.hpp>
#include <ql/types.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/pricingengines/ucva/mcucvaengine.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/termstructures/voltermstructure.hpp>
#include <ql/instruments/swap/crosscurrencyswap.hpp>
#include <ql/processes/cirprocess.hpp>

//! \file mccrosscurrencyswapucvaengine.hpp
//! \brief Monte Carlo pricing engine for UCVA calculation of fixed-fixed cross currency swap 
//!

namespace QuantLib {
	//! \brief Monte Carlo pricing engine for UCVA calculation of fixed-fixed cross currency swap 
	/*!	This class is used to calculate UCVA of fixed-fixed cross currency swap via Monte Carlo
		simulation, giving the information of the Counterparty and the underlying instrument.
		
		The kernel calculation part for UCVA is implemented in its pathpricer class.

		\see MCUCVAEngine
		\see CrossCurrencySwapUCVAPathPricer
		\see MCVanillaSwapUCVAEngine
		\see MCEquitySwapUCVAEngine
	*/
	class _RISKANALYSISTOOL_CALCULATION_API MCCrossCurrencySwapUCVAEngine
		: public MCUCVAEngine < CrossCurrencySwap >
		, public boost::enable_shared_from_this < MCCrossCurrencySwapUCVAEngine > {
	public:
		//! \name Constructors & Destructors
		//@{
		MCCrossCurrencySwapUCVAEngine(
			const Calendar &calender, const DayCounter &dayCounter, Date referenceDate,
			const Handle<YieldTermStructure> &riskFreeTermStructure,
			const boost::shared_ptr<const CrossCurrencySwap> &swap,
			const boost::shared_ptr<const Counterparty> &issuer,
			Rate domesticRate, Real domesticRateSpeed, Real domesticRateMean, Volatility domesticRateVol,
			Rate foreignRate, Real foreignRateSpeed, Real foreignRateMean, Volatility foreignRateVol, const Matrix &corr,
			Size timeStepsPerYear = 360, bool antitheticVariate = true,
			Size requiredSamples = 50000, Real requiredTolerance = 0.0001, Size maxSamples = QL_MAX_INTEGER,
			BigNatural seed = SeedGenerator::instance().get()
			);
		~MCCrossCurrencySwapUCVAEngine();
		//@}
	protected:
		//! \name Protected interface
		//@{
		//! Override function in base class. Return the processes of interest, which will be used by path generator. 
		virtual std::vector<boost::shared_ptr<StochasticProcess1D>> instrumentProcess() const {
			std::vector<boost::shared_ptr<StochasticProcess1D>> processes(3);
			processes[1] = boost::make_shared<CIRprocess>(domesticRateMean_, domesticRateSpeed_, domesticRateVol_, domesticRate_);
			processes[2] = boost::make_shared<CIRprocess>(foreignRateMean_, foreignRateSpeed_, foreignRateVol_, foreignRate_);

			return processes;
		}

		//! Return shared pointer of path pricer
		virtual boost::shared_ptr<path_pricer_type> pathPricer() const;
		//@}
	private:
		Rate domesticRate_;
		Real domesticRateSpeed_;
		Real domesticRateMean_;
		Volatility domesticRateVol_;
		Rate foreignRate_;
		Real foreignRateSpeed_;
		Real foreignRateMean_;
		Volatility foreignRateVol_;
	};
}