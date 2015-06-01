/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Tian, XIE

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
#include <ql/instruments/swap/equityswap.hpp>
#include <ql/pricingengines/ucva/mcucvaengine.hpp>
#include <ql/credit/counterparty.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/types.hpp>

//! \file mcequityswapucvaengine.hpp
//! \brief Monte Carlo pricing engine for UCVA calculation of equity swap
//!

namespace QuantLib {
	//! \brief Monte Carlo pricing engine for UCVA calculation of equity swap
	/*!	This class is used to calculate UCVA of equity swap via Monte Carlo simulation,
		giving the information of the Counterparty and the underlying instrument.
		
		The kernel calculation part for UCVA is implemented in its pathpricer class.

		\see MCUCVAEngine
		\see EquitySwapUCVAPathPricer
		\see MCCrossCurrencySwapUCVAEngine
		\see MCVanillaSwapUCVAEngine
	*/
	class _RISKANALYSISTOOL_CALCULATION_API MCEquitySwapUCVAEngine
		: public MCUCVAEngine < EquitySwap >
		, public boost::enable_shared_from_this < MCEquitySwapUCVAEngine > {
	public:

		//! \name Constructor
		MCEquitySwapUCVAEngine(
			Handle<YieldTermStructure> riskFreeTermStructure,
			const boost::shared_ptr<const EquitySwap> &swap,
			const boost::shared_ptr<const Counterparty> &issuer,
			const Matrix& corr,
			Size timeStepsPerYear = 360, bool antitheticVariate = true,
			Size requiredSamples = 50000, Real requiredTolerance = 0.0001, Size maxSamples = QL_MAX_INTEGER,
			BigNatural seed = QuantLib::SeedGenerator::instance().get()
			);

		//! \name Destructor
		~MCEquitySwapUCVAEngine();
	protected:

		virtual boost::shared_ptr<path_pricer_type> pathPricer() const;

		//! return the underlying stock process 
		std::vector<boost::shared_ptr<StochasticProcess1D>> instrumentProcess() const;

	private:
		boost::shared_ptr<Counterparty> C;
		double rho;
		size_t n;
		mutable boost::shared_ptr<QuantLib::BlackScholesMertonProcess> process;
	};
}