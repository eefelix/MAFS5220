/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Yan, CHEN

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
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/instruments/vanillaswap.hpp>
#include <ql/pricingengines/ucva/mcucvaengine.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/processes/cirprocess.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>

//! \file mcvanillaswapucvaengine.hpp
//! \brief Monte Carlo pricing engine for UCVA calculation of vanilla swap 
//!

namespace QuantLib {
	//! \brief Monte Carlo pricing engine for UCVA calculation of vanilla swap
	/*!	This class is used to calculate UCVA of vanilla swap via Monte Carlo simulation,
		giving the information of the Counterparty and the underlying instrument.
		
		The kernel calculation part for UCVA is implemented in its pathpricer class.

		\see MCUCVAEngine
		\see VanillaSwapUCVAPathPricer
		\see MCCrossCurrencySwapUCVAEngine
		\see MCEquitySwapUCVAEngine
	*/
	class _RISKANALYSISTOOL_CALCULATION_API MCVanillaSwapUCVAEngine
		: public MCUCVAEngine < VanillaSwap >
		, public boost::enable_shared_from_this < MCVanillaSwapUCVAEngine > {
	public:
		//! \name Constructors & Destructors
		//@{
		MCVanillaSwapUCVAEngine(
			const Calendar &calender, const DayCounter &dayCounter, Date referenceDate,
			const Handle<YieldTermStructure> &riskFreeTermStructure,
			const boost::shared_ptr<const VanillaSwap> &swap,
			const boost::shared_ptr<const Counterparty> &issuer,
			const boost::shared_ptr<const CoxIngersollRoss> &cirModel,
			const Matrix &corr,
			Size timeStepsPerYear = 360, bool antitheticVariate = true,
			Size requiredSamples = 50000, Real requiredTolerance = 0.0001, Size maxSamples = QL_MAX_INTEGER,
			BigNatural seed = SeedGenerator::instance().get()
			);
		MCVanillaSwapUCVAEngine::~MCVanillaSwapUCVAEngine();
		//@}
		//! \name Public interface
		//@{
		//! Return the cir model 
		const boost::shared_ptr<const CoxIngersollRoss> &Model() const {
			return cirModel_;
		}
		//@}

	protected:
		//! \name Protected interface
		//@{
		//! Return the processes of instruments
		std::vector<boost::shared_ptr<StochasticProcess1D>> instrumentProcess() const {
			std::vector<boost::shared_ptr<StochasticProcess1D>> process;
			Real theta = cirModel_->params()[0];
			Real k = cirModel_->params()[1];
			Real sigma = cirModel_->params()[2];
			Real r0 = cirModel_->params()[3];
			process.push_back(
				boost::shared_ptr<CIRprocess>(new CIRprocess(theta, k, sigma, r0)));
			return process;
		}
		//! Return the shared pointer of path_pricer_type
		virtual boost::shared_ptr<path_pricer_type> pathPricer() const;
		//@}
	private:
		const boost::shared_ptr<const CoxIngersollRoss> cirModel_;
	};
}