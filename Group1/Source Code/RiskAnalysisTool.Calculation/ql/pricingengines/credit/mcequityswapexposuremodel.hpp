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
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ql/instruments/swap/equityswap.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/pricingengines/credit/mcexposuremodel.hpp>
#include <ql/models/shortrate/onefactormodel.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/methods/montecarlo/multipath.hpp>
//! \file mcequityswapexposuremodel.hpp
//! \brief exposure calculation model for equity swap.
//!
namespace QuantLib {
	//! \brief TVA model for equity swap.
	/*! This class implements tva model for equity swap, working as tva calculation helper, responsible for
		providing information of underlying instrument and calculation the exposure of cash flow at given time.
		It's used in the calculation of CVA, DVA and FVA of underlring instrument under Monte Carlo Simulation
		framework.

		\remark The exposure of equity swap at given time t, NPV(t), is calculated as the discounted cash flow.
		The payoff of floating leg at time \f$t_i\f$ is
		\f[
		\frac {q(E(S_{t_i})-E(S_{t_{i-1}}))} {r-q}
		\f]
		where \f$E(S_{t_i}) = S_0e^{(r-q)t_i}\f$.

		\see MCVanillaSwapExposureModel
		\see MCCrossCurrencySwapExposureModel
		*/
	class _RISKANALYSISTOOL_CALCULATION_API MCEquitySwapExposureModel : public MCExposureModel
	{
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of equity swap model using the information of the contract including the
		//! calendar, day counter, reference date, yield term structure, and the instrument.
		MCEquitySwapExposureModel(const Calendar &calendar, const DayCounter &dayCounter, Date referenceDate,
			/*const boost::shared_ptr<const OneFactorAffineModel> &shortRateModel,*/
			const boost::shared_ptr<StochasticProcess1D>& stockProcess,
			const boost::shared_ptr<const EquitySwap> &instrument)
			:MCExposureModel(calendar, dayCounter, referenceDate
			, std::vector<boost::shared_ptr<StochasticProcess1D>>({ stockProcess }),instrument)
			/*shortRateModel_(shortRateModel)*/ {
				arguments_ = boost::make_shared<EquitySwap::arguments>();
				instrument_->setupArguments(arguments_.get());
		};
		//@}
		//! \name Public interface
		//{@
		//! Return the exposure of cash flow at given time and paths information.
		Real exposure(const MultiPath& path, const Time issuerDefaultTime, const Time investorDefaultTime,
			const Handle<YieldTermStructure> disTS) const;
		//@}
	private:
		/*const boost::shared_ptr<const OneFactorAffineModel> shortRateModel_;*/
		boost::shared_ptr<EquitySwap::arguments> arguments_;
	};
}