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

#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/processes/geometricbrownianprocess.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/termstructures/credit/probabilitytraits.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>
#include <ql/stochasticprocess.hpp>

//! \file deterministicdefaultmodel.hpp
//! \brief intensity default model based on deterministic default intensity
//!

namespace QuantLib {
	//! \brief Intensity default model based on deterministic default intensity
	//!
	//! The model uses deterministic function to describe defalut intensity, it is boostrapped from CDS data.
	class _RISKANALYSISTOOL_CALCULATION_API DeterministicDefaultModel : public DefaultModel
	{
	public:
		//! \name Constructors & Destructors
		//{@
		DeterministicDefaultModel(){}
		~DeterministicDefaultModel(){}
		//@}
		//! \name Inspectors
		//{@
		void setHazardRateStructure(
			boost::shared_ptr < QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> > newHarzardRateStructure)
		{
			hazardRateStructure_ = newHarzardRateStructure;
		}

	public:
		//!	Calculate the default probability
		virtual const Probability defaultProbability(Time t) const override{ return hazardRateStructure_->defaultProbability(t, true); }
		//!	Return the default time of the path based on Monte Carlo method
		virtual const Time defaultTime(const Path& path) const override;
		//! Return a geometric brownian motion with zero drift and zero volatility
		virtual boost::shared_ptr<StochasticProcess1D> process() const override{ return boost::make_shared<GeometricBrownianMotionProcess>(0, 0, 0); }
		//@}
	private:
		boost::shared_ptr < QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> >
			hazardRateStructure_;
	};
}
