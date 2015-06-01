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
#include <ql/models/defaultmodel.hpp>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ql/termstructures/volatility/equityfx/blackvoltermstructure.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>
#include <ql/methods/montecarlo/path.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/types.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/models/calibrationhelper.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/method.hpp>
#include <ql/time/period.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/date.hpp>
#include <ql/math/array.hpp>
#include <ql/models/model.hpp>

//! \file at1pdefaultmodel.hpp
//! \brief Analytic tractable first passage model
//! 

namespace QuantLib{

	//! \brief Analytic tractable first default model
	/*! The class support the calculations relating to
		the analytical tractable first passage time model
		under the AT1P setting the volatility of the
		the firm value, is piecewise constant.
	*/

	class _RISKANALYSISTOOL_CALCULATION_API AT1Pmodel : public DefaultModel, public CalibratedModel
	{
	public:

		//! \name destructors
		~AT1Pmodel(){};

		//! \name Constructors
		//! @{
		AT1Pmodel(
			boost::shared_ptr<const std::vector<QuantLib::Period>> periods,
			boost::shared_ptr<QuantLib::YieldTermStructure> discountcurve,
			const QuantLib::Real dividend = 0,
			const QuantLib::Real H_V = 0.4,
			const QuantLib::Real B = 0,
			const QuantLib::Real x0 = 1,
			boost::shared_ptr<QuantLib::Array> volatility_ = nullptr);
		//! @}

		//!	Calculate the default probability according to AT1P model
		const QuantLib::Probability defaultProbability(QuantLib::Time T) const override;

		//!	Return the default time of the path according to AT1P model
		const QuantLib::Time defaultTime(const Path &path) const override;

		//!	Return the process of counterparty after calibration
		boost::shared_ptr<StochasticProcess1D> process() const override;

		//!	Calibrate the at1p model volatility using sets of CDSs as calibration helper
		void calibrate(
			const std::vector<boost::shared_ptr<CalibrationHelper> >&,
			OptimizationMethod& method,
			const EndCriteria& endCriteria,
			const Constraint& constraint = Constraint(),
			const std::vector<Real>& weights = std::vector<Real>(),
			const std::vector<bool>& fixParameters = std::vector<bool>()) override;

	private:
		QuantLib::Real defaultBarrier(QuantLib::Time t) const;

		boost::shared_ptr<QuantLib::BlackVolTermStructure> getVolTermSturcture() const;

		QuantLib::Integer betaT(QuantLib::Time T) const;
		QuantLib::Real SigmaT(QuantLib::Time T, QuantLib::Array sigma) const;
		QuantLib::Real SigmaInt(QuantLib::Time T, QuantLib::Array sigma) const;
		QuantLib::Real QtauGreaterThanT(QuantLib::Time T, QuantLib::Array sigma) const;

		const QuantLib::Real H_V;
		const QuantLib::Real B;
		const QuantLib::Real dividend;
		const QuantLib::Real x0;
		Parameter& Volatility;

		boost::shared_ptr<std::vector<QuantLib::Time>> Cdsmaturities;
		boost::shared_ptr<std::vector<Date>> CDSmaturityDates;
		boost::shared_ptr<const std::vector<QuantLib::Period>> Periods;
		boost::shared_ptr<QuantLib::YieldTermStructure> DiscountCurve;
		mutable boost::shared_ptr<QuantLib::BlackVolTermStructure> VolTermStructure;

		mutable boost::shared_ptr<BlackScholesMertonProcess> process_;
	};
}