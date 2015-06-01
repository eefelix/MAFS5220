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
#include <ql/models/defaultmodel.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <ql/processes/cirprocess.hpp>
#include <ql/methods/montecarlo/path.hpp>

//! \file cirdefaultmodel.hpp
//! \brief Intensity default model based on the CIR process
//!

namespace QuantLib {
	//! \brief Intensity default model based on the CIR process
	//!
	//! The model uses CIR process to describe defalut intensity, it should be calibrated by CDS data.
	class _RISKANALYSISTOOL_CALCULATION_API	CirDefaultModel : public DefaultModel, public CoxIngersollRoss
	{
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of CIR default model by setting parameters as default initial values
		CirDefaultModel(
			Rate r0 = 0.05,
			Real theta = 0.1,
			Real k = 0.1,
			Real sigma = 0.1);
		~CirDefaultModel(){};
		//@}
	public:
		//! \name Public interface
		//{@
		//!	Calculate the default probability according to CIR dedault model
		const Probability defaultProbability(Time t) const override;
		//!	Return the default time of the path according to CIR dedault model,
		//! based on Monte Carlo method
		const Time defaultTime(const Path& path) const override;
		//!	Return the process of counterparty after calibration
		boost::shared_ptr<StochasticProcess1D> process() const override;

		//!	Calibrate the CIR default model volatility using sets of CDSs as calibration helper
		void calibrate(
			const std::vector<boost::shared_ptr<CalibrationHelper> >&,
			OptimizationMethod& method,
			const EndCriteria& endCriteria,
			const Constraint& constraint = Constraint(),
			const std::vector<Real>& weights = std::vector<Real>(),
			const std::vector<bool>& fixParameters = std::vector<bool>());
		//@}
	private:
		mutable boost::shared_ptr<CIRprocess> process_;
	};
}
