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

#include <Calculation\Calculation.h>
#include <boost/shared_ptr.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/models/calibrationhelper.hpp>
#include <ql/models/model.hpp>
//! \file defaultcdshelper.hpp
//! \brief CDS helpers to calibrate default model
//!
namespace QuantLib
{
	//! \brief CDS helpers to calibrate default model
	//!
	//! Calculate default probability based on default and market data, it will be used 
	//! to calibrate default model.
	class _RISKANALYSISTOOL_CALCULATION_API DefaultCdsHelper : public CalibrationHelper{
	public:
		//! \name Constructors & Destructors
		//{@
		DefaultCdsHelper(
			Time cdsTime,
			const Handle<Quote>& impliedDefaultProb,
			const Handle<YieldTermStructure>& termStructure,
			const boost::shared_ptr<CalibratedModel> & model
			);
		~DefaultCdsHelper();
		//@}

		//! \name Public interface
		//{@
		//! Calculate default probability based on default model
		Real modelValue() const override;
		//! Return default probability boostrapped from CDS data
		Real blackPrice(Volatility volatility) const override;
		//@}

	private:
		Time cdsTime_;
		boost::shared_ptr<CalibratedModel> model_;
		void addTimesTo(std::list<Time>& times) const override;
	};
}