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
#include <ql/models/calibrationhelper.hpp>
#include <ql/instruments/bonds/zerocouponbond.hpp>

//! \file zerocouponbondhelper.hpp
//! \brief Zero coupon bond calibration helper
//!

namespace QuantLib
{
	//! \brief Zero coupon bond calibration helper
	/*! One factor affine short rate model calibration helper using zero coupon bond.
		It can be use to calibrate CIR model and vasicek model, etc.
	*/
	class _RISKANALYSISTOOL_CALCULATION_API ZerocouponbondHelper : public CalibrationHelper {
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of zero coupon bond helper using the information of the bond
		ZerocouponbondHelper(
			Natural settlementDays,
			const Calendar& calendar,
			Real faceAmount,
			const Date& maturityDate,
			BusinessDayConvention paymentConvention,
			Real redemption,
			const Date& issueDate,
			const Handle<Quote>& volatility,
			const Handle<YieldTermStructure>& termStructure,
			CalibrationErrorType calibrationErrorType
			= RelativePriceError);

		//@}
		//! \name Public interface
		//{@
		//! This method is to be used with tree-based methods
		virtual void addTimesTo(std::list<Time>& times) const;
		//!
		//! Calculate model value of zero coupon bond
		virtual Real modelValue() const;
		//!
		//! Calculate black value of zero coupon bond
		virtual Real blackPrice(Volatility volatility) const;
		//@}
	private:
		void performCalculations() const;

		mutable boost::shared_ptr<ZeroCouponBond> zcbond_;

	};
}