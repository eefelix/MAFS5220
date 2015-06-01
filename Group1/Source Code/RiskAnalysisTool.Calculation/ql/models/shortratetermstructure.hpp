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

#include <boost/shared_ptr.hpp>

#include <ql/models/shortrate/onefactormodel.hpp>
#include <ql/termstructure.hpp>

//! \file shortratetermstructure.hpp
//! \brief Interest-rate term structure class based on short-rate model
//!

namespace QuantLib {
	//! \brief Interest-rate term structure class based on short-rate model
	/*!	This class is used to generate yield term structure based on the given
		One factor & Affine short-rate model. The discount factors in the term
		structure are all expectations and calculated by the corresponding
		formula of given model.
	*/

	class _RISKANALYSISTOOL_CALCULATION_API ShortRateTermStructure
	{
	public:
		//! \name Constructors & Destructors
		//! @{
		//!	Creates an instance of yield term structure giving the model
		ShortRateTermStructure(
			const Calendar& calendar,
			const DayCounter& daycounter,
			const boost::shared_ptr<OneFactorAffineModel> &model);

		~ShortRateTermStructure();
		//! @}

		//! \name Term structure
		//!
		//! These methods return the yield term structure from reference date
		//! given max date or max time and the short rate at reference date.
		//! @{
		boost::shared_ptr<YieldTermStructure> GetTermStructure(
			Date referenceDate, Date maxDate, Real rate, Frequency freq = Quarterly);

		boost::shared_ptr<YieldTermStructure> GetTermStructure(
			Date referenceDate, Time Tmax, Real rate, Frequency freq = Quarterly);
		//! @}

	private:
		Calendar calendar_;
		DayCounter daycounter_;
		boost::shared_ptr<OneFactorAffineModel> model_;
	};

}