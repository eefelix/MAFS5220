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
#include <ql/instruments/swap/equityswap.hpp>

//! \file analyticequityswapengine.hpp
//! \brief Analytic engine of equity swap.
//! 

namespace QuantLib{

	//! \brief Analytic engine of equity swap.
	/*! This class is the analytic engine of equity swap.
		The cash expected cash flow from the equity leg 
		and the fixed coupon leg is discounted back to
		reference date.
	*/

	class _RISKANALYSISTOOL_CALCULATION_API AnalyticESEngine : public EquitySwap::engine
	{
	public:
		//! \name Constructors
		AnalyticESEngine(){};

		//! Calculate the npv of the equity swap.
		/*! The swap leg is already constructed 
			the calculation is just discounting 
			the cash flow back to the current time.
		*/
		void calculate() const;
	};
}