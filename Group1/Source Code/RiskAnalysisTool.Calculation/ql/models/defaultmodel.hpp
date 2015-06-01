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
#include <ql/models/model.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/methods/montecarlo/path.hpp>
#include <ql/types.hpp>
#include <ql/termstructures/defaulttermstructure.hpp>

//! \file defaultmodel.hpp
//! \brief Abstract base class for default model
//!

namespace QuantLib {
	//! \brief abstract base class for default model
	/*!	Abstract base class for default model, define pure virtual functions need
		to be implemented in derived class.
	*/
	class _RISKANALYSISTOOL_CALCULATION_API DefaultModel /*: public CalibratedModel*/ {
	public:
		//! \name Constructors & Destructors
		//{@
		DefaultModel(){}
		virtual ~DefaultModel(){}
		//@}

		//! \name Publice interface
		//@{
		//!	Calculate the default probability at time t according to the given dedault model
		virtual const Probability defaultProbability(Time t) const = 0;

		//!	Return the default time of the path according to the given dedault model,
		//! based on Monte Carlo method
		virtual const Time defaultTime(const Path& path) const = 0;

		//! Return stochastic process describing default state variable
		virtual boost::shared_ptr<StochasticProcess1D> process() const = 0;
		//@}
	};
}
