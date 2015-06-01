/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2015 Fang Xiong
 Copyright (C) 2015 Xiang Gao

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
#include <ql/processes/squarerootprocess.hpp>
//! \file cirprocess.hpp
//! \brief CIR process class
//!
namespace QuantLib {
	//! \brief CIR process class
	//!
	/*! This class describes a square-root process governed by
		\f[
		dx = a (b - x_t) dt + \sigma \sqrt{x_t} dW_t.
		\f]
		The process is used to model CIR process which provides more precise evolve method.

		\remark The Implict Milstein discretization scheme is recommended to use for efficient simulation.
		Implict Milstein is got from classical Milstein, in which we replace the drift term \f$a(b-r_t)dt\f$
		with \f$a(b-r_{t+dt})dt\f$. We then obtain:
		\f[
		x_{t+dt} = \frac {x_t + ab dt + \sigma \sqrt{x_t}\sqrt{dt} dw + \frac{1}{4} \sigma^2 dt({dw}^2 - 1)} {1 + a dt}
		\f]
	*/
	class _RISKANALYSISTOOL_CALCULATION_API CIRprocess : public StochasticProcess1D {
	public:
		//! An enum type for Discretization method
		enum Discretization {
			Euler,	//!< Euler Discretization Scheme
			Milstein,	//!< Milstein Discretization Scheme
			ImplicitMilstein,	//!< Implicit Milstein Discretization Scheme
			NonCentralChiSquareVariance	//!< Non-Central Chi-Sqaure Distribution
		};
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of CIR process using the four parameters
		//! long-term mean, revert speed, volatility and initial value
		CIRprocess(Real mean, Real speed, Volatility sigma, Real x0 = 0.0, Discretization discretization = ImplicitMilstein);
		~CIRprocess() { }
		//@}

		//! \name Inspectors
		//{@
		Real x0() const;
		Real revertSpeed() const;
		Real revertLevel() const;
		Volatility volatility() const;
		//@}

		//! \name Public interface
		//{@
		//! Return the drift term of the process
		Real drift(Time t, Real x) const;
		//! Return the diffusion of the process
		Real diffusion(Time t, Real x) const;
		//! Return simulation value of the process at time t0+dts
		Real evolve(Time t0, Real x0, Time dt, Real dw) const;
		//@}
	private:
		Real derivDifussion(Time t, Real x) const;
	private:
		Real x0_, revertSpeed_, revertLevel_;
		Volatility volatility_;
		Discretization discretization_;
	};
}