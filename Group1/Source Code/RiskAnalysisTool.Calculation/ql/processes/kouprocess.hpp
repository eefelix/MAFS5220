/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Tian XIE

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
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/processes/eulerdiscretization.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
//! \file kouprocess.hpp
//! \brief Kou process class
//!
namespace QuantLib
{
	//!	\brief Kou process class
	/*!	This class describes a double exponential jump process, initiated by Steven KOU,
		which is a compromise between reality and tractability. It gives an explanation
		of the two empirical phenomena which received much attention in financial markets:
		the asymmetric leptokurtic feature and the volatility smile. It permits to obtain 
		analytical solutions to the prices of many derivatives : European call and put options;
		interest rate derivatives, such as swaptions, caps, floors, and bond options;
		as well as path-dependant options, such as perpetual American options, barrier,
		and lookback options .
	*/
	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type = __URng_Poisson_Type>
	class /*_RISKANALYSISTOOL_CALCULATION_API*/ GeneralizedKouProcess : public BlackScholesMertonProcess
	{
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of Kou process using the four parameters: initial value, 
		//! risk-free term structure, dividend term structure, black volatility term structure
		//! and the additional jump information
		GeneralizedKouProcess(
			const Handle<Quote>& x0,
			const Handle<YieldTermStructure>& dividendTS,
			const Handle<YieldTermStructure>& riskFreeTS,
			const Handle<BlackVolTermStructure>& blackVolTS,
			const Real jumpIntensity,
			const Real posProbability,
			const Real posJumpMean,
			const Real negJumpMean,
			const __URng_Poisson_Type& URng_Poisson = __URng_Poisson_Type(SeedGenerator::instance().get()),
			const __URng_DoubleExpDist__Type& URng_DoubleExpDist = __URng_Poisson_Type(SeedGenerator::instance().get()),
			const boost::shared_ptr<discretization>& d = boost::shared_ptr<discretization>(new EulerDiscretization)
			);
		//@}

		//! \name Public interface
		//{@
		//! Return simulation value of the process at time t0+dt
		Real evolve(Time t0, Real x0, Time dt, Real dw) const override;
		const Real jumpIntensity() const;
		const Real posProbability() const;
		const Real posJumpMean() const;
		const Real negJumpMean() const;
		//@}
	private:
		// Additional jump pameters
		Real jumpIntensity_;
		Real posProbability_;
		Real posJumpMean_;
		Real negJumpMean_;
		// Uniform random numbers generators
		__URng_DoubleExpDist__Type URng_Poisson_;
		__URng_DoubleExpDist__Type URng_DoubleExpDist_;

		//	adjustment for jump, called by evolve
		Real drift_adj(Time dt) const;
		Real diffusion_adj(Time dt) const;
	};

	typedef GeneralizedKouProcess<MersenneTwisterUniformRng, MersenneTwisterUniformRng> KouProcess;

	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::GeneralizedKouProcess(
		const Handle<Quote>& x0,
		const Handle<YieldTermStructure>& dividendTS,
		const Handle<YieldTermStructure>& riskFreeTS,
		const Handle<BlackVolTermStructure>& blackVolTS,
		const Real jumpIntensity,
		const Real posProbability,
		const Real posJumpMean,
		const Real negJumpMean,
		const __URng_Poisson_Type& URng_Poisson,
		const __URng_DoubleExpDist__Type& URng_DoubleExpDist,
		const boost::shared_ptr<discretization>& d
		) :	
		BlackScholesMertonProcess(x0, dividendTS, riskFreeTS, blackVolTS, d),
		jumpIntensity_(jumpIntensity), posProbability_(posProbability),
		posJumpMean_(posJumpMean), negJumpMean_(negJumpMean),
		URng_Poisson_(URng_Poisson),
		URng_DoubleExpDist_(URng_DoubleExpDist)
	{
	}

	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::drift_adj(Time dt) const
	{
		//	adjustment for drift \lambda*E[J]
		/*
		return jumpIntensity_*dt*(posProbability_*posJumpMean_*(1 - posJumpMean_)
			- (1 - posProbability_)*negJumpMean_*(1 - negJumpMean_));
		*/
		throw std::exception("Not implemented.");
		return 0.0;
	}

	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::diffusion_adj(Time dt) const
	{
		//	adjustment for diffusion
		InverseCumulativePoisson invP(jumpIntensity_*dt);

		//	poisson random number, number of jumps
		int Pt = (int)invP(URng_Poisson_.next().value);
		Real DoubleExp;
		Real totalJump = 0.;

		for (int i = 1; i <= Pt; ++i) {
			DoubleExp = URng_DoubleExpDist_.next().value;

			// inverse uniform random number to get double exponential random number
			if (DoubleExp <= 1 - posProbability_)
				totalJump += (1 - posProbability_)
				*(exp(posJumpMean_*log(DoubleExp / (1 - posProbability_))) - 1);
			else
				totalJump += posProbability_*
				(exp(-negJumpMean_*log((1 - DoubleExp) / posProbability_)) - 1);
		}

		return totalJump;
	}

	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	inline Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::evolve(Time t0, Real x0, Time dt, Real dw) const
	{
		return apply(x0, discretization_->drift(*this, t0, x0, dt) /*- drift_adj(dt)*/
			+ stdDeviation(t0, x0, dt)*dw + diffusion_adj(dt));
	}

	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	inline const Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::jumpIntensity() const
	{
		return jumpIntensity_;
	}
	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	inline const Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::posProbability() const
	{
		return posProbability_;
	}
	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	inline const Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::posJumpMean() const
	{
		return posJumpMean_;
	}
	template <typename __URng_Poisson_Type, typename __URng_DoubleExpDist__Type>
	inline const Real GeneralizedKouProcess<__URng_Poisson_Type, __URng_DoubleExpDist__Type>::negJumpMean() const
	{
		return negJumpMean_;
	}
}