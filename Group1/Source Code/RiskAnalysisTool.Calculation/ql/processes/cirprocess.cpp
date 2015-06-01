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

#include <pch.h>
#include "cirprocess.hpp"
#include <ql/math/distributions/chisquaredistribution.hpp>
#include <ql/math/distributions/normaldistribution.hpp>

using namespace QuantLib;

CIRprocess::CIRprocess(
	Real mean, Real speed, Volatility sigma,
	Real x0 /*= 0.0*/,
	Discretization discretization)
	:StochasticProcess1D(),
	x0_(x0), revertLevel_(mean), revertSpeed_(speed), volatility_(sigma)
	, discretization_(discretization) {
	//QL_REQUIRE(2 * revertLevel_*revertSpeed_ >= volatility_*volatility_,
	//	"CIR stable conditon not be satisfied, zero interest rate will be produced.");
}

Real CIRprocess::x0() const { return x0_; }
Real CIRprocess::revertSpeed() const { return revertSpeed_; }
Real CIRprocess::revertLevel() const { return revertLevel_; }
Volatility CIRprocess::volatility() const { return volatility_; }

//Real CIRprocess::drift(Time t, Real x) const
//{
//	throw std::exception("not implemented.");
//	return -1.;
//}
//
//Real CIRprocess::diffusion(Time t, Real x) const
//{
//	throw std::exception("not implemented.");
//	return -1.;
//}

Real CIRprocess::drift(Time, Real x) const {
	return revertSpeed_*(revertLevel_ - x);
}

Real CIRprocess::diffusion(Time, Real x) const {
	return volatility_*std::sqrt(std::max(x, 0.0));
}

Real CIRprocess::derivDifussion(Time, Real x) const {
	return 0.5 * volatility_ / sqrt(std::max(x, QL_EPSILON));
}

Real CIRprocess::evolve(Time t0, Real x0, Time dt, Real dw) const {
	Real xdt;
	switch (discretization_) {
	case Euler:
		xdt = x0 + drift(t0, x0) * dt + diffusion(t0, x0) * dw;
		break;
	case Milstein:
		xdt = x0 + drift(t0, x0) * dt + diffusion(t0, x0) * dw
			+ 0.5 * diffusion(t0, x0) * derivDifussion(t0, x0) * (dw * dw - dt);
		break;
	case ImplicitMilstein:
		xdt = (x0 + revertLevel_*revertSpeed_*dt + volatility_*std::sqrt(x0)*dw*std::sqrt(dt) +
			0.25*volatility_*volatility_*dt*(dw*dw - 1)) / (1 + revertSpeed_*dt);
		break;
	case NonCentralChiSquareVariance:
	{
		double coef = (1 - std::exp(-revertSpeed_*dt))*volatility_*volatility_ / 4. / revertSpeed_;
		double df = 4 * revertLevel_*revertSpeed_ / volatility_ / volatility_;
		double ncp = x0*std::exp(-revertSpeed_*dt) / coef;

		QuantLib::InverseNonCentralChiSquareDistribution NCCS(df, ncp, 100);
		xdt = coef*NCCS(NormalDistribution(0, 1)(dw));
		break;
	}
	default:
		QL_FAIL("Invalid discretization method.");
	}

	if (xdt <= 0.0) {
		//QL_DEBUG("CIRprocess::evolve producing a non-positive value, returning EPSILON instead.");
		xdt = QL_EPSILON;
	}

	return xdt;
}