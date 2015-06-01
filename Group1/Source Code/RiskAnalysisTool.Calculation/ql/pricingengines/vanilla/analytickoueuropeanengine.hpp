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

#include <boost/shared_ptr.hpp>
#include <ql/instruments/vanillaoption.hpp>
#include <ql/pricingengines/vanilla/mcvanillaengine.hpp>
#include <ql/processes/kouprocess.hpp>

//! \file analytickoueuropeanengine.hpp
//! \brief Analytic engine for european option
//!

namespace QuantLib
{
	//! \brief Analytic engine for european option, whose underlying stock process is a kou process.
	//!
	class _RISKANALYSISTOOL_CALCULATION_API AnalyticKouEuropeanEngine : public VanillaOption::engine
	{
	public:

		//! \name Constructor
		AnalyticKouEuropeanEngine(
			const boost::shared_ptr<KouProcess> &process, 
			const Real &tolerance = 0.0001);
	public:

		//! calculate the european option price
		void calculate() const override;
	private:
		const Real tolerance_;
		boost::shared_ptr<KouProcess> process_;
	};

	//! \brief Kou process helper
	/*! support functions to analytically calculate 
		the option price under kou process assumption
	*/
	class KouHelper{
	public:
		static const double Gamma(
			double mu, double sigma, double lambda,
			double p, double eta1, double eta2, double a,
			double T, double tolerance);
	private:
		static const int combination(int n, int k);
		static const double Hh(int n, double x);
		static const double P(int n, int k, double eta1, double eta2, double p);

		static const double Q(int n, int k, double eta1, double eta2, double p);
		static const double I(int n, double c, double alpha, double beta, double delta);
		static const double KouHelper::Pi(int n, double lambda, double T);

	};
}
