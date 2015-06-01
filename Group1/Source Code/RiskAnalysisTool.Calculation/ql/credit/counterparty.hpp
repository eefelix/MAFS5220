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
#include <boost/shared_ptr.hpp>
#include <ql/termstructures/credit/piecewisedefaultcurve.hpp>
#include <ql/math/interpolations/backwardflatinterpolation.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/dategenerationrule.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/models/defaultmodel.hpp>
#include <ql/methods/montecarlo/path.hpp>

//! \file counterparty.hpp
//! \brief Counterparty	class
//!
namespace QuantLib
{
	//! \brief Counterparty	class
	/*!	This class is used to describe defalut behavior of investor and issuer involved in a transaction,
		based on the default model.
	*/
	class _RISKANALYSISTOOL_CALCULATION_API Counterparty
	{
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of Counterparty using its CDS information and the default model.
		Counterparty(
			int settlementDays,
			double recoveryRate,
			const boost::shared_ptr<const std::vector<double>>& cdsSpreads,
			const boost::shared_ptr<const std::vector<Period>>& cdsTenors,
			const boost::shared_ptr<YieldTermStructure>& discountCurve,
			const boost::shared_ptr<DefaultModel>& model_,
			const QuantLib::Frequency& freq = QuantLib::Quarterly,
			const QuantLib::BusinessDayConvention& busDayConvention = QuantLib::Following,
			const QuantLib::DateGeneration::Rule& dateGenerationRule = QuantLib::DateGeneration::TwentiethIMM,
			const QuantLib::Date& referenceDate = QuantLib::Date(),
			const QuantLib::DayCounter& dayCounter = QuantLib::Actual365Fixed(),
			const QuantLib::Calendar& calendar = QuantLib::TARGET(),
			const boost::shared_ptr<YieldTermStructure>& defaultdiscountCurve = 0 /* for fva */);
		~Counterparty();
		//@}

	public:
		//! \name Public interface
		//{@
		//! Return default probability within time t based on the default model.
		const Real getDefaultProb(const Time& t) const;
		//! Return default time given a path describing the evolution of default state variable.
		const Time getDefaultTime(const Path & path) const;
		//! Return random number of default time 
		const Time defaultTimeGenerator(const Probability& prob, const double& epsilon = 1e-2) const;

		//this interface return default time by copula will be implemented later
		//const Time defaultTime() const;

		//! Calibrate default model
		void modelCalibrate(
			OptimizationMethod &method, const EndCriteria &endCriteria,
			bool needCalibrate = true, const Constraint &constraint = Constraint(),
			const std::vector<Real> &weights = std::vector<Real>(),
			const std::vector<bool> &fixParameters = std::vector<bool>()
			);
		//@}

		//! Return the funding spread at time t
		double fundingSpread(Time t) const;

	public:
		//! \name Inspectors
		//{@
		//! Return shared pointer of default model
		const boost::shared_ptr<const DefaultModel> getModel() const;
		const double getRecoveryRate() const { return recoveryRate_; }
		boost::shared_ptr<YieldTermStructure> discountCurve() const { return discountCurve_; }
		const Date referenceDate() const { return referenceDate_; }
		//! Return shared pointer of process describing default state variable
		boost::shared_ptr<StochasticProcess1D> createDefaultProcess() const { return model_->process(); }
		//@}
	private:
		void bootStrap(
			int settlementDays,
			const boost::shared_ptr<const std::vector<double>>& cdsSpreads,
			const QuantLib::Frequency& freq = QuantLib::Quarterly,
			const QuantLib::BusinessDayConvention& busDayConvention = QuantLib::Following,
			const QuantLib::DateGeneration::Rule& dateGenerationRule = QuantLib::DateGeneration::TwentiethIMM,
			const QuantLib::Calendar& calendar = QuantLib::TARGET());
		void cdsHelperGenerator();

	private:
		boost::shared_ptr < QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> >
			hazardRateStructure_;
		const boost::shared_ptr<DefaultModel> model_;
		std::vector<boost::shared_ptr<CalibrationHelper>> cdsHelper_;

		Date referenceDate_;
		double recoveryRate_;
		boost::shared_ptr<YieldTermStructure> discountCurve_;
		boost::shared_ptr<YieldTermStructure> defaultdiscountCurve_;
		boost::shared_ptr<const std::vector<QuantLib::Period>> cdsTenors_;
		DayCounter dayCounter_;

		mutable std::vector<double> fundingSpread_;
		//boost::shared_ptr<const std::vector<double>> cdsSpreads_;
		/*int settlementDays_;
		QuantLib::Frequency freq_;
		QuantLib::BusinessDayConvention busDayConvention_;
		QuantLib::DateGeneration::Rule dateGenerationRule_;
		QuantLib::Calendar calendar_;*/

		bool CALIBRATION_FLAG_;
		mutable bool FUNDINGSPREAD_FLAG_;
	};
}