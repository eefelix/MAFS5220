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

#include <vector>

#include <ql/instruments/swap.hpp>
#include <ql/currency.hpp>
#include <ql/types.hpp>
#include <ql/time/date.hpp>
#include <ql/time/schedule.hpp>
#include <ql/time/businessdayconvention.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/cashflow.hpp>

//! \file crosscurrencyswap.hpp
//! \brief fixed-fixed cross currency swap class 
//!
namespace QuantLib
{
	//! \brief fixed-fixed cross currency swap class 
	//!
	//! This class is used to describe fixed-fixed cross currency swap
	//! 
	class _RISKANALYSISTOOL_CALCULATION_API CrossCurrencySwap : public QuantLib::Swap
	{
	public:
		/** An enum type for Swap Type.
		*  The documentation block cannot be put after the enum!
		*/
		enum Type{ 
			payDomestic = 0,	//!< pay
			payForeign	//!< receive
		};

		//! argument class of cross currency swap
		class arguments;
		//! result class of cross currency swap
		class results;
		//! engine class of cross currency swap
		class engine;
		//! \name Constructors & Destructors
		//{@
		CrossCurrencySwap(
			Type type,
			Rate fxRate,
			Volatility fxVol,
			const Date& referenceDate,
			Real payNominal,
			Real receiveNominal,
			const Schedule& paySchedule,
			Rate payRate,
			const DayCounter& payDayCount,
			const Schedule& receiveSchedule,
			Rate receiveRate,
			const DayCounter& receiveDayCount,
			boost::shared_ptr<Currency> payCurrency,
			boost::shared_ptr<Currency> receiveCurrency,
			boost::optional<BusinessDayConvention> paymentConvention = boost::none
			);

		~CrossCurrencySwap();
		//@}
		//!	\name Instrument interface
		//{@
		void setupArguments(PricingEngine::arguments* args) const; 
		void fetchResults(const PricingEngine::results* r) const;
		//@}

		//! \name Inspectors
		//{@
		//! Return type of cross currency swap, which should be payDomestic or payForeign
		Type type() const;
		const Rate fxRate() const;
		const Volatility fxVolatility() const;
		const Date& referenceDate() const;
		Real payNominal() const;
		Real receiveNominal() const;
		const Schedule& paySchedule() const;
		Rate payRate() const;
		const DayCounter& payDayCount() const;
		const Schedule& receiveSchedule() const;
		Rate receiveRate() const;
		const DayCounter& receiveDayCount() const;
		BusinessDayConvention paymentConvention() const;
		const Leg& payLeg() const;
		const Leg& receiveLeg() const;
		//@}
	private:
		Date referenceDate_;

		Type type_;
		Rate fxRate_;
		Volatility fxVol_;
		Real payNominal_;
		Real receiveNominal_;
		Schedule paySchedule_;
		Rate payRate_;
		DayCounter payDayCount_;
		Schedule receiveSchedule_;
		Rate receiveRate_;
		DayCounter receiveDayCount_;
		BusinessDayConvention paymentConvention_;

		boost::shared_ptr < Currency > payCurrency_;
		boost::shared_ptr < Currency > receiveCurrency_;
	};

	//! %Arguments for cross currency swap calculation
	class CrossCurrencySwap::arguments : public Swap::arguments{
	public:
		arguments() : payNominal(Null<Real>()), receiveNominal(Null<Real>())
		{
		}
		Type type;
		Date referenceDate;

		Rate fxRate;
		Volatility fxVol;

		Real payNominal;
		Real receiveNominal;

		DayCounter payDayCount;
		DayCounter receiveDayCount;

		std::vector<Date> payDates;
		std::vector<Date> receiveDates;

		std::vector<Real> payCoupons;
		std::vector<Real> receiveCoupons;

		void validate() const;
	};

	//! %Results for cross currency swap calculation
	class CrossCurrencySwap::results : public Swap::results
	{
	};

	//! base class for cross currency swap pricing engine
	class CrossCurrencySwap::engine : public GenericEngine < CrossCurrencySwap::arguments,
		CrossCurrencySwap::results >
	{
	};

	inline const Date& CrossCurrencySwap::referenceDate()const{ return referenceDate_; }
	inline Real CrossCurrencySwap::payNominal()const{ return payNominal_; }
	inline Real CrossCurrencySwap::receiveNominal()const{ return receiveNominal_; }

	inline const Schedule& CrossCurrencySwap::paySchedule() const{ return paySchedule_; }
	inline Rate CrossCurrencySwap::payRate() const{ return payRate_; }
	inline const DayCounter& CrossCurrencySwap::payDayCount() const{ return payDayCount_; }

	inline const Schedule& CrossCurrencySwap::receiveSchedule() const{ return receiveSchedule_; }
	inline Rate CrossCurrencySwap::receiveRate() const{ return receiveRate_; }
	inline const DayCounter& CrossCurrencySwap::receiveDayCount() const{ return receiveDayCount_; }

	inline BusinessDayConvention CrossCurrencySwap::paymentConvention() const{ return paymentConvention_; }

	inline const Leg& CrossCurrencySwap::payLeg() const{ return legs_[0]; }
	inline const Leg& CrossCurrencySwap::receiveLeg() const{ return legs_[1]; }
}