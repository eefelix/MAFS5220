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
#include <ql/instruments/swap.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/schedule.hpp>
#include <ql/pricingengine.hpp>
#include <boost/shared_ptr.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/time/calendars/target.hpp>

//! \file equityswap.hpp
//! \brief equity swap class
//! 
namespace QuantLib {
	//! \brief equity swap class
	/*! This class defines the instrument equity swap
		Equity swap is a contract between two parties,
		which one party pays fixed rate coupon on the notional,
		the other party pays dividend return from the stock.
		At maturity, the notional is exchanged with the stock.
		For simplicity, the dividend yield and the volatility
		of the stock is supposed constant.

		See Counterparty Credit Risk, Collateral and Funding(Brigo 2013)
		page 169 for more detailed descriptions.

		The two legs of the equity swap is initialized using
		the private "initialize" method, whenever an equity swap
		is constructed.
	*/

	class _RISKANALYSISTOOL_CALCULATION_API EquitySwap : public QuantLib::Swap {
	public:
		/** An enum type.
		*  The documentation block cannot be put after the enum!
		*/
		enum Type 
		{
			Receiver = -1, //!< receive fixed coupon
			Payer = 1      //!< pay fixed coupon
		};

		//! argument class of equity swap
		class arguments;

		//! result class of equity swap
		class results;

		//! engine class of equity swap
		class engine;

		//! \name Destructor
		~EquitySwap() {};

		//! \name Constructors
		//! @{
		EquitySwap(
			Type type,
			QuantLib::Date startdate,
			QuantLib::Date referencedate,
			QuantLib::Real startprice,
			QuantLib::Real spotprice,
			QuantLib::Integer amount,
			QuantLib::Real cumulativecoupoon, /* per share */
			QuantLib::Real cumulativedividend,/* per share */
			QuantLib::Rate fixedrate, // fixed rate in swap for the stock dividend
			QuantLib::Rate dividend, // constant dividend yield
			QuantLib::Time maturity, /* in years */
			QuantLib::Real sigma, /* the volatitlity of the stock */
			boost::shared_ptr<QuantLib::YieldTermStructure> discountcurve,
			const QuantLib::Frequency& freq = QuantLib::Semiannual,
			const QuantLib::BusinessDayConvention& busDayConvention = QuantLib::Following,
			const QuantLib::DateGeneration::Rule& rule_ = QuantLib::DateGeneration::Backward,
			const QuantLib::Calendar& calendar = TARGET());

		EquitySwap(const EquitySwap&);
		//! @}

		//! return the number of payoffs 
		size_t Count() const;

		//! return the maturity
		Real getMaturity() const;

		//! return the continuous dividend yield of the stock
		Real getDividendYield() const;

		//! return the current spot price of the stock
		Real getSpotPrice() const;

		//! return the constant volatility of the stock
		Real getVolatility() const;

		//! return the calendar of used in the calculation
		Calendar getCalendar() const;

		bool isExpired() const
		{ 
			if(referenceDate_>=startDate_)
				return false;
			else return true;
		}


	private:
		size_t count_;
		Type type_;
		QuantLib::Date referenceDate_;
		QuantLib::Date startDate_;
		QuantLib::Real startPrice_;
		QuantLib::Real spotPrice_;
		QuantLib::Integer amount_;
		QuantLib::Time maturity_;
		QuantLib::Rate fixedRate_;
		QuantLib::Rate dividend_;
		QuantLib::Real cumulativeCoupoon_;
		QuantLib::Real cumulativeDividend_;
		boost::shared_ptr<QuantLib::YieldTermStructure> discountCurve_;
		QuantLib::Frequency frequency_;
		QuantLib::BusinessDayConvention BusDayConvention_;
		QuantLib::DateGeneration::Rule rule_;
		QuantLib::Calendar calendar_;
		QuantLib::Schedule schedule_;
		std::vector<QuantLib::Time> yearfraction_; /* yearfraction(0) = 0 */
		QuantLib::Volatility sigma_;

		void initialize();
		void setupArguments(QuantLib::PricingEngine::arguments*) const;
		void fetchResults(const QuantLib::PricingEngine::results*) const;
	};

	//! %Arguments for equity swap calculation
	class EquitySwap::arguments : public QuantLib::Swap::arguments{
	public:
		arguments() {
		}
		void validate() const {
			QL_REQUIRE(discountCurve, "Discount Curve Required!!!");
			QL_REQUIRE(type, "Side Required!!!");
			QL_REQUIRE(spotPrice, "Recovery rate Required!!!");
			QL_REQUIRE(startPrice, "Recovery rate Required!!!");
			QL_REQUIRE(amount, "Recovery rate Required!!!");
			QL_REQUIRE(fixedRate, "Interest Rate Required!!!");
			QL_REQUIRE(dividend, "Dividend Required!!!");
			QL_REQUIRE(maturity, "Recovery rate Required!!!");
			QL_REQUIRE(sigma, "Recovery rate Required!!!");
		}
		boost::shared_ptr<QuantLib::YieldTermStructure> discountCurve;
		Type type;
		QuantLib::Date referenceDate;
		QuantLib::Date startDate;
		QuantLib::Real spotPrice;
		QuantLib::Real startPrice;
		QuantLib::Integer amount;
		QuantLib::Rate fixedRate;
		QuantLib::Rate dividend;
		QuantLib::Time maturity;
		QuantLib::Real sigma;
		QuantLib::Real cumulativeCoupoon;
		QuantLib::Real cumulativeDividend;
		QuantLib::Frequency freq;
		boost::shared_ptr<QuantLib::Schedule> schedule;
		size_t count;
		std::vector<QuantLib::Time> yearfraction;
	};

	//! %Results class of equity swap
	class EquitySwap::results : public Swap::results{
	};

	//! equity option engine base class
	class EquitySwap::engine :
		public QuantLib::GenericEngine < EquitySwap::arguments,
		EquitySwap::results >
	{
	};

	// inline functions

	inline Real EquitySwap::getDividendYield() const{
		return dividend_;
	}

	inline size_t EquitySwap::Count() const{
		return count_;
	}

	inline Real EquitySwap::getMaturity() const
	{
		return maturity_;
	}

	inline Real EquitySwap::getSpotPrice() const
	{
		return spotPrice_;
	}

	inline Real EquitySwap::getVolatility() const
	{
		return sigma_;
	}

	inline Calendar EquitySwap::getCalendar() const
	{
		return calendar_;
	}
}
