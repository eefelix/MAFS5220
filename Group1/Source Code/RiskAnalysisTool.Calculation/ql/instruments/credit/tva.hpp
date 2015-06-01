/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Yan, CHEN
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
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/instrument.hpp>
#include <ql/pricingengine.hpp>
#include <ql/credit/counterparty.hpp>
//!	\file tva.hpp
//!	\brief Total value adjustment manager class
namespace QuantLib {
	//! \brief Total value adjustment manager class
	/*!	This class inherits from Instrument class, responsible for calculation total value adjustment (includes
		CVA, DVA, FVA) of a single instrument or portfolio. The arguments are the information of two counterparties
		(issuer and investor) using the Counterparty class and the information of underlying instrument or portfolio
		using their corresponding model class. The Pricingengine uses tvaengine based on Monte Carlo simulation
		framework. The results are including CVA, DVA and FVA.
		\see counterparty.hpp
		\see tvamodel.hpp
		\see tvaengine.hpp
	*/
	class  TVA : public Instrument {
	public:
		//! argument class of TVA manager
		class arguments;
		//! result class of TVA manager
		class results;
		//! engine class of TVA manager
		class engine;
	public:
		virtual ~TVA(){}

	public:
		void setupArguments(PricingEngine::arguments*) const {}
		void fetchResults(const PricingEngine::results*) const;
		bool isExpired() const { return false; }

		const Real CVA() const{
			calculate();
			return CVA_;
		}

		const Real DVA() const{
			calculate();
			return DVA_;
		}

		const Real FVA() const{
			calculate();
			return FVA_;
		}
		
	protected:
		/*std::vector<boost::shared_ptr<const Instrument>> portfolio_;
		std::vector<boost::shared_ptr<const Counterparty>> counterparties_;
		Handle<YieldTermStructure> discountCurve_;*/
		mutable Real CVA_;
		mutable Real DVA_;
		mutable Real FVA_;
	};

	class TVA::arguments : public virtual PricingEngine::arguments{
	public:
		virtual ~arguments() {
		}
		virtual void validate() const {
			/*QL_REQUIRE(portfolio.size() >= 1, "Portfolio should contains at least one instrument.");
			QL_REQUIRE(portfolio.size() == 2, "Exact two counterparties should be specified.");
			QL_REQUIRE(!discountCurve.empty(), "Discount curve is needed");*/
		}
	public:
		/*std::vector<boost::shared_ptr<const Instrument>> portfolio;
		std::vector<boost::shared_ptr<const Counterparty>> counterparties;
		Handle<YieldTermStructure> discountCurve;*/
	};

	class TVA::results : public Instrument::results{
	public:
		virtual ~results() {
		}
		virtual void reset() {
			Instrument::results::reset();
			CVA = Null<Real>();
			DVA = Null<Real>();
			FVA = Null<Real>();
		}
	public:
		Real CVA;
		Real DVA;
		Real FVA;
	};

	class TVA::engine : public GenericEngine < TVA::arguments, TVA::results >{
	};
	
	inline void TVA::fetchResults(const PricingEngine::results* r) const{
		const TVA::results* results =
			dynamic_cast<const TVA::results*>(r);
		QL_ENSURE(results != 0,
			"no results returned from pricing engine");

		CVA_ = results->CVA;
		DVA_ = results->DVA;
		FVA_ = results->FVA;

		additionalResults_ = results->additionalResults;
	}
}
