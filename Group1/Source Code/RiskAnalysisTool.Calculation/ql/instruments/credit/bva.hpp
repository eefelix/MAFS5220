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
#include <ql/types.hpp>
#include <ql/math/matrix.hpp>
#include <ql/instrument.hpp>
#include <ql/pricingengine.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/credit/counterparty.hpp>
//! \file bva.hpp
//! \brief Bilateral value adjustment manager class
//!
namespace QuantLib {
	//! \brief Bilateral value adjustment manager class
	/*!	This class inherits from Instrument class, responsible for calculation bilateral value adjustment (includes
	CVA, DVA) of a single instrument or portfolio. The arguments are the information of two counterparties
	(issuer and investor) using the Counterparty class and the information of underlying instrument or portfolio
	using their corresponding model class. The Pricingengine uses tvaengine based on Monte Carlo simulation
	framework. The results are including CVA, DVA.
	\see counterparty.hpp
	\see mctvamodel.hpp
	\see mcbvaengine.hpp
	*/
	class BVA : public Instrument {
	public:
		//! argument class of BVA manager
		class arguments;
		//! result class of BVA manager
		class results;
		//! engine class of BVA manager
		class engine;
	public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of BVA manager using two counterparties, instruments (single instrument or portfolio),
		//!	correlation matrix of the counterparties and instruments, and the yield term structure.
		BVA(
			/*const boost::shared_ptr<Counterparty> &self,
			const boost::shared_ptr<Counterparty> &counterparty,
			const std::vector<boost::shared_ptr<const Instrument>> &portfolio,
			const Matrix &correlation,
			const Handle<YieldTermStructure> &discountCurve*/
			)
			/*: self_(self)
			, counterparty_(counterparty)
			, portfolio_(portfolio)
			, correlation_(correlation)
			, discountCurve_(discountCurve)*/ {
		}
		//@}
	public:
		//! \name Inspectors
		//{@
		/*boost::shared_ptr<Counterparty> self() {
			return self_;
		}

		boost::shared_ptr<Counterparty> counterparty() {
			return counterparty_;
		}

		Handle<YieldTermStructure> discountCurve() {
			return discountCurve_;
		}*/
		//@}

		//! \name Public interface
		//{@
		//! Calculate the CVA and return the result, or return the result directly if it has calculated.
		const Real CVA() const{
			calculate();
			return CVA_;
		}
		//! Calculate the DVA and return the result, or return the result directly if it has calculated.
		const Real DVA() const{
			calculate();
			return DVA_;
		}

		void setupArguments(PricingEngine::arguments* args) const override{};
		void fetchResults(const PricingEngine::results*) const;
		bool isExpired() const { return false; }
		//@}
	protected:
		/*std::vector<boost::shared_ptr<const Instrument>> portfolio_;
		const boost::shared_ptr<Counterparty> self_;
		const boost::shared_ptr<Counterparty> counterparty_;
		Matrix correlation_;
		Handle<YieldTermStructure> discountCurve_;*/

		mutable Real CVA_;
		mutable Real DVA_;
	};

	class BVA::arguments : public virtual PricingEngine::arguments {
	public:
		virtual ~arguments() {
		}
		void validate() const {
		}
	};

	class BVA::results : public virtual Instrument::results {
	public:
		Real CVA;
		Real DVA;
	public:
		void reset() {
			Instrument::results::reset();
			CVA = Null<Real>();
			DVA = Null<Real>();
		}
	};

	class BVA::engine
		: public GenericEngine < BVA::arguments, BVA::results > {
	};

	inline void BVA::fetchResults(const PricingEngine::results* r) const{
		const BVA::results* results =
			dynamic_cast<const BVA::results*>(r);
		QL_ENSURE(results != 0,
			"no results returned from pricing engine");

		CVA_ = results->CVA;
		DVA_ = results->DVA;

		additionalResults_ = results->additionalResults;
	}
}