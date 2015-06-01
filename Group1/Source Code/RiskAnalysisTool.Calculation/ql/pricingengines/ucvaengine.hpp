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
#include <ql/pricingengine.hpp>

//! \file ucvaengine.hpp
//! \brief Abstract base class for ucva engine. 
//!
namespace QuantLib {
	//! \brief Abstract base class for ucva engine.
	/*!	This class is base class for ucva engine, which just contains information as least as possible.
		Each derived class should implement its own method to calculate ucva.
	*/
    template <typename INST>
	class UCVAEngine : public INST::engine {
	public:
		//! \name Constructors & Destructors
		//{@
		//! Create an instance of UCVA Engine using the information of counterparty and underlying instrument.
		UCVAEngine(
			const boost::shared_ptr<const INST> instrument, const boost::shared_ptr<const Counterparty> issuer,
			const Time &endTime, const Handle<YieldTermStructure> &discountCurve
			) : instrument_(instrument), issuer_(issuer), discountCurve_(discountCurve), endTime_(endTime) {
		}

		virtual ~UCVAEngine() {
		}
		//@}
	public:
		//! \name Inspector
		//{@
		const Handle<YieldTermStructure> & discountCurve() const {
			return discountCurve_;
		}

		Time endTime() const {
			return endTime_;
		}

		const boost::shared_ptr<const Counterparty> & issuer() const {
			return issuer_;
		}

		const boost::shared_ptr<const INST> & instrument() const {
			return instrument_;
		}
		//@}
    protected:
        const Time endTime_;
        const boost::shared_ptr<const INST> instrument_;
        const boost::shared_ptr<const Counterparty> issuer_;
        const Handle<YieldTermStructure> discountCurve_;
    };
}