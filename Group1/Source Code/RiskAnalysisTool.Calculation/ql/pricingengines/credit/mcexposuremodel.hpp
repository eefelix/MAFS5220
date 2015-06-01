/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Yan, CHEN

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
#include <ql/instrument.hpp>
#include <ql/pricingengine.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/methods/montecarlo/multipath.hpp>
#include <ql/time/all.hpp>
//! \file mcexposuremodel.hpp
//! \brief Base class for exposure calculation model.
//!
namespace QuantLib {
	//! \brief base class for exposure calculation model.
	/*!	This class is a base class for exposure calculation model, working as tva calculation helper, responsible
		for providing information of underlying instrument and calculation the exposure of cash flow at given time.
	*/
	class MCExposureModel {
	public:
		//! \name Constructors & Destructors
		//{@
		MCExposureModel(
			const Calendar &calendar, const DayCounter &dayCounter, Date referenceDate,
			const std::vector<boost::shared_ptr<StochasticProcess1D>>& instrumentProcess,
			const boost::shared_ptr<const Instrument> &instrument)
			: calendar_(calendar), dayCounter_(dayCounter), referenceDate_(referenceDate),
			instrumentProcess_(instrumentProcess), instrument_(instrument) {
		}

		virtual ~MCExposureModel() {
		}
		//@}

	public:
		//! \name Public interface
		//{@
		//! Return shared pointer to underlying instrument
		boost::shared_ptr<const Instrument> instrument() const {
			return instrument_;
		}
		//! Return shared pointer to stochastic procss array containing underlying process
		const std::vector<boost::shared_ptr<StochasticProcess1D>>& instrumentProcess() const{
			return instrumentProcess_;
		};

		//! Calculate the exposure of cash flow at given time and multi-paths information
		virtual Real exposure(const MultiPath& path, const Time issuerDefaultTime,
			const Time investorDefaultTime, const Handle<YieldTermStructure> disTS) const = 0;
		//@}
	protected:
		Calendar calendar_;
		DayCounter dayCounter_;
		Date referenceDate_;

		const boost::shared_ptr<const Instrument> instrument_;
		const std::vector<boost::shared_ptr<StochasticProcess1D>> instrumentProcess_;
	};
}