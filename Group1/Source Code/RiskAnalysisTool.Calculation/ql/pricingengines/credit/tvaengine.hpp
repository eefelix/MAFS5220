/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Xiang, GAO

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
#include <ql/credit/counterparty.hpp>
#include <ql/models/shortratetermstructure.hpp>
//! \file tvaengine.hpp
//! \brief Base class for tva engine.

namespace QuantLib {
	//! \brief Base class for bva engine. 
	/*!	This class is base class for tva engine, which is used to set some basic information.
	Each derived class should implement its own method to calculate bva.
	*/
	template<typename INST>
	class TVAEngine : public INST::engine {
	public:
		//! \name Constructors & Destructors
		//{@
		//! Create an instance of TVA Engine using the information of the two counterparties using
		//!	Counterparty class and underlying instruments using instrument model class.
		TVAEngine(
			Time endTime,
			const std::vector<boost::shared_ptr<const Counterparty>> &counterparties,
			const std::vector<boost::shared_ptr<const MCExposureModel>> &tvaModels,
			const Handle<YieldTermStructure> &domseticCurve, const boost::shared_ptr<OneFactorAffineModel>& domesticCIRModel
			) : endTime_(endTime), counterparties_(counterparties), tvaModels_(tvaModels), domesticCurve_(domseticCurve), domesticCIRModel_(domesticCIRModel)
		{
			domseticShortRateModel_ = boost::make_shared<ShortRateTermStructure>(domesticCurve_->calendar(), domesticCurve_->dayCounter(), domesticCIRModel_);
		}

		virtual ~TVAEngine() {
		}
		//@}
	public:
		//! \name Inspector
		//{@
		std::vector<boost::shared_ptr<const Counterparty>> counterparties() const {
			return counterparties_;
		}

		/*std::vector<boost::shared_ptr<const Instrument>> portfolio() {
			return portfolio_;
			}*/

		const Time endTime() const { return endTime_; }
		const Handle<YieldTermStructure> & domesticCurve() const {
			return domesticCurve_;
		}

		const boost::shared_ptr<ShortRateTermStructure>& domesticShortRateModel() const{
			return domseticShortRateModel_;
		}

		const boost::shared_ptr<OneFactorAffineModel>& domesticCIRModel() const{
			return domesticCIRModel_;
		}
		const std::vector<boost::shared_ptr<const MCExposureModel>>& models() const{
			return tvaModels_;
		}
		//@}
	protected:
		//std::vector<boost::shared_ptr<const Instrument>> portfolio_;
		std::vector<boost::shared_ptr<const Counterparty>> counterparties_;
		Handle<YieldTermStructure> domesticCurve_;
		boost::shared_ptr<OneFactorAffineModel> domesticCIRModel_;
		std::vector<boost::shared_ptr<const MCExposureModel>> tvaModels_;
		boost::shared_ptr<ShortRateTermStructure> domseticShortRateModel_;
		Time endTime_;
	};
}