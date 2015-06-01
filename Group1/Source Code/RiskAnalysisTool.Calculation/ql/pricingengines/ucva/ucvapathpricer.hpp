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
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

//! \file ucvapathpricer.hpp
//! \brief Abstract base class used to calculate cash flow at default time of each path under the setting of UCVA. 
//!
namespace QuantLib {
	//! \brief Abstract base class used to calculate cash flow at default time of each path under the setting of UCVA. 
	/*!	This class is used to calculate cash flow at default time of each path.
		To use this class, derived class must implement the pure virtual function defaultNPV.

		\see MCUCVAEngine
	*/
	template <typename ENGINE>
	class UCVAPathPricer : public ENGINE::path_pricer_type {
	public:
		typedef typename ENGINE::path_pricer_type::argument_type argument_type;
		typedef typename ENGINE::path_pricer_type::result_type result_type;
	protected:
		//! \name Constructors 
		//{@
		UCVAPathPricer(const boost::shared_ptr<const ENGINE> &engine)
			: engine_(engine) {
		}
		//@}
	public:
		//! \name Destructors
		//{@
		~UCVAPathPricer() {
		}
		//@}
	public:
		//! \name Operator overload
		//{@
		//!	Calculate the UCVA given the Multi-path
		result_type operator() (const argument_type &multiPath) const {
			auto engine = engine_.lock();
			QL_REQUIRE(engine, "Engine has been destroyed.");


			auto issuer = engine->issuer();
			Handle<YieldTermStructure> discountCurve = engine->discountCurve();

			std::vector<Path> paths;
			for (int i = 1; i < multiPath.assetNumber(); ++i) {
				paths.push_back(multiPath[i]);
			}
			MultiPath underlyingPath(std::move(paths));

			// get default time for issuer and invertor
			QL_REQUIRE(multiPath.assetNumber() >= 1, "missing issuer/investor path information");
			Time issuerDefaultTime = issuer->getDefaultTime(multiPath[0]);

			// find the firt party default
			Real cva = 0.0;

			if (issuerDefaultTime <= engine->endTime()) {
				cva = (1 - issuer->getRecoveryRate()) * std::max(this->defaultNPV(underlyingPath, issuerDefaultTime), 0.0);
			}

			return cva;
		}
		//@}
	protected:
		//! \name Protected interface
		//{@
		//! Return default NPV at each path, must be implemented in derived class
		virtual Real defaultNPV(const MultiPath& path, Time defaultTime) const = 0;
		
	protected:
		//! Return weak pointer of pricing engine
		const boost::weak_ptr<const ENGINE> engine_;
		//@}

	};
}