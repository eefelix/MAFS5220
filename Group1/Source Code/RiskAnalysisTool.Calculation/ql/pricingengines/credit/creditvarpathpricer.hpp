/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Xiang, GAO
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
#include <boost/weak_ptr.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/methods/montecarlo/multipath.hpp>

//! \file creditvarpathpricer.hpp
//! \brief Abstract base class used to calculate default value of each path. 
//!
//! This class is used to calculate default value of each path. 
//! To use this class, derived class must implement the pure virtual function defaultNPV.
//! 
//!	\date 2015-04-15

namespace QuantLib {
	//! \brief Abstract base class used to calculate default value of each path. 
	//!
	//! This class is used to calculate default value of each path. 
	//! To use this class, derived class must implement the pure virtual function defaultNPV.
	template <typename ENGINE>
	class CreditVaRPathPricer : public ENGINE::path_pricer_type {
	public:
		typedef typename ENGINE::path_pricer_type::argument_type argument_type;
		typedef typename ENGINE::path_pricer_type::result_type result_type;

		//! \name Constructors 
		//{@
        CreditVaRPathPricer(const boost::shared_ptr<const ENGINE> &engine)
			: engine_(engine) {
		}
		//@}
	public:
		//! \name Destructors
		//{@
        ~CreditVaRPathPricer() {
		}
		//@}
	public:
		//! \name Operator overload
		//{@
		result_type operator() (const argument_type &multiPath) const {
			const boost::shared_ptr<const ENGINE> engine = engine_.lock();
			QL_REQUIRE(engine, "Engine has been destroyed.");

			double cvar = 0.0;

			const boost::shared_ptr<const Counterparty> issuer = engine->counterparties()[0];
			const boost::shared_ptr<const Counterparty> investor = engine->counterparties()[1];
			const std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>>
				processes = engine->processList();

			//QL_REQUIRE(multiPath.assetNumber() >= 2, "missing issuer/investor path information");

			Time issuerDefaultTime = issuer->getDefaultTime(multiPath[0]);
			//Time investorDefaultTime = investor->getDefaultTime(multiPath[1]);

			std::vector<boost::shared_ptr<const MCExposureModel>> tvaModels = engine->models();
			boost::shared_ptr<ShortRateTermStructure> domesticTSModel = engine->domesticShortRateModel();

			if (issuerDefaultTime <= engine->endTime()) {

				std::vector<Path> instrumentPaths;

				for (auto model : tvaModels){
					instrumentPaths.push_back(multiPath[1]); // path of domestic rate
					std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> process = model->instrumentProcess();
					std::vector<std::size_t> pathIndex;
					for (auto p : process){
						auto path = std::find(processes.begin(), processes.end(), p);
						assert(path != processes.end());
						pathIndex.push_back(std::distance(processes.begin(), path));
					}

					for (auto i : pathIndex) {
						instrumentPaths.push_back(multiPath[i]);
					}

					MultiPath instrumentPath(std::move(instrumentPaths));

					Path domesticInterestRatePath = multiPath[1]; // should initiate the path information
					TimeGrid timeGrid = domesticInterestRatePath.timeGrid();

					auto defaultPos = timeGrid.begin();
					while (*defaultPos < issuerDefaultTime)
						++defaultPos;
					Real shortRate = domesticInterestRatePath[timeGrid.index(*defaultPos)];

					// find default date
					Date defaultDate = engine->domesticCurve()->referenceDate() + (int)issuerDefaultTime * 360;
					// construct short rate term structure from default time
					boost::shared_ptr<YieldTermStructure> shortRateDiscountTS = domesticTSModel->GetTermStructure(
						defaultDate, engine->endTime(), shortRate);

					cvar += model->exposure(instrumentPath, issuerDefaultTime, issuerDefaultTime, Handle<YieldTermStructure>(shortRateDiscountTS));

					instrumentPaths.clear();
				}
				cvar = (1 - issuer->getRecoveryRate()) * std::max(cvar, 0.0);
			}


			Array r(1);
			r[0] = cvar;
			return r;
		}
		//@}
	protected:
		const boost::weak_ptr<const ENGINE> engine_;
	};
}