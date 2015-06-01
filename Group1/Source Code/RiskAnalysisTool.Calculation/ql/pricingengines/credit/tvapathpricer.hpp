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
#include <boost/weak_ptr.hpp>
#include <ql/credit/counterparty.hpp>

//! \file tvapathpricer.hpp
//! \brief Path pricer class used to calculate the TVA on each Multi-path for TVAEngine class. 
//!
namespace QuantLib {
	//! \brief Path pricer class used to calculate the TVA on each Multi-path for TVAEngine class.
	/*!	This class is the core part for TVAEngine class, which is used to calculate the TVA (includes CVA, DVA and FVA)
		on each Multi-path given the instruments (single instrument or portfolio), based on the definition of TVA
		(includes CVA, DVA and FVA) in Brigo's book: "Counterparty Credit Risk, Collateral and Funding".

		If the underlying is a portfolio, netting is considered.

		\see MCTVAEngine
		*/
	template <typename ENGINE>
	class TVAPathPricer : public ENGINE::path_pricer_type {
	public:
		typedef typename ENGINE::path_pricer_type::argument_type argument_type;
		typedef typename ENGINE::path_pricer_type::result_type result_type;

		//! \name Constructors 
		//{@
		TVAPathPricer(const boost::shared_ptr<const ENGINE> &engine)
			: engine_(engine) {
		}
		//@}
	public:
		//! \name Destructors
		//{@
		~TVAPathPricer() {
		}
		//@}
	public:
		//! \name Operator overload
		//{@
		//! Calculate the TVA (includes CVA, DVA and FVA) given the Multi-path
		result_type operator() (const argument_type &multiPath) const {
			const boost::shared_ptr<const ENGINE> engine = engine_.lock();
			QL_REQUIRE(engine, "Engine has been destroyed.");

			/*double cva = 0.0;
			double dva = 0.0;*/
			double fva = 0.0;

			const boost::shared_ptr<const Counterparty> issuer = engine->counterparties()[0];
			const boost::shared_ptr<const Counterparty> investor = engine->counterparties()[1];
			const std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>>
				processes = engine->processList();

			QL_REQUIRE(multiPath.assetNumber() >= 2, "missing issuer/investor path information");

			Time issuerDefaultTime = issuer->getDefaultTime(multiPath[0]);
			Time investorDefaultTime = investor->getDefaultTime(multiPath[1]);

			std::vector<boost::shared_ptr<const MCExposureModel>> tvaModels = engine->models();
			boost::shared_ptr<ShortRateTermStructure> domesticTSModel = engine->domesticShortRateModel();

			double defaultTime = std::min({ investorDefaultTime, issuerDefaultTime, engine->endTime() });
			if (defaultTime == engine->endTime()){
				investorDefaultTime = engine->endTime();
				issuerDefaultTime = engine->endTime();
			}
			TimeGrid tGrid = multiPath[0].timeGrid();
			Time t = 0;
			int i = 0;
			while (t < defaultTime){
				t = tGrid[++i];

				Path domesticInterestRatePath = multiPath[2]; // should initiate the path information
				auto tPos = tGrid.index(t);
				Real shortRate = domesticInterestRatePath[tPos];

				// calculate discount factor D(0, defaultTime) from reference time to the default time
				Real discountFactor = 1.0;
				Time previousTime = tGrid[0];
				for (auto k = 0; k <= tPos; ++k) {
					discountFactor *= std::exp(-domesticInterestRatePath.value(k) * (tGrid[k] - previousTime));
					previousTime = tGrid[k];
				}
				// find target date
				Date targetDate = engine->domesticCurve()->referenceDate()+(int)t*360;
				//Date defaultDate = engine->domesticCurve()->referenceDate();
				/*while (engine->domesticCurve()->dayCounter().yearFraction(referenceDate, defaultDate) < defaultTime) {
					defaultDate += 1;
				}*/


				// construct short rate term structure from default time
				boost::shared_ptr<YieldTermStructure> shortRateDiscountTS = domesticTSModel->GetTermStructure(
					targetDate, engine->endTime(), shortRate);

				double exposure = 0.;

				std::vector<Path> instrumentPaths;
				for (auto model : tvaModels){
					instrumentPaths.push_back(domesticInterestRatePath);
					std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> process = model->instrumentProcess();
					std::vector<std::size_t> pathIndex;
					for (auto p : process){
						auto path = std::find(processes.begin(),processes.end(),p);
						assert(path != processes.end());
						pathIndex.push_back(std::distance(processes.begin(), path));
					}

					for (auto i : pathIndex) {
						instrumentPaths.push_back(multiPath[i]);
					}

					MultiPath instrumentPath(std::move(instrumentPaths));

					exposure += model->exposure(instrumentPath, t, t, Handle<YieldTermStructure>(shortRateDiscountTS)) * discountFactor;
					instrumentPaths.clear();
				}
				fva += issuer->fundingSpread(t)*exposure;
			}

			fva *= tGrid[1]-tGrid[0];

			Array r(1);
			//r[0] = cva;
			//r[1] = dva;
			r[0] = fva;
			return r;
		}
		//@}
	protected:
		const boost::weak_ptr<const ENGINE> engine_;
	};
}