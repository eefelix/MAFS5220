/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Xiang, GAO
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
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/models/shortratetermstructure.hpp>

//! \file bvapathpricer.hpp
//! \brief Path pricer class used to calculate the BVA on each Multi-path for BVAEngine class.
//!
namespace QuantLib {
	//! \brief Path pricer class used to calculate the BVA on each Multi-path for BVAEngine class. 
	/*!	This class is the core part for BVAEngine class, which is used to calculate the BVA (includes CVA, DVA)
		on each Multi-path given the instruments (single instrument or portfolio), based on the definition of BVA
		(includes CVA, DVA) in Brigo's book: "Counterparty Credit Risk, Collateral and Funding".

		If the underlying is a portfolio, netting is considered.
		
		\see MCBVAEngine
	*/
 
	template <typename ENGINE>
	class BVAPathPricer : public ENGINE::path_pricer_type {
	public:
		typedef typename ENGINE::path_pricer_type::argument_type argument_type;
		typedef typename ENGINE::path_pricer_type::result_type result_type;

		//! \name Constructors 
		//{@
		BVAPathPricer(const boost::shared_ptr<const ENGINE> &engine)
			: engine_(engine) {
		}
		//@}
	public:
		//! \name Destructors
		//{@
		~BVAPathPricer() {
		}
		//@}
	public:
		//! \name Operator overload
		//{@
		//! Calculate the TVA (includes CVA, DVA and FVA) given the Multi-path
		result_type operator() (const argument_type &multiPath) const {
			const boost::shared_ptr<const ENGINE> engine = engine_.lock();
			QL_REQUIRE(engine, "Engine has been destroyed.");

			double cva = 0.0;
			double dva = 0.0;

			const boost::shared_ptr<const Counterparty> issuer = engine->counterparties()[0];
			const boost::shared_ptr<const Counterparty> investor = engine->counterparties()[1];
			const std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>>
				processes = engine->processList();

			QL_REQUIRE(multiPath.assetNumber() >= 2, "missing issuer/investor path information");

			Time issuerDefaultTime = issuer->getDefaultTime(multiPath[0]);
			Time investorDefaultTime = investor->getDefaultTime(multiPath[1]);

			std::vector<boost::shared_ptr<const MCExposureModel>> tvaModels = engine->models();
			boost::shared_ptr<ShortRateTermStructure> domesticTSModel = engine->domesticShortRateModel();

			if (issuerDefaultTime <= investorDefaultTime && issuerDefaultTime <= engine->endTime()) {
				std::vector<Path> instrumentPaths;

				// get domestic short rate at default time
				Path domesticInterestRatePath = multiPath[2]; // should initiate the path information
				TimeGrid timeGrid = domesticInterestRatePath.timeGrid();

				// calculate discount factor D(0, defaultTime) from reference time to the default time
				// and short rate at default time
				Real discountFactor = 1.0;
				auto defaultPos = timeGrid.begin();
				while (*defaultPos < issuerDefaultTime){
					discountFactor *= std::exp(-domesticInterestRatePath[timeGrid.index(*defaultPos)] * (timeGrid[1] - timeGrid[0]));
					++defaultPos;
				}	
				Real shortRate = domesticInterestRatePath[timeGrid.index(*defaultPos)];
				
				// find default date
				Date defaultDate = engine->domesticCurve()->referenceDate() + (int)issuerDefaultTime * 360;
				/*Date defaultDate = engine->domesticCurve()->referenceDate();
				while (engine->domesticCurve()->dayCounter().yearFraction(referenceDate, defaultDate) < issuerDefaultTime) {
					defaultDate += 1;
				}*/
				// construct short rate term structure from default time
				boost::shared_ptr<YieldTermStructure> shortRateDiscountTS = domesticTSModel->GetTermStructure(
					defaultDate, engine->endTime(), shortRate);

				for (auto model : tvaModels){
					instrumentPaths.push_back(domesticInterestRatePath); // path of domestic rate
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

					cva += model->exposure(instrumentPath, issuerDefaultTime, investorDefaultTime, Handle<YieldTermStructure>(shortRateDiscountTS));
					instrumentPaths.clear();
				}
				cva = (1 - issuer->getRecoveryRate()) * std::max(cva, 0.0) * discountFactor;
			}

			if (investorDefaultTime <= issuerDefaultTime && investorDefaultTime <= engine->endTime()) {
				std::vector<Path> instrumentPaths;

				// get domestic short rate at default time
				Path domesticInterestRatePath = multiPath[2]; // should initiate the path information
				TimeGrid timeGrid = domesticInterestRatePath.timeGrid();
				
				// calculate discount factor D(0, defaultTime) from reference time to the default time
				// and short rate at default time
				Real discountFactor = 1.0;
				auto defaultPos = timeGrid.begin();
				while (*defaultPos < investorDefaultTime){
					discountFactor *= std::exp(-domesticInterestRatePath[timeGrid.index(*defaultPos)] * (timeGrid[1] - timeGrid[0]));
					++defaultPos;
				}
				Real shortRate = domesticInterestRatePath[timeGrid.index(*defaultPos)];

				// find default date
				Date defaultDate = engine->domesticCurve()->referenceDate() + (int)investorDefaultTime*360;
				/*Date defaultDate = engine->domesticCurve()->referenceDate();
				while (engine->domesticCurve()->dayCounter().yearFraction(referenceDate, defaultDate) < issuerDefaultTime) {
					defaultDate += 1;
				}*/
				// construct short rate term structure from default time
				boost::shared_ptr<YieldTermStructure> shortRateDiscountTS = domesticTSModel->GetTermStructure(
					defaultDate, engine->endTime(), shortRate);

				for (auto model : tvaModels){
					instrumentPaths.push_back(domesticInterestRatePath);

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

					dva += model->exposure(instrumentPath, issuerDefaultTime, investorDefaultTime, Handle<YieldTermStructure>(shortRateDiscountTS));
					instrumentPaths.clear();
				}
				dva = -(1 - investor->getRecoveryRate()) * std::min(dva, 0.0) * discountFactor;
			}

			Array r(2);
			r[0] = cva;
			r[1] = dva;
			return r;
		}
		//@}
	protected:
		const boost::weak_ptr<const ENGINE> engine_;
	};
}