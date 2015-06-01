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
#include <boost/make_shared.hpp>
#include <ql/qldefines.hpp>
#include <ql/types.hpp>
#include <ql/math/matrix.hpp>
#include <ql/methods/montecarlo/mctraits.hpp>
#include <ql/math/statistics/arraystatistics.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/pricingengines/mcsimulation.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/pricingengines/ucvaengine.hpp>
#include <ql/methods/montecarlo/multivaluemctraits.hpp>
#include <ql/processes/stochasticprocessarray.hpp>
//! \file mcucvaengine.hpp
//! \brief Abstract base class for UCVA calculation via Monte Carlo Simulation. 
//!
namespace QuantLib {
	//! \brief Abstract base class for UCVA calculation via Monte Carlo Simulation. 
	/*!	This abstract base class is the base Monte Carlo pricing engine class of instrument for UCVA
		calculation.
		
		The kernel calculation part is its path pricer class.

		The user who wants to calculate the UCVA of a instrument should implement the UCVA engine class
		and the path pricer class for the instrument.

		\see UCVAPathPricer
	*/
	template <typename INST, typename RNG = PseudoRandom, typename S = Statistics>
	class MCUCVAEngine : public UCVAEngine<INST>, protected McSimulation < MultiVariate, RNG, S > {
	public:
		typedef typename MultiVariate<RNG>::path_type path_type;
		typedef typename McSimulation<MultiVariate, RNG, S>::stats_type stats_type;
		typedef typename McSimulation<MultiVariate, RNG, S>::path_pricer_type path_pricer_type;
		typedef typename McSimulation<MultiVariate, RNG, S>::path_generator_type path_generator_type;
	public:
		//! \name Constructors & Destructors
		//{@
		MCUCVAEngine(const boost::shared_ptr<const INST> &instrument,
			const boost::shared_ptr<const Counterparty> &issuer,// const boost::shared_ptr<const Counterparty> &investor,
			const Matrix &correlation,
			Time endTime, const Handle<YieldTermStructure>& discountCurve,
			Size timeSteps, Size timeStepsPerYear, bool antitheticVariate,
			Size requiredSamples, Real requiredTolerance, Size maxSamples, BigNatural seed
			) : UCVAEngine(instrument, issuer, endTime, discountCurve), McSimulation(antitheticVariate, false),
			correlation_(correlation),
			timeSteps_(timeSteps), timeStepsPerYear_(timeStepsPerYear),
			requiredSamples_(requiredSamples), requiredTolerance_(requiredTolerance), maxSamples_(maxSamples),
			seed_(seed) {
		}

		virtual ~MCUCVAEngine() {
		}
		//@}
	public:
		//! \name Public interface
		//{@
		//! calculate UCVA by monte carlo simulation and store results in pricing engine
		void calculate() const {
			McSimulation::calculate(requiredTolerance_, requiredSamples_, maxSamples_);
			const S& stats = this->mcModel_->sampleAccumulator();
			Real mean = stats.mean();
			Real errorEstimate = stats.errorEstimate();
			results_.additionalResults["UCVA"] = mean;
			if (RNG::allowsErrorEstimate) {
				results_.additionalResults["UCVA_ErrorEstimate"] = errorEstimate;
				results_.errorEstimate = errorEstimate;
			}
		}
		//! Return shared pointer to stochastic procss array containing underlying process as well as default process 
		boost::shared_ptr<StochasticProcessArray> process() const {
			if (!process_) {
				std::vector<boost::shared_ptr<StochasticProcess1D>> processes;
				processes.push_back(issuer_->createDefaultProcess());
				for (boost::shared_ptr<StochasticProcess1D> p : this->instrumentProcess()) {
					processes.push_back(p);
				}

				process_.reset(new StochasticProcessArray(processes, correlation_));
			}
			return process_;
		}
		//@}
	protected:
		//! \name Protected interface
		//{@
		//! Return time grid used by Monte Carlo Simulation
		virtual TimeGrid timeGrid() const {
			Time maturity = endTime_;

			if (timeSteps_ != Null<Size>()) {
				return TimeGrid(maturity, timeSteps_);
			}
			else if (timeStepsPerYear_ != Null<Size>()) {
				Size steps = (Size)std::round(timeStepsPerYear_* maturity);
				return TimeGrid(maturity, std::max<Size>(static_cast<std::size_t>(std::round(maturity)), steps));
			}
			else {
				QL_FAIL("time steps not specified");
			}
		}

		//! Return shared pointer to path generator who generates sample pathes in simulation
		virtual boost::shared_ptr<path_generator_type> pathGenerator() const {
			Size factors = this->process()->factors();
			TimeGrid grid = timeGrid();
			Size steps = grid.size() - 1;
			typename RNG::rsg_type gen = RNG::make_sequence_generator(factors * steps, seed_);
			return boost::make_shared<path_generator_type>(this->process(), grid, gen, false);
		}
		//@}
		virtual boost::shared_ptr<path_pricer_type> pathPricer() const = 0;
		virtual std::vector<boost::shared_ptr<StochasticProcess1D>> instrumentProcess() const = 0;
	private:
		mutable boost::shared_ptr<StochasticProcessArray> process_;
		const Matrix correlation_;
		const Size timeSteps_, timeStepsPerYear_;
		const Size requiredSamples_, maxSamples_;
		const Real requiredTolerance_;
		mutable BigNatural seed_;
	};
}