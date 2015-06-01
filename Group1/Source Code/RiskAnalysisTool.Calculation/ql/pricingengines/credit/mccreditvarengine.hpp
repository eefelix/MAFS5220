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
#include <ql/methods/montecarlo/multivaluemctraits.hpp>
#include <ql/pricingengines/credit/creditvarengine.hpp>
#include <ql/pricingengines/credit/creditvarpathpricer.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <ql/instruments/credit/var.hpp>
#include <ql/processes/stochasticprocessarray.hpp>
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/models/shortratetermstructure.hpp>
#include <ql/models/shortrate/onefactormodel.hpp>

//! \file mccreditvarengine.hpp
//! \brief Abstract base class for Credit VaR calculation via Monte Carlo Simulation. 
//!
//! This class is used to calculate Credit VaR via Monte Carlo Simulation.
//! 
//!	\date 2015-04-15
namespace QuantLib {
	//! \brief Abstract base class for Credit VaR calculation via Monte Carlo Simulation. 
	//!
	//! This class is used to calculate Credit VaR via Monte Carlo Simulation.
	template <typename RNG = PseudoRandom, typename S = GenericArrayStatistics<GaussianStatistics>>
	class MCCreditVaREngine : public CreditVaREngine<VAR> /*TVA::engine*/
		, public McSimulation < MultiValueMultiVariate, RNG, S >, public boost::enable_shared_from_this < MCCreditVaREngine<RNG, S> > {
	public:
		typedef typename MultiValueMultiVariate<RNG>::path_type path_type;
		typedef typename McSimulation<MultiValueMultiVariate, RNG, S>::stats_type stats_type;
		typedef typename McSimulation<MultiValueMultiVariate, RNG, S>::path_pricer_type path_pricer_type;
		typedef typename McSimulation<MultiValueMultiVariate, RNG, S>::path_generator_type path_generator_type;
	public:
		//! \name Constructors & Destructors
		//{@
		MCCreditVaREngine(
			const std::vector<boost::shared_ptr<const Counterparty>> &counterparties,
			const std::vector<boost::shared_ptr<const MCExposureModel>> &instrumentModels,
			const Matrix &correlationMatrix,
			const std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>>& processes,
			const Handle<YieldTermStructure> &domseticCurve,
			const boost::shared_ptr<OneFactorAffineModel>& domesticCIRModel,
			Size timeSteps, Size timeStepsPerYear, bool antitheticVariate,
			Size requiredSamples, Real requiredTolerance, Size maxSamples,
			Time endTime, BigNatural seed = SeedGenerator::instance().get()
			) : CreditVaREngine(endTime, counterparties, instrumentModels, domseticCurve, domesticCIRModel), McSimulation(antitheticVariate, false),
			correlation_(correlationMatrix), processes_(processes),
			timeSteps_(timeSteps), timeStepsPerYear_(timeStepsPerYear),
			requiredSamples_(requiredSamples), requiredTolerance_(requiredTolerance), maxSamples_(maxSamples),
			seed_(seed) {

		}

		virtual ~MCCreditVaREngine() {
		}
		//@}

	public:
		//! \name Public interface
		//{@
		//! calculate bva by monte carlo simulation and store results in pricing engine
		void calculate() const {
			McSimulation::calculate(requiredTolerance_, requiredSamples_, maxSamples_);
			const S& stats = this->mcModel_->sampleAccumulator();

			Array cvar = stats.percentile(0.99);
			Array errorEstimate = stats.errorEstimate();

			results_.VAR = cvar[0];

			results_.value = results_.VAR;

			if (RNG::allowsErrorEstimate) {
				results_.errorEstimate = errorEstimate[0];
			}
		}

		//! Return shared pointer to stochastic procss array containing underlying process as well as default process
		boost::shared_ptr<StochasticProcessArray> process() const {
			if (!process_) {
				/*std::vector<boost::shared_ptr<StochasticProcess1D>> processes;
				processes.push_back(counterparties_[0]->createDefaultProcess());
				processes.push_back(counterparties_[1]->createDefaultProcess());
				boost::shared_ptr<QuantLib::CIRprocess> domesticRateProcess = boost::make_shared<QuantLib::CIRprocess>(
				domesticCIRModel_->params()[0], domesticCIRModel_->params()[1], domesticCIRModel_->params()[2], domesticCIRModel_->params()[3]);
				processes.push_back(domesticRateProcess);
				int i = 3;
				for (boost::shared_ptr<StochasticProcess1D> p : this->instrumentProcess()) {
				processes.push_back(p);
				if (processToPath_.find(p) == processToPath_.end()){
				processToPath_.emplace(std::make_pair(p, i));
				++i;
				}
				}
				*/
				process_.reset(new StochasticProcessArray(processes_, correlation_));
			}
			return process_;
		}

		const std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>>& processList() const{
			return processes_;
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
				return TimeGrid(maturity, std::max<Size>(static_cast<Size>(std::round(maturity)), steps));
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

		virtual boost::shared_ptr<path_pricer_type> pathPricer() const {
			return boost::make_shared<CreditVaRPathPricer<MCCreditVaREngine<RNG, S>>>(shared_from_this());
		}

		//! Return shared pointer to stochastic procss array containing underlying process
		/*virtual std::vector<boost::shared_ptr<StochasticProcess1D>> instrumentProcess() const{
			std::vector<boost::shared_ptr<StochasticProcess1D>> processes;
			for (boost::shared_ptr<const MCExposureModel> model : CreditVarModels_)
				for (boost::shared_ptr<StochasticProcess1D> p : model->instrumentProcess())
					processes.push_back(p);


			return processes;
		}*/
		//@}
	private:
		mutable boost::shared_ptr<StochasticProcessArray> process_;
		std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> processes_;
		const QuantLib::Matrix correlation_;
		const Size timeSteps_, timeStepsPerYear_;
		const Size requiredSamples_, maxSamples_;
		const Real requiredTolerance_;
		mutable BigNatural seed_;
	};
}