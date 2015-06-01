#pragma once

#include <Calculation/Calculation.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/instrument.hpp>
#include <ql/pricingengine.hpp>
#include <ql/credit/counterparty.hpp>

namespace QuantLib {
	class  VAR : public Instrument {
	public:
		class arguments;
		class results;
		class engine;
	public:
		virtual ~VAR(){}

		const Real CreditVAR() const{
			calculate();
			return VAR_;
		}

	public:
		void setupArguments(PricingEngine::arguments*) const {}
		void fetchResults(const PricingEngine::results*) const;
		bool isExpired() const { return false; }


	protected:
		/*std::vector<boost::shared_ptr<const Instrument>> portfolio_;
		std::vector<boost::shared_ptr<const Counterparty>> counterparties_;
		Handle<YieldTermStructure> discountCurve_;*/
		mutable Real VAR_;
	};

	class VAR::arguments : public virtual PricingEngine::arguments{
	public:
		virtual ~arguments() {
		}
		virtual void validate() const {
			/*QL_REQUIRE(portfolio.size() >= 1, "Portfolio should contains at least one instrument.");
			QL_REQUIRE(portfolio.size() == 2, "Exact two counterparties should be specified.");
			QL_REQUIRE(!discountCurve.empty(), "Discount curve is needed");*/
		}
	public:
		/*std::vector<boost::shared_ptr<const Instrument>> portfolio;
		std::vector<boost::shared_ptr<const Counterparty>> counterparties;
		Handle<YieldTermStructure> discountCurve;*/
	};

	class VAR::results : public Instrument::results{
	public:
		virtual ~results() {
		}
		virtual void reset() {
			Instrument::results::reset();
			VAR = Null<Real>();
		}
	public:
		Real VAR;
	};

	class VAR::engine : public GenericEngine < VAR::arguments, VAR::results >{
	};

	inline void VAR::fetchResults(const PricingEngine::results* r) const{
		const VAR::results* results =
			dynamic_cast<const VAR::results*>(r);
		QL_ENSURE(results != 0,
			"no results returned from pricing engine");

		VAR_ = results->VAR;

		additionalResults_ = results->additionalResults;
	}
}
