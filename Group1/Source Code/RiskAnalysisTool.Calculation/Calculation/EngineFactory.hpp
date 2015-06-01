#pragma once
#include <Calculation/Utilities/DynamicVisitor.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <Calculation/Utilities/cliext.hpp>

namespace QuantLib {
	class PricingEngine;
	class MCExposureModel;
}

namespace RiskAnalysisTool {
	namespace Calculation {
#ifdef _M_CEE
		public ref class EngineFactory : public Utilities::DynamicVisitor < RiskAnalysisTool::Requests::ComputationRequest ^ > {
		public:
			EngineFactory(boost::shared_ptr<QuantLib::YieldTermStructure> * curve, std::vector<boost::shared_ptr<const QuantLib::Counterparty>>* counterparties,
				std::vector<boost::shared_ptr<const QuantLib::MCExposureModel>>* models, boost::shared_ptr<QuantLib::CoxIngersollRoss>* domesticModel);
			virtual ~EngineFactory();
			!EngineFactory();
		public:
			boost::shared_ptr<QuantLib::PricingEngine> CreateEngine(RiskAnalysisTool::Requests::ComputationRequest ^request);
		protected:
			void OnVisit(RiskAnalysisTool::Requests::TvaRequest ^tvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::BvaRequest ^bvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::UCvaRequest ^ucvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::CreditVaRRequest ^cvarRequest);
		protected:
			cliext::shared_ptr<QuantLib::PricingEngine> engine_;
			boost::shared_ptr<QuantLib::YieldTermStructure> * curve_;
			std::vector<boost::shared_ptr<const QuantLib::Counterparty>>* counterparties_;
			std::vector<boost::shared_ptr<const QuantLib::MCExposureModel>>* models_;
			boost::shared_ptr<QuantLib::CoxIngersollRoss>* domesticModel_;
		};
		#endif
	}
}
