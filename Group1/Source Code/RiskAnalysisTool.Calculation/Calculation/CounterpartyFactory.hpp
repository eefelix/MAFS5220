#pragma once
#include <Calculation/Utilities/DynamicVisitor.hpp>

namespace QuantLib {
	class Counterparty;
}

namespace RiskAnalysisTool {
	namespace Calculation {
#ifdef _M_CEE
		public ref class CounterpartyFactory : public Utilities::DynamicVisitor < RiskAnalysisTool::Requests::ComputationRequest ^ > {
		public:
			CounterpartyFactory(boost::shared_ptr<QuantLib::YieldTermStructure> * curve);
			virtual ~CounterpartyFactory();
			!CounterpartyFactory();
		public:
			std::vector<boost::shared_ptr<QuantLib::Counterparty>> *CreateCounterparties(RiskAnalysisTool::Requests::ComputationRequest ^request);
		protected:
			void OnVisit(RiskAnalysisTool::Requests::TvaRequest ^tvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::BvaRequest ^bvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::UCvaRequest ^ucvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::CreditVaRRequest ^cvarRequest);
		protected:
			boost::shared_ptr<QuantLib::YieldTermStructure> * curve_;
			std::vector<boost::shared_ptr<QuantLib::Counterparty>>* counterparties_;
		};
#endif
	}
}
