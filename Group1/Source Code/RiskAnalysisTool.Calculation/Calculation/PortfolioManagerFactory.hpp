#pragma once
#include <Calculation/Utilities/DynamicVisitor.hpp>

namespace QuantLib {
	class Instrument;
}

namespace RiskAnalysisTool {
	namespace Calculation {
//#ifdef _M_CEE
		public ref class PortfolioManagerFactory : public Utilities::DynamicVisitor < RiskAnalysisTool::Requests::ComputationRequest ^ > {
		public:
			PortfolioManagerFactory();
			virtual ~PortfolioManagerFactory();
			!PortfolioManagerFactory();
		public:
			QuantLib::Instrument *CreatePortfolioManager(RiskAnalysisTool::Requests::ComputationRequest ^request);
		protected:
			void OnVisit(RiskAnalysisTool::Requests::TvaRequest ^tvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::BvaRequest ^bvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::UCvaRequest ^ucvaRequest);
			void OnVisit(RiskAnalysisTool::Requests::CreditVaRRequest ^cvarRequest);
		protected:
			QuantLib::Instrument *portfolioManager_;
		};
//#endif
	}
}
