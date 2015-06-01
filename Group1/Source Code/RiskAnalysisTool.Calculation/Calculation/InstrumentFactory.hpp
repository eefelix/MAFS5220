#pragma once
#include <Calculation/Utilities/DynamicVisitor.hpp>

namespace QuantLib {
	class Instrument;
}

namespace RiskAnalysisTool {
	namespace Calculation {
#ifdef _M_CEE
		public ref class InstrumentFactory : public Utilities::DynamicVisitor<RiskAnalysisTool::Instruments::Instrument ^> {
		public:
			InstrumentFactory(RiskAnalysisTool::Requests::ComputationRequest^ request, boost::shared_ptr<QuantLib::YieldTermStructure> *curve);
			virtual ~InstrumentFactory();
			!InstrumentFactory();
		public:
			QuantLib::Instrument *CreateQLInstrument(RiskAnalysisTool::Instruments::Instrument ^instrument);
		protected:
			void OnVisit(RiskAnalysisTool::Instruments::InterestRateSwap ^instrument);
			void OnVisit(RiskAnalysisTool::Instruments::EquitySwap ^instrument);
			void OnVisit(RiskAnalysisTool::Instruments::CrossCurrencySwap ^instrument);
		protected:
			QuantLib::Instrument *instrument_;
			RiskAnalysisTool::Requests::ComputationRequest^ request_;
			boost::shared_ptr<QuantLib::YieldTermStructure> *curve_;
		};
#endif
	}
}
