#pragma once
#include <Calculation/Utilities/DynamicVisitor.hpp>

namespace QuantLib {
	class Instrument;
	class MCExposureModel;
}

namespace RiskAnalysisTool {
	namespace Calculation {
#ifdef _M_CEE
		public ref class ModelFactory : public Utilities::DynamicVisitor < RiskAnalysisTool::Instruments::Instrument ^ > {
		public:
			ModelFactory(RiskAnalysisTool::Requests::ComputationRequest^ request, 
				boost::shared_ptr<QuantLib::YieldTermStructure>* curve, boost::shared_ptr<QuantLib::Instrument>* inst);
			virtual ~ModelFactory();
			!ModelFactory();
		public:
			QuantLib::MCExposureModel *CreateQLModel(RiskAnalysisTool::Instruments::Instrument ^instrument);
		protected:
			void OnVisit(RiskAnalysisTool::Instruments::InterestRateSwap ^instrument);
			void OnVisit(RiskAnalysisTool::Instruments::EquitySwap ^instrument);
			void OnVisit(RiskAnalysisTool::Instruments::CrossCurrencySwap ^instrument);
		protected:
			QuantLib::MCExposureModel *model_;
			RiskAnalysisTool::Requests::ComputationRequest^ request_;
			boost::shared_ptr<QuantLib::YieldTermStructure> *curve_;
			boost::shared_ptr<QuantLib::Instrument> *instrument_;
		};
#endif
	}
}
