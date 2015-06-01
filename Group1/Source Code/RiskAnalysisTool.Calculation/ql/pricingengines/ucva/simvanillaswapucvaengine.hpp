//#pragma once
//
//#include <Calculation/Calculation.h>
//#include <ql/instruments/vanillaswap.hpp>
//#include <ql/pricingengine.hpp>
//#include <ql/credit/counterparty.hpp>
//#include <ql/processes/cirprocess.hpp>
//
//namespace QuantLib {
//
//    class _RISKANALYSISTOOL_CALCULATION_API IRSCVAEngine : public QuantLib::VanillaSwap::engine {
//    public:
//		IRSCVAEngine(
//			DayCounter &daycounter,
//			Date &npvDate,
//			boost::shared_ptr<Counterparty> Counterparty,
//			boost::shared_ptr<CIRprocess> Underlyingprocess,
//			long numberofpath);
//
//        void calculate() const;
//
//    private:
//		boost::shared_ptr<CIRprocess> Underlyingprocess_;
//		boost::shared_ptr<Counterparty> Counterparty_;
//		QuantLib::Date npvDate_;
//		QuantLib::DayCounter daycounter_;
//		long numberofpath_;
//    };
//
//}