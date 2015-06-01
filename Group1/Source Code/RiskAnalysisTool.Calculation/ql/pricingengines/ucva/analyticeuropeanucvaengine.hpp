//#pragma once
//
//#include <Calculation/Calculation.h>
//#include <ql/instruments/europeanoption.hpp>
//#include <ql/processes/blackscholesprocess.hpp>
//#include <ql/credit/counterparty.hpp>
//
//namespace QuantLib {
//
//    class _RISKANALYSISTOOL_CALCULATION_API AnalyticEuropeanCVAengine : public QuantLib::EuropeanOption::engine {
//    public:
//        AnalyticEuropeanCVAengine(
//            boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process,
//            boost::shared_ptr<QuantLib::YieldTermStructure> discountcurve,
//            boost::shared_ptr<Counterparty> c,
//            double rho_);
//        void calculate() const;
//    private:
//        boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process_;
//        boost::shared_ptr<QuantLib::YieldTermStructure> discountCurve;
//		boost::shared_ptr<Counterparty> C;
//        double rho;
//
//    };
//}