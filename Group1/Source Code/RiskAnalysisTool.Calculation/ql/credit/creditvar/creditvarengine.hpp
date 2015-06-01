#pragma once
//
//#include <Calculation/Calculation.h>
//#include <ql/errors.hpp>
//#include <memory>
//#include <vector>
//#include <ql/stochasticprocess.hpp>
//#include <ql/instrument.hpp>
//#include <ql/payoff.hpp>
//#include <ql/time/calendar.hpp>
//#include <ql/time/daycounter.hpp>
//#include <ql/processes/kouprocess.hpp>
//
//namespace QuantLib {
//
//    //template<typename __UnderlyingProcess_Type>
//    class _RISKANALYSISTOOL_CALCULATION_API CreditVaRComputationEngine {
//    public:
//        CreditVaRComputationEngine();
//        virtual ~CreditVaRComputationEngine();
//
//        CreditVaRComputationEngine& setRiskFreeRate(double riskFreeRate);
//        CreditVaRComputationEngine& setDefaultBarrier(double newDefaultBarrier);
//        CreditVaRComputationEngine& setConfidence(double newConfidence = 0.999);
//        CreditVaRComputationEngine& setCorrelation(double newCorrelation);
//        CreditVaRComputationEngine& setLGD(double newLGD);
//        CreditVaRComputationEngine& setSteps(long newSteps);
//        CreditVaRComputationEngine& setNumberOfPath(long numberOfPath);
//        CreditVaRComputationEngine& setReferenceDate(QuantLib::Date newReferenceDate);
//        CreditVaRComputationEngine& setVarDate(QuantLib::Date newVarDate);
//        CreditVaRComputationEngine& setMaturityDate(QuantLib::Date newMaturityDate);
//        CreditVaRComputationEngine& setDayCounter(QuantLib::DayCounter newDayCounter);
//        CreditVaRComputationEngine& setCalendar(QuantLib::Calendar newCalendar);
//        CreditVaRComputationEngine& setUnderlyingProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newUnderlyingProcess);
//        CreditVaRComputationEngine& setUnderlyingQProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newUnderlyingQProcess);
//        CreditVaRComputationEngine& setDefaultProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newDefaultProcess);
//        CreditVaRComputationEngine& setDefaultQProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newDefaultQProcess);
//        //CreditVaRComputationEngine& setInstrument(const QuantLib::Instrument& newInstrument);
//        CreditVaRComputationEngine& setPayoff(std::shared_ptr<const QuantLib::Payoff> newPayoff);
//
//        // to do
//        void calculate(
//            const std::vector<QuantLib::Time>& cdsMaturity,
//            char* defaultModel,
//            double K = 100 /*strike price for option*/,
//            QuantLib::Option::Type Type = QuantLib::Option::Call  /*Type for option*/);
//        double getVaR() const;
//        double getES() const;
//        //
//        boost::shared_ptr<std::vector<double>> lossDistribution_;
//        /*
//        double(*DefaultTimeDistFunc_)(const QuantLib::Handle<QuantLib::BlackVolTermStructure>& vol,
//        double from, double to, double initialState, double barrier);
//        */
//    private:
//        double riskFreeRate_;
//        double confidence_;
//        double defaultBarrier_;
//        double correlation_;
//        double LGD_;
//        long steps_;
//        long numberOfPath_;
//
//        QuantLib::Date referenceDate_;
//        QuantLib::Date maturityDate_;
//        QuantLib::Date varDate_;
//        QuantLib::DayCounter dayCounter_;
//        QuantLib::Calendar calendar_;
//
//        //std::shared_ptr<const QuantLib::Instrument> instrument_;
//        std::shared_ptr<const QuantLib::Payoff> payoff_;
//        //every instrument has its own special payoff inherit from Payoff Class
//        std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> underlyingProcess_;
//        std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> underlyingQProcess_;
//        std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> defaultProcess_;
//        std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> defaultQProcess_;
//
//        double defaultTime(
//            const std::vector<QuantLib::Time>& cdsMaturity,
//            double defaultProb,
//            const QuantLib::Handle<QuantLib::BlackVolTermStructure>& vol,
//            double from,
//            double initialState,
//            double epsilon = 1e-6) const;
//        std::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> CreateNewUnderlyingProcess(
//            std::string processName,
//            double newState,
//            QuantLib::Date newDate,
//            double interestRate,
//            double dividend,
//            double volatility,
//            double jumpIntensity,
//            double posProbability,
//            double posJumpMean,
//            double negJumpMean
//            ) const;
//    };
//}
