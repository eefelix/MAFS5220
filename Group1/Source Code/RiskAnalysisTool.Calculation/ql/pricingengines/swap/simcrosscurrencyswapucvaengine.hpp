//#pragma once
//
//#include <Calculation/Calculation.h>
//#include <ql/instruments/swap/crosscurrencyswap.hpp>
//#include <ql/credit/counterparty.hpp>
//
//namespace QuantLib
//{
//	class _RISKANALYSISTOOL_CALCULATION_API CCSEngine :public CrossCurrencySwap::engine
//	{
//	public:
//		CCSEngine(
//			double exchangeRate,
//			const Date& npvDate,
//			double rd,
//			double rf,
//			boost::shared_ptr<YieldTermStructure>	discountCurve,
//			long nPath = 50000,
//			long nSteps = 1000);
//		~CCSEngine();
//		void calculate() const;
//		static int defaultNumbers;
//	protected:
//		void pathPricer(
//			double & receiveNPV,
//			double & payNPV,
//			const double & fxRate,
//			const Date& startDate) const;
//
//		//boost::shared_ptr<Counterparty>	counterparty_;
//		double exchangeRate_;
//		Date npvDate_;
//		long nPath_, nSteps_;
//		double domesticRate_, foreignRate_;
//		boost::shared_ptr<YieldTermStructure> discountCurve_;
//	};
