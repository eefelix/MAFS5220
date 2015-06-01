#include <pch.h>
#include "analyticeuropeanucvaengine.hpp"
//
//using namespace QuantLib;
//
//
//AnalyticEuropeanCVAengine::AnalyticEuropeanCVAengine(
//	boost::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> process, 
//	boost::shared_ptr<QuantLib::YieldTermStructure> discountcurve,
//	boost::shared_ptr<Counterparty> c, double rho_) :
//process_(process), discountCurve(discountcurve), C(c), rho(rho_)
//{}
//
//void AnalyticEuropeanCVAengine::calculate() const
//{
//	QL_REQUIRE(arguments_.exercise->type() == QuantLib::Exercise::European,
//		"not an European option");
//
//	boost::shared_ptr<QuantLib::StrikedTypePayoff> payoff =
//		boost::dynamic_pointer_cast<QuantLib::StrikedTypePayoff>(arguments_.payoff);
//	QL_REQUIRE(payoff, "non-striked payoff given");
//
//	long seed1 = QuantLib::SeedGenerator::instance().get();
//	long seed2 = QuantLib::SeedGenerator::instance().get();
//
//	QuantLib::MersenneTwisterUniformRng rnd1(seed1);
//	QuantLib::MersenneTwisterUniformRng rnd2(seed2);
//
//	QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> normrnd1(rnd1);
//	QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> normrnd2(rnd2);
//
//	QuantLib::DayCounter daycounter = process_->riskFreeRate()->dayCounter();
//	double t = daycounter.yearFraction(process_->riskFreeRate()->referenceDate(),arguments_.exercise->lastDate());
//	double t0 = 0.0;
//	double U_0 = process_->x0();
//	double C_0 = C->process_->x0();
//	double C_t = C_0;
//	double U_t = U_0;
//	double h_C = 0;
//	double dw1 = 0, dw2 = 0;
//	double dw2_ = 0;
//	double dt = 0.004; // 0.004=1/250
//
//	double S = 0;
//	size_t n = 10000;
//	std::vector<double> cva(n, 0);
//
//	boost::shared_ptr<FlatForward> constinterest =
//		boost::dynamic_pointer_cast<FlatForward>(discountCurve);
//	double interestrate = discountCurve->forwardRate(t0, t0 + 0.000001, QuantLib::Continuous);
//
//	boost::shared_ptr<FlatForward> constdividend =
//		boost::dynamic_pointer_cast<FlatForward>(*(process_->dividendYield()));
//	double dividend = process_->dividendYield()->forwardRate(t0, t0 + 0.000001, QuantLib::Continuous);
//
//	boost::shared_ptr<BlackConstantVol> constvol =
//		boost::dynamic_pointer_cast<BlackConstantVol>(*(process_->blackVolatility()));
//	double sigma = process_->blackVolatility()->blackVol(0, 10);
//	
//
//	for (size_t i = 0; i < n; ++i)
//	{
//		t0 = 0.0;
//		U_t = U_0;
//
//		do{
//			dw1 = normrnd1.next().value;
//			dw2 = normrnd2.next().value;
//			dw2_ = dw1*rho + sqrt(1 - rho*rho)*dw2;
//
//			C_t = C->process_->evolve(t0, C_t, std::min(dt, t - t0), dw1);
//			U_t = process_->evolve(t0, U_t, std::min(dt, t - t0), dw2_);
//			t0 = std::min(t, dt + t0);
//
//			h_C = C_0*C->process_->dividendYield()->discount(t0) / discountCurve->discount(t0);
//
//			if (C_t < h_C){
//				// counterparty default
//				Date defaultdate = C->process_->riskFreeRate()->referenceDate() + ((int)(t0 * 360));
//				boost::shared_ptr<QuantLib::YieldTermStructure> newdiscountCurve =
//					boost::shared_ptr<QuantLib::YieldTermStructure>(
//					new QuantLib::FlatForward(defaultdate, interestrate, discountCurve->dayCounter()));
//				
//				if (constinterest){
//					boost::shared_ptr<QuantLib::BlackScholesMertonProcess> process(new QuantLib::BlackScholesMertonProcess(
//						QuantLib::Handle<QuantLib::Quote>(new QuantLib::SimpleQuote(U_t)),
//						Handle<YieldTermStructure>(
//						new QuantLib::FlatForward(defaultdate, dividend, discountCurve->dayCounter())),
//						Handle<YieldTermStructure>(
//						new QuantLib::FlatForward(defaultdate, interestrate, discountCurve->dayCounter())),
//						Handle<BlackVolTermStructure>(
//						new QuantLib::BlackConstantVol(defaultdate,TARGET(),sigma,discountCurve->dayCounter()))));
//
//					boost::shared_ptr<QuantLib::AnalyticEuropeanEngine> engine(new QuantLib::AnalyticEuropeanEngine(process));
//
//					QuantLib::EuropeanOption option(payoff, arguments_.exercise);
//					/* the risk free rate should be adapted to the time level tau!!!!!!!
//					for the time being the interest rate is assumed constant, so is ok */
//					option.setPricingEngine(engine);
//
//					// cva is the positive part of NPV
//					cva[i] = std::max(option.NPV(), 0.0);
//					S += cva[i] * discountCurve->discount(t0);
//				}
//			}
//		} while (t0 < t);
//	}
//	results_.value = S / n;
//}