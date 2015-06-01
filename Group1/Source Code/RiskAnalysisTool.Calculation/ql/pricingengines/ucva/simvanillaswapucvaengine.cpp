#include <pch.h>
//#include <algorithm>
//#include "simvanillaswapucvaengine.hpp"
//
//using namespace QuantLib;
//
//
//namespace{
//	double CIR_h(double a, double sigma){
//		return sqrt(a*a + 2. * sigma*sigma);
//	}
//	double CIR_A(double a, double b, double sigma, double t0, double T){
//		return pow(2. * CIR_h(a, sigma)*exp((a + CIR_h(a, sigma))*(T - t0) / 2.) /
//			(2. * CIR_h(a, sigma) + (a + CIR_h(a, sigma))*(exp((T - t0)*CIR_h(a, sigma)) - 1.)), 2. * a*b / (sigma*sigma));
//	}
//	double CIR_B(double a, double b, double sigma, double t0, double T){
//		return (2. * (exp((T - t0)*CIR_h(a, sigma)) - 1.)) /
//			(2. * CIR_h(a, sigma) + (a + CIR_h(a, sigma))*(exp((T - t0)*CIR_h(a, sigma)) - 1.));
//	}
//	double CIR_Price(double a, double b, double sigma, double t0, double T, double interestRate){
//		return CIR_A(a, b, sigma, t0, T)*exp(-CIR_B(a, b, sigma, t0, T)*interestRate);
//	}
//}
//
//IRSCVAEngine::IRSCVAEngine(
//	DayCounter &daycounter,
//	Date &npvDate,
//	boost::shared_ptr<Counterparty> Counterparty,
//	boost::shared_ptr<CIRprocess> Underlyingprocess,
//	long numberofpath) :
//	daycounter_(daycounter), npvDate_(npvDate),
//	Counterparty_(Counterparty), Underlyingprocess_(Underlyingprocess),
//	numberofpath_(numberofpath)
//{}
//
//void IRSCVAEngine::calculate() const{
//	QL_REQUIRE(Counterparty_ != nullptr,
//		"Counterparty is empty");
//	results_.value = 0.0;
//	results_.valuationDate = npvDate_;
//
//	std::vector<QuantLib::Date> fixedDate_ = arguments_.fixedPayDates;
//	std::vector<QuantLib::Date> floatDate_ = arguments_.floatingPayDates;
//	// transfer settlement Dates to Times
//	std::vector<Time> fixedTime_;
//	for (int i = 0; i != fixedDate_.size(); ++i){
//		fixedTime_.push_back(daycounter_.yearFraction(npvDate_, fixedDate_[i]));
//	}
//	std::vector<Time> floatTime_;
//	for (int j = 0; j != floatDate_.size(); ++j){
//		floatTime_.push_back(daycounter_.yearFraction(npvDate_, floatDate_[j]));
//	}
//
//	QuantLib::BigInteger seed1 = QuantLib::SeedGenerator::instance().get();
//	QuantLib::MersenneTwisterUniformRng URng1(seed1);
//	QuantLib::BigInteger seed2 = QuantLib::SeedGenerator::instance().get();
//	QuantLib::MersenneTwisterUniformRng URng2(seed2);
//	QuantLib::BoxMullerGaussianRng < QuantLib::MersenneTwisterUniformRng > dffusionRng(URng2);
//
//	for (int i = 0; i != numberofpath_; ++i){
//		Probability defaultProb = URng1.next().value;
//		Time tau = Counterparty_->defaultTimeGenerator(defaultProb);
//		if (tau <= fixedTime_.back()){
//			double interestRate = 0.;
//			//determine the remaining float/fixed leg
//			//construct the discount factor
//			int remaining_float = 0;
//			for (int j = 0; j != floatTime_.size(); ++j){
//				if (tau <= floatTime_[j]){
//					remaining_float = j;
//					break;
//				}
//			}
//			int remaining_fixed = 0;
//			for (int j = 0; j != fixedTime_.size(); ++j){
//				if (tau <= fixedTime_[j]){
//					remaining_fixed = j;
//					break;
//				}
//			}
//			//evolue interestRate to default time tau
//			interestRate = Underlyingprocess_->evolve(0, Underlyingprocess_->x0(), tau, dffusionRng.next().value);
//
//			std::vector<Real> discountfactor1_;
//			std::vector<Real> discountfactor2_;
//			for (int l = remaining_float; l != floatTime_.size(); ++l){
//				discountfactor1_.push_back(CIR_Price(Underlyingprocess_->revertSpeed(), Underlyingprocess_->revertLevel(),
//					Underlyingprocess_->volatility(), tau, floatTime_[l], interestRate));
//			}
//			for (int l = remaining_fixed; l != fixedTime_.size(); ++l){
//				discountfactor2_.push_back(CIR_Price(Underlyingprocess_->revertSpeed(), Underlyingprocess_->revertLevel(),
//					Underlyingprocess_->volatility(), tau, fixedTime_[l], interestRate));
//			}
//
//			//calculate the NPV(tau)+
//			double NPVtau = 0.;
//			NPVtau += arguments_.nominal * (1. - discountfactor1_[0]);
//			NPVtau += arguments_.nominal * arguments_.floatingSpreads[0] * discountfactor1_[0];
//			for (int m = 1; m != discountfactor1_.size(); ++m){
//				NPVtau += arguments_.nominal * (discountfactor1_[m - 1] - discountfactor1_[m]);
//				NPVtau += arguments_.nominal * arguments_.floatingSpreads[m] * discountfactor1_[m];
//			}
//			for (int m = 0; m != discountfactor2_.size(); ++m){
//				NPVtau -= (arguments_.nominal* arguments_.fixedCoupons[m]) * discountfactor2_[m];
//			}
//
//			if (arguments_.type == QuantLib::VanillaSwap::Payer){
//				NPVtau = std::max(NPVtau, 0.);
//			}
//			else{
//				NPVtau = std::max(-NPVtau, 0.);
//			}
//			results_.value += (1 - Counterparty_->getRecoveryRate()) * NPVtau *
//				CIR_Price(Underlyingprocess_->revertSpeed(), Underlyingprocess_->revertLevel(),
//				Underlyingprocess_->volatility(), 0, tau, Underlyingprocess_->x0());
//		}
//	}
//	results_.value = results_.value / numberofpath_;
//}
