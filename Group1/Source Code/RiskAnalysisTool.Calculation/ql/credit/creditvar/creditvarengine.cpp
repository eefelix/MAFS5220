#include <pch.h>
//#include <algorithm>
//#include "creditvarengine.hpp"
//
//using namespace QuantLib;
//
//
//namespace{
//	double VarianceIntegral(
//		const std::vector<QuantLib::Time>& cdsMaturity, 
//		const QuantLib::Handle<QuantLib::BlackVolTermStructure>& vol, 
//		double from, double to);
//	
//	double AT1PProb(
//		const std::vector<QuantLib::Time>& cdsMaturity, 
//		const QuantLib::Handle<QuantLib::BlackVolTermStructure>& vol,
//		double from,double to, double initialState, double barrier){
//		//double variance = 0.04*(to-from);//std::pow(vol->blackVol(vol->maxTime(), vol->minStrike(), true), 2)*(to-from);
//		double variance = VarianceIntegral(cdsMaturity,vol, from, to);
//		double probability = 0.;
//		double d1 = (std::log(initialState / barrier) - 0.5 * variance) / std::sqrt(variance);
//		double d2 = -d1 - std::sqrt(variance);
//		QuantLib::CumulativeNormalDistribution normalDist(0,1);
//		return probability = 1. - normalDist(d1) + initialState / barrier*normalDist(d2);
//	}
//
//	double VarianceIntegral(
//		const std::vector<QuantLib::Time>& cdsMaturity, 
//		const QuantLib::Handle<QuantLib::BlackVolTermStructure>& vol,
//		double from, double to){
//
//		if (from == to)
//			return 0;
//
//		double sum = 0.0;
//		int indF = 0, indT = 0;
//		// find the first cds maturity greater than from
//		while (cdsMaturity[indF] <= from)
//			++indF;
//
//		// find the last cds maturity less than to
//		indT = indF;
//		while (cdsMaturity[indT] < from)
//			++indT;
//		--indT;
//		
//		// from and to in the same interval
//		if (indT<indF)
//			return (to - from)*std::pow(vol->blackVol(cdsMaturity[indF], 25, true), 2);
//
//		// add the area in the ends
//		sum += (cdsMaturity[indF] - from)*std::pow(vol->blackVol(cdsMaturity[indF], 25, true), 2);
//		sum += (to-cdsMaturity[indT])*std::pow(vol->blackVol(cdsMaturity[indT+1], 25, true), 2);
//
//		for (int i = indF+1; i <= indT; ++i){
//			sum += 
//				(cdsMaturity[i] - cdsMaturity[i-1])
//				*std::pow(vol->blackVol(cdsMaturity[i], 25, true), 2);
//		}
//
//		return sum;
//	}
//}
//CreditVaRComputationEngine::CreditVaRComputationEngine() {
//	riskFreeRate_ = 0.02;
//	defaultBarrier_ = -1.;
//	confidence_ = 0.;
//	steps_ = 500;
//	numberOfPath_ = 500;
//
//	referenceDate_ = QuantLib::Date();
//	maturityDate_ = QuantLib::Date();
//	varDate_ = QuantLib::Date();
//	dayCounter_ = QuantLib::Business252();
//	calendar_ = QuantLib::TARGET();
//
//	underlyingProcess_ = nullptr;
//	underlyingQProcess_ = nullptr;
//	defaultProcess_ = nullptr;
//	defaultQProcess_ = nullptr;
//	payoff_ = nullptr;
//	lossDistribution_ = boost::make_shared<std::vector<double>>();
//}
//
//CreditVaRComputationEngine::~CreditVaRComputationEngine() {
//	// Empty
//}
//
//void CreditVaRComputationEngine::calculate(
//	const std::vector<QuantLib::Time>& cdsMaturity, 
//	char* defalutModel,
//	double K,
//	QuantLib::Option::Type type) {
//	
//	std::string ProcessName = typeid(*underlyingProcess_).name();
//	std::string KouProcessStr = typeid(QuantLib::KouProcess).name();
//
//
//	lossDistribution_->reserve(numberOfPath_);
//	double varTime_ = dayCounter_.yearFraction(referenceDate_, varDate_);
//	double productMaturity_ = dayCounter_.yearFraction(referenceDate_, maturityDate_);
//
//	QuantLib::BigInteger uSeed = QuantLib::SeedGenerator::instance().get();
//	QuantLib::MersenneTwisterUniformRng uURng(uSeed);
//	QuantLib::BoxMullerGaussianRng < QuantLib::MersenneTwisterUniformRng > uRng(uURng);
//
//	QuantLib::BigInteger dSeed = QuantLib::SeedGenerator::instance().get();
//	QuantLib::MersenneTwisterUniformRng dURng(dSeed);
//
//	// get jump parameters
//
//	double jInt = 2.0;
//	double pJM = 3.0;
//	double nJM = 1.0;
//	double pProb = 0.5;
//
//	if (ProcessName == KouProcessStr){
//		std::shared_ptr<const QuantLib::KouProcess> underlyingProcess
//			= std::dynamic_pointer_cast<const QuantLib::KouProcess>(underlyingProcess_);
//
//		jInt = underlyingProcess->jumpIntensity();
//		pJM = underlyingProcess->posJumpMean();
//		nJM = underlyingProcess->negJumpMean();
//		pProb = underlyingProcess->posProbability();
//	}
//
//	// underlying model parameter
//	double interestRate = underlyingQProcess_->riskFreeRate()->zeroRate(1, QuantLib::Continuous);
//	double dividend = underlyingQProcess_->dividendYield()->zeroRate(1, QuantLib::Continuous);
//	double volatility = underlyingQProcess_->blackVolatility()->blackVol(1, K);
//
//	QuantLib::Handle<QuantLib::BlackVolTermStructure> dVol = defaultQProcess_->blackVolatility();
//	double PriceAtVaRTime = 0.;
//	double PriceAtCurrent = 0.;
//
//	// calculate price at current time
//	for (int i = 0; i < 10000; ++i){
//		double dProb = dURng.next().value;
//		double V0 = defaultQProcess_->x0();
//		double tau = defaultTime(cdsMaturity,dProb,dVol,0,V0);
//
//		if (tau < productMaturity_){
//			double state = underlyingQProcess_->evolve(0, underlyingProcess_->x0(), tau, tau*uRng.next().value);
//
//			QuantLib::Date defaultDate = varDate_ + 1;
//			while (dayCounter_.yearFraction(varDate_, defaultDate) < tau){
//				++defaultDate;
//			}
//
//			std::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> subProcess
//				= CreateNewUnderlyingProcess(
//				ProcessName,
//				state,
//				defaultDate,
//				interestRate,
//				dividend,
//				volatility,
//				jInt,
//				pProb,
//				pJM,
//				nJM);
//
//			std::shared_ptr<QuantLib::AnalyticEuropeanEngine> engine
//				= std::make_shared<QuantLib::AnalyticEuropeanEngine>(RiskAnalysisTool::Calculation::make_boost_ptr(subProcess));
//			std::shared_ptr<QuantLib::VanillaOption> option = std::make_shared<QuantLib::VanillaOption>(
//				boost::make_shared<QuantLib::PlainVanillaPayoff>(type, K),
//				boost::make_shared<QuantLib::EuropeanExercise>(defaultDate)
//				);
//			option->setPricingEngine(RiskAnalysisTool::Calculation::make_boost_ptr(engine));
//
//			PriceAtVaRTime += option->NPV()*std::exp(-interestRate*tau)*(1 - LGD_);
//		}
//		else{
//			// generate the state of underlying asset without jump
//			double state =
//				underlyingQProcess_->evolve(0, underlyingQProcess_->x0(), productMaturity_, productMaturity_*uRng.next().value);
//			PriceAtCurrent
//				+= (*payoff_)(state)*std::exp(-interestRate*productMaturity_);
//		}	
//	}
//
//	PriceAtCurrent /= 10000;
//
//	for (int path = 0; path < numberOfPath_; ++path){
//
//		// generate defalut time tau
//		double dProb = dURng.next().value;
//		double V0 = defaultQProcess_->x0();
//		double tau = defaultTime(cdsMaturity,dProb, dVol, 0, V0);
//
//		// generate the state of underlying asset without jump
//		double t = tau > varTime_ ? varTime_ : tau;
//		double state = underlyingProcess_->evolve(0, underlyingProcess_->x0(), t, t*uRng.next().value);
//
//		QuantLib::Date newDate = referenceDate_ + 1;
//		while (dayCounter_.yearFraction(referenceDate_, newDate) < t){
//			++newDate;
//		}
//
//		std::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> newProcess
//			= CreateNewUnderlyingProcess(
//			ProcessName,
//			state,
//			newDate,
//			interestRate,
//			dividend,
//			volatility,
//			jInt,
//			pProb,
//			pJM,
//			nJM);
//
//		if (tau < varTime_){
//			std::shared_ptr<QuantLib::AnalyticEuropeanEngine> engine
//				= std::make_shared<QuantLib::AnalyticEuropeanEngine>(RiskAnalysisTool::Calculation::make_boost_ptr(newProcess));
//			std::shared_ptr<QuantLib::VanillaOption> option = std::make_shared<QuantLib::VanillaOption>(
//				boost::make_shared<QuantLib::PlainVanillaPayoff>(type, K),
//				boost::make_shared<QuantLib::EuropeanExercise>(newDate)
//				);
//			option->setPricingEngine(RiskAnalysisTool::Calculation::make_boost_ptr(engine));
//
//			PriceAtVaRTime = option->NPV()*std::exp(interestRate*(varTime_ - tau))*(1 - LGD_);
//		}
//		else{
//			for (int i = 0; i < steps_; ++i){
//				double newDefaultProb = dURng.next().value;
//				double newV = defaultProcess_->evolve(0, defaultProcess_->x0(), varTime_, uRng.next().value*varTime_);
//				double temp = defaultTime(cdsMaturity,dProb, dVol, 0, defaultQProcess_->x0());
//				while (temp < varTime_){
//					newV = defaultProcess_->evolve(0, defaultProcess_->x0(), varTime_, uRng.next().value*varTime_);
//					temp = defaultTime(cdsMaturity,dProb, dVol, 0, defaultProcess_->x0());
//				}
//				double newTau = defaultTime(cdsMaturity,newDefaultProb, dVol, varTime_, newV);
//
//				if (newTau < productMaturity_){
//					double newState = newProcess->evolve(0, newProcess->x0(), newTau, newTau*uRng.next().value);
//
//					QuantLib::Date defaultDate = varDate_ + 1;
//					while (dayCounter_.yearFraction(varDate_, defaultDate) < t){
//						++defaultDate;
//					}
//
//					std::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> subProcess
//						= CreateNewUnderlyingProcess(
//						ProcessName,
//						newState,
//						defaultDate,
//						interestRate,
//						dividend,
//						volatility,
//						jInt,
//						pProb,
//						pJM,
//						nJM);
//
//					std::shared_ptr<QuantLib::AnalyticEuropeanEngine> engine
//						= std::make_shared<QuantLib::AnalyticEuropeanEngine>(RiskAnalysisTool::Calculation::make_boost_ptr(subProcess));
//					std::shared_ptr<QuantLib::VanillaOption> option = std::make_shared<QuantLib::VanillaOption>(
//						boost::make_shared<QuantLib::PlainVanillaPayoff>(type, K),
//						boost::make_shared<QuantLib::EuropeanExercise>(defaultDate)
//						);
//					option->setPricingEngine(RiskAnalysisTool::Calculation::make_boost_ptr(engine));
//
//					PriceAtVaRTime += option->NPV()*std::exp(interestRate*(varTime_ - tau))*(1 - LGD_);
//				}
//				else{
//					// generate the state of underlying asset without jump
//					double newState =
//						newProcess->evolve(0, newProcess->x0(), productMaturity_ - varTime_, (productMaturity_ - varTime_)*uRng.next().value);
//					PriceAtVaRTime
//						+= (*payoff_)(newState)*std::exp(interestRate*(varTime_ - productMaturity_));
//				}
//			}
//			PriceAtVaRTime /= steps_;
//		}
//		lossDistribution_->push_back(PriceAtCurrent - PriceAtVaRTime);
//	}
//	
//	std::stable_sort(lossDistribution_->begin(), lossDistribution_->end());
//
//}
//
//double CreditVaRComputationEngine::getVaR() const {
//	// calculate VaR according to confidence
//	int quantile = (int)numberOfPath_*confidence_;
//
//	return lossDistribution_->at(quantile);
//}
//
//double CreditVaRComputationEngine::getES() const {
//	// calculate Excepted shortfall
//	int quantile = (int)numberOfPath_*confidence_;
//
//	std::vector<double>::iterator lossDistItr = lossDistribution_->begin();
//	for (int i = 0; i <= quantile; ++i)
//		++lossDistItr;
//	return std::accumulate(lossDistItr,
//		lossDistribution_->end(), 0.0)
//		/ (numberOfPath_ - quantile - 1);
//}
//
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setRiskFreeRate(double newRiskFreeRate) {
//	riskFreeRate_ = newRiskFreeRate;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setDefaultBarrier(double newDefaultBarrier) {
//	QL_REQUIRE(newDefaultBarrier >= 0 && newDefaultBarrier < 1, "Default Barrier must fall into (0,1)!");
//	defaultBarrier_ = newDefaultBarrier;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setConfidence(double newConfidence) {
//	QL_REQUIRE(newConfidence > 0 && newConfidence < 1, "Default Barrier must fall into (0,1)!");
//	confidence_ = newConfidence;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setReferenceDate(QuantLib::Date newReferenceDate) {
//	referenceDate_ = newReferenceDate;
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setVarDate(QuantLib::Date newVarDate) {
//	varDate_ = newVarDate;
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setMaturityDate(QuantLib::Date newMaturityDate) {
//	maturityDate_ = newMaturityDate;
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setDayCounter(QuantLib::DayCounter newDayCounter) {
//	dayCounter_ = newDayCounter;
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setCalendar(QuantLib::Calendar newCalendar) {
//	calendar_ = newCalendar;
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setCorrelation(double newCorrelation) {
//	QL_REQUIRE(newCorrelation <= 1 && newCorrelation >= -1, "error correlation!");
//	correlation_ = newCorrelation;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setLGD(double newLGD) {
//	QL_REQUIRE(newLGD >= 0 && newLGD <= 1, "error LGD!");
//	LGD_ = newLGD;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setNumberOfPath(long newNumberOfPath) {
//	QL_REQUIRE(newNumberOfPath > 0, "Product Maturity must be positive!");
//	numberOfPath_ = newNumberOfPath;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setSteps(long newSteps) {
//	QL_REQUIRE(newSteps > 0, "Product Maturity must be positive!");
//	steps_ = newSteps;
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setUnderlyingProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newUnderlyingProcess) {
//	underlyingProcess_ = newUnderlyingProcess;//reset(newUnderlyingProcess);
//	QL_REQUIRE(underlyingProcess_ != nullptr, "Underlying process is null!");
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setUnderlyingQProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newUnderlyingQProcess) {
//	underlyingQProcess_ = newUnderlyingQProcess;//reset(newUnderlyingProcess);
//	QL_REQUIRE(underlyingQProcess_ != nullptr, "Underlying process is null!");
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setDefaultProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newDefaultProcess) {
//	defaultProcess_ = newDefaultProcess;//reset(newUnderlyingProcess);
//	QL_REQUIRE(defaultProcess_ != nullptr, "Underlying process is null!");
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setDefaultQProcess(std::shared_ptr<const QuantLib::GeneralizedBlackScholesProcess> newDefaultQProcess) {
//	defaultQProcess_ = newDefaultQProcess;//reset(newUnderlyingQProcess);
//	QL_REQUIRE(defaultQProcess_ != nullptr, "Underlying process is null!");
//
//	return *this;
//}
//
//inline CreditVaRComputationEngine& CreditVaRComputationEngine::setPayoff(std::shared_ptr<const QuantLib::Payoff> newPayoff) {
//	payoff_ = newPayoff;//reset(newUnderlyingProcess);
//	QL_REQUIRE(payoff_ != nullptr, "Underlying process is null!");
//
//	return *this;
//}
//
//double CreditVaRComputationEngine::defaultTime(
//	const std::vector<QuantLib::Time>& cdsMaturity, 
//	double defaultProb,
//	const QuantLib::Handle<QuantLib::BlackVolTermStructure>& vol,
//	double from, 
//	double initialState,
//	double epsilon/* = 1e-6*/) const {
//	double t = 0.0;
//	double tempT1 = 0.000001;
//	double tempT2 = vol->maxTime();
//	double tempProb = 0.0;
//	double gridient = 0.;
//
//	// calculate default time by bisection method
//	if (AT1P(cdsMaturity, vol, from, tempT2, initialState, defaultBarrier_) <= defaultProb)
//		return tempT2;
//	else if (AT1P(cdsMaturity, vol, from, tempT1, initialState, defaultBarrier_) >= defaultProb)
//		return tempT1;
//
//	while (std::abs(
//		AT1P(cdsMaturity, vol, from, (tempT1+tempT2)*0.5, initialState, defaultBarrier_)
//		- defaultProb) > epsilon){
//		if (AT1P(cdsMaturity, vol, from, (tempT1 + tempT2)*0.5, initialState, defaultBarrier_) > defaultProb){
//			tempT2 = (tempT1 + tempT2)*0.5;
//		}
//		else{
//			tempT1 = (tempT1 + tempT2)*0.5;
//		}
//	}
//	return t = (tempT1+tempT2)*0.5;
//}
//
//std::shared_ptr<QuantLib::GeneralizedBlackScholesProcess> CreditVaRComputationEngine::CreateNewUnderlyingProcess(
//	std::string processName,
//	double newState,
//	QuantLib::Date newDate,
//	double interestRate,
//	double dividend,
//	double volatility,
//	double jInt,
//	double pProb,
//	double pJM,
//	double nJM
//	) const {
//	if (processName == typeid(QuantLib::KouProcess).name()){
//		return std::make_shared<QuantLib::KouProcess>(
//			QuantLib::Handle<QuantLib::Quote>(boost::make_shared<QuantLib::SimpleQuote>(newState)),
//			QuantLib::Handle<QuantLib::YieldTermStructure>(boost::make_shared<QuantLib::FlatForward>(newDate, dividend, dayCounter_)),
//			QuantLib::Handle<QuantLib::YieldTermStructure>(boost::make_shared<QuantLib::FlatForward>(newDate, interestRate, dayCounter_)),
//			QuantLib::Handle<QuantLib::BlackVolTermStructure>(boost::make_shared<QuantLib::BlackConstantVol>(newDate, calendar_, volatility, dayCounter_)),
//			jInt,
//			pProb,
//			pJM,
//			nJM
//			);
//	}
//	else{
//		return std::make_shared<QuantLib::BlackScholesProcess>(
//			QuantLib::Handle<QuantLib::Quote>(boost::make_shared<QuantLib::SimpleQuote>(newState)),
//			QuantLib::Handle<QuantLib::YieldTermStructure>(boost::make_shared<QuantLib::FlatForward>(newDate, dividend, dayCounter_)),
//			//QuantLib::Handle<QuantLib::YieldTermStructure>(boost::make_shared<QuantLib::FlatForward>(newDate, interestRate, dayCounter_)),
//			QuantLib::Handle<QuantLib::BlackVolTermStructure>(boost::make_shared<QuantLib::BlackConstantVol>(newDate, calendar_, volatility, dayCounter_))
//			);
//	}
//}