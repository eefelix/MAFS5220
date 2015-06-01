#include <pch.h>
#include <ql/calibration/kouprocesscalibrator.hpp>
#include <ql/pricingengines/vanilla/analytickoueuropeanengine.hpp>

//! \file kouprocesscalibrator.cpp
//!
using namespace QuantLib;

KouProcessCalibrator::KouProcessCalibrator(
	const QuantLib::Date& referenceDate, const QuantLib::Calendar& calendar, const QuantLib::DayCounter& dayCounter,
	double riskFreeRate, double spotPrice, double dividend,
	std::shared_ptr<std::vector<double>> strike,
	std::shared_ptr<std::vector<QuantLib::Date>> maturityDate,
	std::shared_ptr<std::vector<double>> optionPrice,
	std::shared_ptr<QuantLib::EndCriteria> endcriteria)
	:referenceDate_(referenceDate), calendar_(calendar), dayCounter_(dayCounter),
	riskFreeRate_(riskFreeRate), spotPrice_(spotPrice), dividend_(dividend),
	strike_(strike), maturityDate_(maturityDate), optionPrice_(optionPrice),
	endcriteria_(endcriteria){

	volImplied_ = 0.1; //volatility
	posJumpMeanImplied_ = 1.0; // positive jump mean
	negJumpMeanImplied_ = 1.0; // negative jump mean
	posProbability_ = 0.5; // positive jump probability
	jumpIntensity_ = 2.0; //jump intensity
};

void KouProcessCalibrator::calibrate(){
	QuantLib::Array initialValue(5);
	initialValue[0] = 0.1; //volatility
	initialValue[1] = 2.0; // positive jump mean
	initialValue[2] = 1.0; // negative jump mean
	initialValue[3] = 0.5; // positive jump probability
	initialValue[4] = 2.0; //jump intensity

	QuantLib::Array lowBoundaryConstraints(5), highBoundaryConstraints(5);
	lowBoundaryConstraints[0] = 0.01;
	lowBoundaryConstraints[1] = 1.001;
	lowBoundaryConstraints[2] = 0.001;
	lowBoundaryConstraints[3] = 0.0001;
	lowBoundaryConstraints[4] = 0.001;

	highBoundaryConstraints[0] = 5.0;
	highBoundaryConstraints[1] = 100.0;
	highBoundaryConstraints[2] = 100.0;
	highBoundaryConstraints[3] = 0.9999;
	highBoundaryConstraints[4] = 100.0;

	QuantLib::NonhomogeneousBoundaryConstraint boundaryConstraint(lowBoundaryConstraints, highBoundaryConstraints);
	QuantLib::Problem modelCalibrate(*this, boundaryConstraint, initialValue);
	QuantLib::Simplex solver(0.05);
	QuantLib::EndCriteria::Type solvedCriteria = solver.minimize(modelCalibrate, *endcriteria_);

	QuantLib::Array currentValue = modelCalibrate.currentValue();

	volImplied_ = modelCalibrate.currentValue().at(0);
	posJumpMeanImplied_ = modelCalibrate.currentValue().at(1);
	negJumpMeanImplied_ = modelCalibrate.currentValue().at(2);
	posProbability_ = modelCalibrate.currentValue().at(3);
	jumpIntensity_ = modelCalibrate.currentValue().at(4);

}


QuantLib::Real KouProcessCalibrator::value(const QuantLib::Array& x) const {
	// calculate square error
	QuantLib::Array tempRes = values(x);
	QuantLib::Real residual = 0.;
	for (auto i = tempRes.begin(); i != tempRes.end(); ++i)
		residual += (*i)*(*i);

	return residual;
}

QuantLib::Disposable<QuantLib::Array> KouProcessCalibrator::values(const QuantLib::Array& x) const {
	double sigma = x[0];
	double eta1 = x[1];
	double eta2 = x[2];
	double p = x[3];
	double lambda = x[4];

	double r = riskFreeRate_ - dividend_;

	// copy from file analytickoueuropeansengine.cpp
	double zeta = p*eta1 / (eta1 - 1) + (1 - p)*eta2 / (eta2 + 1) - 1;
	double adjustEta1 = eta1 - 1;
	double adjustEta2 = eta2 + 1;
	double adjustLambda = lambda*(zeta + 1);
	double adjustP = p / (1 + zeta)*eta1 / (eta1 - 1);

	//
	std::vector<double>	maturity(maturityDate_->size());

	for (int i = 0; i <maturity.size(); ++i){
		maturity[i] = dayCounter_.yearFraction(referenceDate_, maturityDate_->at(i));
	}
	// calculate price difference
	QuantLib::Array residual(optionPrice_->size());
	double tolerance_ = 1e-4;
	for (int i = 0; i < residual.size(); ++i)
		residual[i]
		= spotPrice_*std::exp(-dividend_*maturity[i])
		* KouHelper::Gamma(r + 0.5*sigma*sigma - lambda*zeta, sigma, adjustLambda, adjustP, adjustEta1, adjustEta2, std::log(strike_->at(i) / spotPrice_), maturity[i], tolerance_)
		- strike_->at(i)*exp(-riskFreeRate_* maturity[i])
		* KouHelper::Gamma(r - 0.5*sigma*sigma - lambda*zeta, sigma, lambda, p, eta1, eta2, log(strike_->at(i) / spotPrice_), maturity[i], tolerance_)
		- optionPrice_->at(i);

	return residual;
}

inline QuantLib::BlackConstantVol KouProcessCalibrator::getVol() {
	return QuantLib::BlackConstantVol(referenceDate_, calendar_, volImplied_, dayCounter_);
}

inline double KouProcessCalibrator::getPosJumpMean() {
	return posJumpMeanImplied_;
}

inline double KouProcessCalibrator::getNegJumpMean() {
	return negJumpMeanImplied_;
}

inline double KouProcessCalibrator::getPosProb() {
	return posProbability_;
}

inline double KouProcessCalibrator::getJumpIntensity() {
	return jumpIntensity_;
}

