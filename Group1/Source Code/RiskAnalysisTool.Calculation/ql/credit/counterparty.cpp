#include <pch.h>
#include <ql/credit/counterparty.hpp>
#include <ql/models/default/deterministicdefaultmodel.hpp>
#include <ql/errors.hpp>
#include <boost/cast.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/models/defaultcdshelper.hpp>

//! \file counterparty.cpp
using namespace QuantLib;

Counterparty::Counterparty(
	int settlementDays,
	double recoveryRate,
	const boost::shared_ptr<const std::vector<double>>& cdsSpreads,
	const boost::shared_ptr<const std::vector<Period>>& cdsTenors,
	const boost::shared_ptr<YieldTermStructure>& discountCurve,
	const boost::shared_ptr<DefaultModel>& model,
	const QuantLib::Frequency& freq /*= QuantLib::Quarterly*/,
	const QuantLib::BusinessDayConvention& busDayConvention /*= QuantLib::Following*/,
	const QuantLib::DateGeneration::Rule& dateGenerationRule /*= QuantLib::DateGeneration::TwentiethIMM*/,
	const QuantLib::Date& referenceDate/* = QuantLib::Date()*/,
	const QuantLib::DayCounter& dayCounter /*= QuantLib::Actual365Fixed()*/,
	const QuantLib::Calendar& calendar/* = QuantLib::TARGET()*/,
	const boost::shared_ptr<YieldTermStructure>& defaultdiscountCurve/* for fva */)
	: model_(model), referenceDate_(referenceDate), cdsTenors_(cdsTenors)
	, recoveryRate_(recoveryRate), dayCounter_(dayCounter), discountCurve_(discountCurve), defaultdiscountCurve_(defaultdiscountCurve),
	CALIBRATION_FLAG_(false), FUNDINGSPREAD_FLAG_(false) {
	bootStrap(settlementDays, cdsSpreads, freq, busDayConvention, dateGenerationRule, calendar);
	cdsHelperGenerator();
}

Counterparty::~Counterparty() {
}

void Counterparty::bootStrap(
	int settlementDays,
	const boost::shared_ptr<const std::vector<double>>& cdsSpreads,
	const QuantLib::Frequency& freq /*= QuantLib::Quarterly*/,
	const QuantLib::BusinessDayConvention& busDayConvention /*= QuantLib::Following*/,
	const QuantLib::DateGeneration::Rule& dateGenerationRule /*= QuantLib::DateGeneration::TwentiethIMM*/,
	const QuantLib::Calendar& calendar /*= QuantLib::TARGET()*/) {
	//QuantLib::Settings::instance().evaluationDate() = referenceDate_;

	std::vector<boost::shared_ptr<QuantLib::DefaultProbabilityHelper> > instruments;
	for (int i = 0; i < cdsTenors_->size(); i++) {
		instruments.push_back(
			boost::make_shared<QuantLib::SpreadCdsHelper>(
			QuantLib::Handle<QuantLib::Quote>(
			boost::make_shared<QuantLib::SimpleQuote>(cdsSpreads->at(i))),
			cdsTenors_->at(i),
			settlementDays,
			calendar,
			freq,
			busDayConvention,
			dateGenerationRule,
			dayCounter_,
			recoveryRate_,
			QuantLib::Handle<QuantLib::YieldTermStructure>(discountCurve_)));
	}

	hazardRateStructure_ =
		boost::make_shared<QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat>>(
		referenceDate_,
		instruments,
		dayCounter_);
}

void Counterparty::cdsHelperGenerator() {
	std::vector<Time> time(cdsTenors_->size());

	for (auto itr = cdsTenors_->begin(); itr != cdsTenors_->end(); ++itr) {
		time[itr - cdsTenors_->begin()] = dayCounter_.yearFraction(referenceDate_, referenceDate_ + *itr);
	}

	cdsHelper_.reserve(time.size());

	boost::shared_ptr<CalibratedModel> calibModel
		= boost::dynamic_pointer_cast<CalibratedModel>(model_);
	if (calibModel) {
		for (auto itr = time.begin(); itr != time.end(); ++itr) {
			cdsHelper_.push_back(
				boost::make_shared<DefaultCdsHelper>(*itr,
				Handle<Quote>(boost::make_shared<SimpleQuote>(hazardRateStructure_->defaultProbability(*itr, true))),
				Handle<YieldTermStructure>(discountCurve_), calibModel));
		}
	}
}

void Counterparty::modelCalibrate(
	OptimizationMethod &method, const EndCriteria &endCriteria, bool needCalibrate /*= true*/,
	const Constraint &constraint/* = Constraint()*/,
	const std::vector<Real> &weights /*= std::vector<Real>()*/,
	const std::vector<bool> &fixParameters /*= std::vector<bool>()*/) {
	boost::shared_ptr<DeterministicDefaultModel> determinModel
		= boost::dynamic_pointer_cast<DeterministicDefaultModel> (model_);
	if (needCalibrate) {
		if (!determinModel.get()) {
			boost::shared_ptr<CalibratedModel> calibModel
				= boost::dynamic_pointer_cast<CalibratedModel> (model_);
			if (calibModel.get()) {
				calibModel->calibrate(cdsHelper_, method, endCriteria, constraint, weights, fixParameters);
			} else {
				QL_FAIL("This model doesn't provide calibration method.");
			}
		} else {
			determinModel->setHazardRateStructure(hazardRateStructure_);
		}
	}
	CALIBRATION_FLAG_ = true;
}

const Real Counterparty::getDefaultProb(const Time& t) const {
	QL_REQUIRE(CALIBRATION_FLAG_,
		"please calibrate model first.");
	return model_->defaultProbability(t);
}

const Time Counterparty::getDefaultTime(const Path & path) const {
	QL_REQUIRE(CALIBRATION_FLAG_,
		"please calibrate model first");
	return model_->defaultTime(path);
}

const Time Counterparty::defaultTimeGenerator(const Probability& prob, const double& epsilon) const {
	double t1 = 1e-6;
	// the last date of CDS
	double t2 = dayCounter_.yearFraction(referenceDate_, referenceDate_ + cdsTenors_->back());

	if (std::abs(getDefaultProb(t1) - prob) < epsilon)
		return t1;

	if (std::abs(getDefaultProb(t2) - prob) < epsilon)
		return t2;

	if (getDefaultProb(t1) > prob) {
		t1 /= 2.0;
		while (getDefaultProb(t1) > prob) {
			if (std::abs(getDefaultProb(t1) - prob) < epsilon)
				return t1;
			t2 = t1;
			t1 /= 2.0;
		}
	} else if (getDefaultProb(t2) < prob) {
		t2 *= 2.0;
		while (getDefaultProb(t2) < prob) {
			if (std::abs(getDefaultProb(t2) - prob) < epsilon)
				return t2;
			t1 = t2;
			t2 *= 2.0;
		}
	}

	while (std::abs(getDefaultProb(0.5*(t1 + t2)) - prob) > epsilon) {
		if (getDefaultProb(0.5*(t1 + t2)) > prob)
			t2 = 0.5*(t1 + t2);
		else if (getDefaultProb(0.5*(t1 + t2)) < prob)
			t1 = 0.5*(t1 + t2);
	}

	return 0.5*(t1 + t2);
}


const boost::shared_ptr<const DefaultModel> Counterparty::getModel() const {
	return model_;
}

double Counterparty::fundingSpread(Time t) const {
	if (defaultdiscountCurve_ == nullptr) {
		QL_FAIL("Cannot calculate funding spread!!! Missing default discount curve!!!");
	}

	std::size_t n = static_cast<std::size_t>(ceil(dayCounter_.yearFraction(referenceDate_, referenceDate_ + (cdsTenors_->back())) * 360) + 1);

	if (FUNDINGSPREAD_FLAG_ == false) {
		fundingSpread_.assign(n, 0);
		for (std::size_t i = 1; i <= n - 1; ++i) {
			fundingSpread_[i] = defaultdiscountCurve_->forwardRate(1. / 360 * i, 1. / 360 * i + 0.0001, Continuous) -
				discountCurve_->forwardRate(1. / 360 * i, 1. / 360 * i + 0.0001, Continuous) -
				(1 - recoveryRate_)*hazardRateStructure_->hazardRate(1. / 360 * i);
		}
		FUNDINGSPREAD_FLAG_ = true;
	}
	if (floor(t * 360 + 0.5) < fundingSpread_.size())
		return fundingSpread_[static_cast<std::size_t>(std::floor(t * 360 + 0.5))];
	else return fundingSpread_.back();
}