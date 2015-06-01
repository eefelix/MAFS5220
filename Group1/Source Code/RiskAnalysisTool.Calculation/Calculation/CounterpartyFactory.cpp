#include "pch.h"
#include <Calculation/Utilities/DateHelper.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/models/default/cirdefaultmodel.hpp>
#include <ql/models/default/deterministicdefaultmodel.hpp>
#include <vector>

#include "CounterpartyFactory.hpp"
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include "Product.hpp"

using namespace System;
using namespace System::Collections::Generic;
using namespace System::Linq;
using namespace RiskAnalysisTool::Requests;
using namespace RiskAnalysisTool::Instruments;
using namespace RiskAnalysisTool::Calculation;
using namespace RiskAnalysisTool::Calculation::Utilities;

CounterpartyFactory::CounterpartyFactory(boost::shared_ptr<QuantLib::YieldTermStructure> * curve) {
	curve_ = curve;
}

CounterpartyFactory::~CounterpartyFactory() {
	this->!CounterpartyFactory();
}
CounterpartyFactory::!CounterpartyFactory() {
	if (counterparties_)
	{
		delete counterparties_;
		counterparties_ = nullptr;
	}
}


std::vector<boost::shared_ptr<QuantLib::Counterparty>> *CounterpartyFactory::CreateCounterparties(RiskAnalysisTool::Requests::ComputationRequest ^request) {
	this->Visit(request);
	auto r = counterparties_;
	counterparties_ = nullptr;
	return r;
}


void CounterpartyFactory::OnVisit(RiskAnalysisTool::Requests::TvaRequest ^tvaRequest) {
	if (counterparties_) {
		delete counterparties_;
		counterparties_ = nullptr;
	}

	counterparties_ = new std::vector<boost::shared_ptr<QuantLib::Counterparty>>();
	
	boost::shared_ptr<std::vector<double>> issuerCdsSpread(new std::vector<double>), investorCdsSpread(new std::vector<double>);
	boost::shared_ptr<std::vector<QuantLib::Period>> issuerCdsTenor(new std::vector<QuantLib::Period>), investorCdsTenor(new std::vector<QuantLib::Period>);

	for each (auto cds in Enumerable::OfType<CreditDefaultSwap ^>(tvaRequest->MarketData)) {
		String^ symbol = cds->Symbol;
		if (symbol->Contains("ISSUER")) {
			issuerCdsSpread->push_back(cds->Spread);
			issuerCdsTenor->push_back(Utilities::DateHelper::ToQLPeriod(cds->Tenor));
		}
		else {
			investorCdsSpread->push_back(cds->Spread);
			investorCdsTenor->push_back(Utilities::DateHelper::ToQLPeriod(cds->Tenor));
		}
	}

	boost::shared_ptr<QuantLib::DefaultModel> issuerModel
		= boost::make_shared<QuantLib::AT1Pmodel>(issuerCdsTenor, *curve_);

	boost::shared_ptr<QuantLib::DefaultModel> investorModel
		= boost::make_shared<QuantLib::AT1Pmodel>(investorCdsTenor, *curve_);

	boost::shared_ptr<QuantLib::Counterparty> issuer (new QuantLib::Counterparty(
		0, tvaRequest->CounterpartyRecoveryRate, issuerCdsSpread
		, issuerCdsTenor, *curve_
		, issuerModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET(), *curve_));

	issuer->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	boost::shared_ptr<QuantLib::Counterparty> investor = boost::make_shared<QuantLib::Counterparty>(
		0, tvaRequest->InvestorRecoveryRate, investorCdsSpread
		, investorCdsTenor, *curve_
		, investorModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET(), *curve_);
	investor->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	counterparties_->push_back(issuer);
	counterparties_->push_back(investor);
}

void CounterpartyFactory::OnVisit(RiskAnalysisTool::Requests::BvaRequest ^bvaRequest) {
	if (counterparties_) {
		delete counterparties_;
		counterparties_ = nullptr;
	}

	counterparties_ = new std::vector<boost::shared_ptr<QuantLib::Counterparty>>();

	boost::shared_ptr<std::vector<double>> issuerCdsSpread(new std::vector<double>), investorCdsSpread(new std::vector<double>);
	boost::shared_ptr<std::vector<QuantLib::Period>> issuerCdsTenor(new std::vector<QuantLib::Period>), investorCdsTenor(new std::vector<QuantLib::Period>);
	QuantLib::Settings::instance().evaluationDate() = (*curve_)->referenceDate();
	for each (auto cds in Enumerable::OfType<CreditDefaultSwap ^>(bvaRequest->MarketData)) {
		String^ symbol = cds->Symbol;
		if (symbol->Contains("ISSUER")) {
			issuerCdsSpread->push_back(cds->Spread);
			issuerCdsTenor->push_back(Utilities::DateHelper::ToQLPeriod(cds->Tenor));
		}
		else {
			investorCdsSpread->push_back(cds->Spread);
			investorCdsTenor->push_back(Utilities::DateHelper::ToQLPeriod(cds->Tenor));
		}
	}
	
	boost::shared_ptr<QuantLib::DefaultModel> issuerModel
		= boost::make_shared<QuantLib::AT1Pmodel>(issuerCdsTenor, *curve_);

	boost::shared_ptr<QuantLib::DefaultModel> investorModel
		= boost::make_shared<QuantLib::AT1Pmodel>(investorCdsTenor, *curve_);
	
	double issuerRC = bvaRequest->CounterpartyRecoveryRate;
	double investorRC = bvaRequest->InvestorRecoveryRate;
	boost::shared_ptr<QuantLib::Counterparty> issuer = boost::make_shared<QuantLib::Counterparty>(
		0, issuerRC, issuerCdsSpread
		, issuerCdsTenor, *curve_
		, issuerModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), 
		(*curve_)->dayCounter(), 
		QuantLib::TARGET()
		);

	issuer->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	boost::shared_ptr<QuantLib::Counterparty> investor = boost::make_shared<QuantLib::Counterparty>(
		0, investorRC, investorCdsSpread
		, investorCdsTenor, *curve_
		, investorModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET()
		);

	investor->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	counterparties_->push_back(issuer);
	counterparties_->push_back(investor);
}

void CounterpartyFactory::OnVisit(RiskAnalysisTool::Requests::UCvaRequest ^ucvaRequest) {
	if (counterparties_) {
		delete counterparties_;
		counterparties_ = nullptr;
	}

	counterparties_ = new std::vector<boost::shared_ptr<QuantLib::Counterparty>>();

	boost::shared_ptr<std::vector<double>> issuerCdsSpread(new std::vector<double>)/*, investorCdsSpread(new std::vector<double>)*/;
	boost::shared_ptr<std::vector<QuantLib::Period>> issuerCdsTenor(new std::vector<QuantLib::Period>)/*, investorCdsTenor(new std::vector<QuantLib::Period>)*/;

	for each (auto cds in Enumerable::OfType<CreditDefaultSwap ^>(ucvaRequest->MarketData)) {
		issuerCdsSpread->push_back(cds->Spread);
		issuerCdsTenor->push_back(Utilities::DateHelper::ToQLPeriod(cds->Tenor));
	}


	boost::shared_ptr<QuantLib::DefaultModel> issuerModel
		= boost::make_shared<QuantLib::AT1Pmodel>(issuerCdsTenor, *curve_);

	boost::shared_ptr<QuantLib::DefaultModel> investorModel
		= boost::make_shared<QuantLib::AT1Pmodel>(issuerCdsTenor, *curve_, 0., 0.);

	boost::shared_ptr<QuantLib::Counterparty> issuer = boost::make_shared<QuantLib::Counterparty>(
		0, ucvaRequest->RecoveryRate, issuerCdsSpread
		, issuerCdsTenor, *curve_
		, issuerModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET(), *curve_);
	issuer->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	boost::shared_ptr<QuantLib::Counterparty> investor = boost::make_shared<QuantLib::Counterparty>(
		0, 1, issuerCdsSpread
		, issuerCdsTenor, *curve_
		, investorModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET(), *curve_);
	investor->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	counterparties_->push_back(issuer);
	counterparties_->push_back(investor);
}

void CounterpartyFactory::OnVisit(RiskAnalysisTool::Requests::CreditVaRRequest ^cvarRequest) {
	if (counterparties_) {
		delete counterparties_;
		counterparties_ = nullptr;
	}

	counterparties_ = new std::vector<boost::shared_ptr<QuantLib::Counterparty>>();

	boost::shared_ptr<std::vector<double>> issuerCdsSpread(new std::vector<double>)/*, investorCdsSpread(new std::vector<double>)*/;
	boost::shared_ptr<std::vector<QuantLib::Period>> issuerCdsTenor(new std::vector<QuantLib::Period>)/*, investorCdsTenor(new std::vector<QuantLib::Period>)*/;

	for each (auto cds in Enumerable::OfType<CreditDefaultSwap ^>(cvarRequest->MarketData)) {
		issuerCdsSpread->push_back(cds->Spread);
		issuerCdsTenor->push_back(Utilities::DateHelper::ToQLPeriod(cds->Tenor));
	}

	boost::shared_ptr<QuantLib::DefaultModel> issuerModel
		= boost::make_shared<QuantLib::AT1Pmodel>(issuerCdsTenor, *curve_);

	boost::shared_ptr<QuantLib::DefaultModel> investorModel
		= boost::make_shared<QuantLib::AT1Pmodel>(issuerCdsTenor, *curve_, 0., 0.);

	boost::shared_ptr<QuantLib::Counterparty> issuer = boost::make_shared<QuantLib::Counterparty>(
		0, cvarRequest->RecoveryRate, issuerCdsSpread
		, issuerCdsTenor, *curve_
		, issuerModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET(), *curve_);
	issuer->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	boost::shared_ptr<QuantLib::Counterparty> investor = boost::make_shared<QuantLib::Counterparty>(
		0, 1, issuerCdsSpread
		, issuerCdsTenor, *curve_
		, investorModel, QuantLib::Quarterly, QuantLib::Following, QuantLib::DateGeneration::TwentiethIMM
		, (*curve_)->referenceDate(), (*curve_)->dayCounter(), QuantLib::TARGET(), *curve_);
	investor->modelCalibrate(QuantLib::Simplex(0.001), QuantLib::EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

	counterparties_->push_back(issuer);
	counterparties_->push_back(investor);
}

