#include "pch.h"
#include <Calculation/Utilities/DateHelper.hpp>
#include <ql/pricingengines/credit/mctvaengine.hpp>
#include <ql/pricingengines/credit/mcbvaengine.hpp>
#include <ql/pricingengines/credit/mccreditvarengine.hpp>
#include "Product.hpp"
#include "EngineFactory.hpp"
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>

using namespace System;
using namespace System::Linq;
using namespace RiskAnalysisTool::Requests;
using namespace RiskAnalysisTool::Instruments;
using namespace RiskAnalysisTool::Calculation;
using namespace RiskAnalysisTool::Calculation::Utilities;

EngineFactory::EngineFactory(boost::shared_ptr<QuantLib::YieldTermStructure> * curve, std::vector<boost::shared_ptr<const QuantLib::Counterparty>>* counterparties,
	std::vector < boost::shared_ptr<const QuantLib::MCExposureModel>>* models, boost::shared_ptr<QuantLib::CoxIngersollRoss>* domesticModel) {
	curve_ = curve;
	counterparties_ = counterparties;
	models_ = models;
	domesticModel_ = domesticModel;
}

EngineFactory::~EngineFactory() {
	this->!EngineFactory();
}
EngineFactory::!EngineFactory() {
	/*if (engine_)
	{
	delete engine_;
	engine_ = nullptr;
	}*/
}


boost::shared_ptr<QuantLib::PricingEngine> EngineFactory::CreateEngine(RiskAnalysisTool::Requests::ComputationRequest ^request) {
	this->Visit(request);
	auto r = engine_;
	engine_.reset();
	return r.native();
}


void EngineFactory::OnVisit(RiskAnalysisTool::Requests::TvaRequest ^tvaRequest) {
	/*if (engine_) {
		delete engine_;
		engine_ = nullptr;
		}*/
	engine_.reset();

	double endTime = Utilities::EndTimeHelper::GetEndTime((*curve_)->dayCounter(), tvaRequest);
	std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> processes;
	processes.reserve(Product::symToProcess.size());
	std::vector<std::string> label(
	{ "ISSUER", "INVESTOR", Utilities::StringHelper::ToNativeString(tvaRequest->DomesticCurrency) });
	processes.push_back(Product::symToProcess[label[0]]);
	processes.push_back(Product::symToProcess[label[1]]);
	processes.push_back(Product::symToProcess[label[2]]);
	for (auto itr = Product::symToProcess.begin(); itr != Product::symToProcess.end(); ++itr) {
		if (std::find(label.begin(), label.end(), itr->first) == label.end()){
			label.push_back(itr->first);
			processes.push_back(itr->second);
		}
	}

	QuantLib::Matrix correlationMatrix = QuantLib::Matrix(Product::symToProcess.size(), Product::symToProcess.size(), 0.);
	for (QuantLib::Size i = 0; i < correlationMatrix.rows(); ++i)
		correlationMatrix[i][i] = 1.;

	for each (auto corr in Enumerable::OfType<Correlation ^>(tvaRequest->MarketData)) {
		std::string str1 = Utilities::StringHelper::ToNativeString(corr->Symbol1),
			str2 = Utilities::StringHelper::ToNativeString(corr->Symbol2);
		QuantLib::Size i = 0, j = 0;
		auto pos = std::find(label.begin(), label.end(), str1);
		i = std::distance(label.begin(), pos);
		pos = std::find(label.begin(), label.end(), str2);
		j = std::distance(label.begin(), pos);

		correlationMatrix[i][j] = corr->Value;
		correlationMatrix[j][i] = corr->Value;
	}

	engine_ = boost::make_shared<QuantLib::MCTVAEngine<>>(
		*counterparties_, *models_, correlationMatrix, processes
		, QuantLib::Handle<QuantLib::YieldTermStructure>(*curve_), *domesticModel_
		, QuantLib::Null<QuantLib::Size>(), 360, false, 300, QuantLib::Null<QuantLib::Real>(), 50000, endTime);
}

void EngineFactory::OnVisit(RiskAnalysisTool::Requests::BvaRequest ^bvaRequest) {
	/*if (engine_) {
		delete engine_;
		engine_ = nullptr;
		}*/
	engine_.reset();

	double endTime = Utilities::EndTimeHelper::GetEndTime((*curve_)->dayCounter(), bvaRequest);
	std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> processes;
	processes.reserve(Product::symToProcess.size());
	std::vector<std::string> label(
	{ "ISSUER", "INVESTOR", Utilities::StringHelper::ToNativeString(bvaRequest->DomesticCurrency) });
	processes.push_back(Product::symToProcess[label[0]]);
	processes.push_back(Product::symToProcess[label[1]]);
	processes.push_back(Product::symToProcess[label[2]]);
	for (auto itr = Product::symToProcess.begin(); itr != Product::symToProcess.end(); ++itr) {
		if (std::find(label.begin(), label.end(), itr->first) == label.end()){
			label.push_back(itr->first);
			processes.push_back(itr->second);
		}
	}

	QuantLib::Matrix correlationMatrix = QuantLib::Matrix(Product::symToProcess.size(), Product::symToProcess.size(), 0.);
	for (QuantLib::Size i = 0; i < correlationMatrix.rows(); ++i)
		correlationMatrix[i][i] = 1.;

	for each (auto corr in Enumerable::OfType<Correlation ^>(bvaRequest->MarketData)) {
		std::string str1 = Utilities::StringHelper::ToNativeString(corr->Symbol1),
			str2 = Utilities::StringHelper::ToNativeString(corr->Symbol2);
		QuantLib::Size i = 0, j = 0;
		auto pos = std::find(label.begin(), label.end(), str1);
		i = std::distance(label.begin(), pos);
		pos = std::find(label.begin(), label.end(), str2);
		j = std::distance(label.begin(), pos);

		correlationMatrix[i][j] = corr->Value;
		correlationMatrix[j][i] = corr->Value;
	}
	assert(processes[0] == Product::symToProcess["ISSUER"]);
	assert(processes[1] == Product::symToProcess["INVESTOR"]);
	assert(processes[2] == Product::symToProcess["$USD"]);

	engine_ = boost::make_shared<QuantLib::MCBVAEngine<>>(
		*counterparties_, *models_, correlationMatrix, processes, QuantLib::Handle<QuantLib::YieldTermStructure>(*curve_)
		, *domesticModel_, QuantLib::Null<QuantLib::Size>()
		, 360, false, 10000, QuantLib::Null<QuantLib::Real>(), 50000, endTime);
}

void EngineFactory::OnVisit(RiskAnalysisTool::Requests::UCvaRequest ^ucvaRequest) {
	/*if (engine_) {
		delete engine_;
		engine_ = nullptr;
		}*/
	engine_.reset();
	double endTime = Utilities::EndTimeHelper::GetEndTime((*curve_)->dayCounter(), ucvaRequest);
	std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> processes;
	processes.reserve(Product::symToProcess.size());
	std::vector<std::string> label(
	{ "ISSUER", "INVESTOR", Utilities::StringHelper::ToNativeString(ucvaRequest->DomesticCurrency) });
	processes.push_back(Product::symToProcess[label[0]]);
	processes.push_back(Product::symToProcess[label[1]]);
	processes.push_back(Product::symToProcess[label[2]]);
	for (auto itr = Product::symToProcess.begin(); itr != Product::symToProcess.end(); ++itr) {
		if (std::find(label.begin(), label.end(), itr->first) == label.end() && itr->first!="INVESTOR"){
			label.push_back(itr->first);
			processes.push_back(itr->second);
		}
	}

	QuantLib::Matrix correlationMatrix = QuantLib::Matrix(Product::symToProcess.size(), Product::symToProcess.size(), 0.);
	for (std::size_t i = 0; i < correlationMatrix.rows(); ++i)
		correlationMatrix[i][i] = 1.;

	for each (auto corr in Enumerable::OfType<Correlation ^>(ucvaRequest->MarketData)) {
		std::string str1 = Utilities::StringHelper::ToNativeString(corr->Symbol1),
			str2 = Utilities::StringHelper::ToNativeString(corr->Symbol2);
		std::size_t i = 0, j = 0;
		auto pos = std::find(label.begin(), label.end(), str1);
		i = std::distance(label.begin(), pos);
		pos = std::find(label.begin(), label.end(), str2);
		j = std::distance(label.begin(), pos);

		correlationMatrix[i][j] = corr->Value;
		correlationMatrix[j][i] = corr->Value;
	}

	engine_ = boost::make_shared<QuantLib::MCBVAEngine<>>(
		*counterparties_, *models_, correlationMatrix, processes, QuantLib::Handle<QuantLib::YieldTermStructure>(*curve_)
		, *domesticModel_, QuantLib::Null<QuantLib::Size>()
		, 360, false, 50, QuantLib::Null<QuantLib::Real>(), 50000, endTime);
}

void EngineFactory::OnVisit(RiskAnalysisTool::Requests::CreditVaRRequest ^cvarRequest) {
	/*if (engine_) {
		delete engine_;
		engine_ = nullptr;
		}*/
	engine_.reset();
	double endTime = Utilities::EndTimeHelper::GetEndTime((*curve_)->dayCounter(), cvarRequest);
	std::vector<boost::shared_ptr<QuantLib::StochasticProcess1D>> processes;
	processes.reserve(Product::symToProcess.size());
	std::vector<std::string> label(
	{ "ISSUER", /*"INVESTOR",*/ Utilities::StringHelper::ToNativeString(cvarRequest->DomesticCurrency) });
	processes.push_back(Product::symToProcess[label[0]]);
	processes.push_back(Product::symToProcess[label[1]]);
	//processes.push_back(Product::symToProcess[label[2]]);
	for (auto itr = Product::symToProcess.begin(); itr != Product::symToProcess.end(); ++itr) {
		if (std::find(label.begin(), label.end(), itr->first) == label.end() && itr->first != "INVESTOR"){
			label.push_back(itr->first);
			processes.push_back(itr->second);
		}
	}

	QuantLib::Matrix correlationMatrix = QuantLib::Matrix(Product::symToProcess.size()-1, Product::symToProcess.size()-1, 0.);
	for (QuantLib::Size i = 0; i < correlationMatrix.rows(); ++i)
		correlationMatrix[i][i] = 1.;

	for each (auto corr in Enumerable::OfType<Correlation ^>(cvarRequest->MarketData)) {
		std::string str1 = Utilities::StringHelper::ToNativeString(corr->Symbol1),
			str2 = Utilities::StringHelper::ToNativeString(corr->Symbol2);
		QuantLib::Size i = 0, j = 0;
		auto pos = std::find(label.begin(), label.end(), str1);
		i = std::distance(label.begin(), pos);
		pos = std::find(label.begin(), label.end(), str2);
		j = std::distance(label.begin(), pos);

		correlationMatrix[i][j] = corr->Value;
		correlationMatrix[j][i] = corr->Value;
	}

	engine_ = boost::make_shared<QuantLib::MCCreditVaREngine<>>(
		*counterparties_, *models_, correlationMatrix, processes, QuantLib::Handle<QuantLib::YieldTermStructure>(*curve_)
		, *domesticModel_, QuantLib::Null<QuantLib::Size>()
		, 360, false, 50, QuantLib::Null<QuantLib::Real>(), 50000, 1.);
}

