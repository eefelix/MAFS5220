#include "pch.h"
#include <ql/instruments/credit/tva.hpp>
#include <ql/instruments/credit/bva.hpp>
#include <ql/instruments/credit/var.hpp>

#include "PortfolioManagerFactory.hpp"

using namespace RiskAnalysisTool::Requests;
using namespace RiskAnalysisTool::Calculation;
using namespace RiskAnalysisTool::Calculation::Utilities;

PortfolioManagerFactory::PortfolioManagerFactory() {
}

PortfolioManagerFactory::~PortfolioManagerFactory() {
	this->!PortfolioManagerFactory();
}
PortfolioManagerFactory::!PortfolioManagerFactory() {
	if (portfolioManager_)
	{
		delete portfolioManager_;
		portfolioManager_ = nullptr;
	}
}


QuantLib::Instrument *PortfolioManagerFactory::CreatePortfolioManager(RiskAnalysisTool::Requests::ComputationRequest ^request) {
	this->Visit(request);
	auto r = portfolioManager_;
	portfolioManager_ = nullptr;
	return r;
}


void PortfolioManagerFactory::OnVisit(RiskAnalysisTool::Requests::TvaRequest ^tvaRequest) {
	if (portfolioManager_) {
		delete portfolioManager_;
		portfolioManager_ = nullptr;
	}

	portfolioManager_ = new QuantLib::TVA;
}

void PortfolioManagerFactory::OnVisit(RiskAnalysisTool::Requests::BvaRequest ^tvaRequest) {
	if (portfolioManager_) {
		delete portfolioManager_;
		portfolioManager_ = nullptr;
	}

	portfolioManager_ = new QuantLib::BVA;
}

void PortfolioManagerFactory::OnVisit(RiskAnalysisTool::Requests::UCvaRequest ^tvaRequest) {
	if (portfolioManager_) {
		delete portfolioManager_;
		portfolioManager_ = nullptr;
	}

	portfolioManager_ = new QuantLib::BVA;
}

void PortfolioManagerFactory::OnVisit(RiskAnalysisTool::Requests::CreditVaRRequest ^tvaRequest) {
	if (portfolioManager_) {
		delete portfolioManager_;
		portfolioManager_ = nullptr;
	}

	portfolioManager_ = new QuantLib::VAR;
}

