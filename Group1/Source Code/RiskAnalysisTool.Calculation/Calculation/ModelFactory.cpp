#include "pch.h"
#include <ql/pricingengines/credit/mccrosscurrencyswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mcequityswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mcvanillaswapexposuremodel.hpp>

#include <Calculation/Utilities/DateHelper.hpp>
#include "Product.hpp"
#include "ModelFactory.hpp"
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/processes/cirprocess.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <ql/pricingengines/cirbondengine.hpp>
#include <ql/models/shortrate/calibrationhelpers/zerocouponbondhelper.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/termstructures/yield/flatforward.hpp>

using namespace System;
using namespace System::Linq;
using namespace RiskAnalysisTool::Instruments;
using namespace RiskAnalysisTool::Calculation;
using namespace RiskAnalysisTool::Calculation::Utilities;

ModelFactory::ModelFactory(
	RiskAnalysisTool::Requests::ComputationRequest^ request,
	boost::shared_ptr<QuantLib::YieldTermStructure>* curve,
	boost::shared_ptr<QuantLib::Instrument>* inst) {
	curve_ = curve;
	instrument_ = inst;
	request_ = request;
}

ModelFactory::~ModelFactory() {
	this->!ModelFactory();
}
ModelFactory::!ModelFactory() {
	if (model_)
	{
		delete model_;
		model_ = nullptr;
	}
}


QuantLib::MCExposureModel *ModelFactory::CreateQLModel(RiskAnalysisTool::Instruments::Instrument ^instrument) {
	this->Visit(instrument);
	auto r = model_;
	model_ = nullptr;
	return r;
}


void ModelFactory::OnVisit(RiskAnalysisTool::Instruments::InterestRateSwap ^instrument) {
	if (model_) {
		delete model_;
		model_ = nullptr;
	}
	QuantLib::Date referenceDate
		= Utilities::DateHelper::ToQLDate(request_->EvaluationData);
	model_ = new QuantLib::MCVanillaSwapExposureModel(
		QuantLib::TARGET(), (*curve_)->dayCounter(), referenceDate
		, boost::dynamic_pointer_cast<QuantLib::VanillaSwap>(*instrument_));
}

void ModelFactory::OnVisit(RiskAnalysisTool::Instruments::EquitySwap ^instrument) {
	if (model_) {
		delete model_;
		model_ = nullptr;
	}
	QuantLib::Date referenceDate
		= Utilities::DateHelper::ToQLDate(request_->EvaluationData);

	boost::shared_ptr<QuantLib::StochasticProcess1D> instrumentProcess;

	String^ stockSym = instrument->PaySymbol;
	if (stockSym->Contains(request_->DomesticCurrency))
		stockSym = instrument->ReceiveSymbol;

	std::string stockStr = Utilities::StringHelper::ToNativeString(stockSym);
	auto pos = Product::symToProcess.find(stockStr);

	if (pos != Product::symToProcess.end()){
		instrumentProcess = pos->second;
	}
	else{
		instrumentProcess =
			boost::make_shared<QuantLib::BlackScholesMertonProcess>(
			QuantLib::Handle<QuantLib::Quote>(boost::make_shared<QuantLib::SimpleQuote>(instrument->SpotPrice))
			, QuantLib::Handle<QuantLib::YieldTermStructure>(boost::make_shared<QuantLib::FlatForward>(
			referenceDate, instrument->DividendYield, (*curve_)->dayCounter()))
			, QuantLib::Handle<QuantLib::YieldTermStructure>(*curve_)
			, QuantLib::Handle<QuantLib::BlackVolTermStructure>(
				boost::make_shared<QuantLib::BlackConstantVol>(referenceDate, (*curve_)->calendar(), instrument->Volatility, (*curve_)->dayCounter())));
		Product::symToProcess.emplace(std::make_pair(stockStr, instrumentProcess));
	}

	model_ = new QuantLib::MCEquitySwapExposureModel(
		QuantLib::TARGET(), (*curve_)->dayCounter(), referenceDate
		, instrumentProcess, boost::dynamic_pointer_cast<QuantLib::EquitySwap>(*instrument_));
}

void ModelFactory::OnVisit(RiskAnalysisTool::Instruments::CrossCurrencySwap ^instrument) {
	if (model_) {
		delete model_;
		model_ = nullptr;
	}

	//boost::shared_ptr<QuantLib::CrossCurrencySwap> ccsPtr(dynamic_cast<QuantLib::CrossCurrencySwap*>(instrument_));
	QuantLib::Date referenceDate
		= Utilities::DateHelper::ToQLDate(request_->EvaluationData);

	String ^paySymbol = instrument->PaySymbol,
		^receiveSymbol = instrument->ReceiveSymbol;
	std::vector<QuantLib::Date> fBondDate;
	std::vector<double> fdiscountor;

	for each (auto bond in Enumerable::OfType<Bond ^>(request_->MarketData)) {
		String^ symbol = bond->Symbol;
		if ((!symbol->Contains(request_->DomesticCurrency)
			&& (symbol->Contains(paySymbol) || symbol->Contains(receiveSymbol))))
		{
			fBondDate.push_back(Utilities::DateHelper::ToQLDate(bond->MaturityDate));
			fdiscountor.push_back(bond->Price / bond->Principle);
		}
	}

	if (std::find(fBondDate.begin(), fBondDate.end(), referenceDate) == fBondDate.end()){
		fBondDate.push_back(referenceDate);
		fdiscountor.push_back(1.);
	}

	//assert(foreignBondDate.size() == foreignDiscountor.size());
	auto fCurve
		= boost::make_shared<QuantLib::DiscountCurve>(fBondDate, fdiscountor, QuantLib::Actual360(), QuantLib::TARGET());
	auto fCIR = boost::make_shared<QuantLib::CoxIngersollRoss>();

	std::vector<boost::shared_ptr<QuantLib::CalibrationHelper>> zcbonds;
	boost::shared_ptr<QuantLib::PricingEngine> engine(
		new QuantLib::CIRBondEngine(fCIR, QuantLib::Handle<QuantLib::YieldTermStructure>(fCurve)));
	for (int i = 0; i < fBondDate.size() - 1; ++i){
		if (fBondDate[i] != referenceDate){
			zcbonds.push_back(boost::shared_ptr<QuantLib::ZerocouponbondHelper>(new QuantLib::ZerocouponbondHelper(
				0,
				fCurve->calendar(),
				1,
				fBondDate[i],
				QuantLib::Following,
				100.0,
				referenceDate,
				QuantLib::Handle<QuantLib::Quote>(boost::make_shared<QuantLib::SimpleQuote>(0.1)),
				QuantLib::Handle<QuantLib::YieldTermStructure>(fCurve))));
			zcbonds.back()->setPricingEngine(engine);
		}
	}
	QuantLib::Simplex solver(0.001);
	const QuantLib::Size maxIteration = 1000;
	const QuantLib::Size minStatIteration = 50;
	const QuantLib::Real rootEpsilon = 1e-8;
	const QuantLib::Real FunctionEpsilon = 1e-8;
	const QuantLib::Real gradientNormEpsilon = 1e-8;
	const QuantLib::EndCriteria endcriteria = QuantLib::EndCriteria(maxIteration, minStatIteration, rootEpsilon, FunctionEpsilon, gradientNormEpsilon);

	fCIR->calibrate(zcbonds, solver, endcriteria, *(fCIR->constraint()));

	boost::shared_ptr<QuantLib::StochasticProcess1D> foreignRateProcess;

	String^ foreignSym = instrument->PaySymbol;
	if (foreignSym->Contains(request_->DomesticCurrency))
		foreignSym = instrument->ReceiveSymbol;

	std::string foreignStr = Utilities::StringHelper::ToNativeString(foreignSym);
	auto pos = Product::symToProcess.find(foreignStr);

	if (pos != Product::symToProcess.end()){
		foreignRateProcess = pos->second;
	}
	else {
		foreignRateProcess
			= boost::make_shared<QuantLib::CIRprocess>(
			fCIR->params()[0], fCIR->params()[1], fCIR->params()[2], fCIR->params()[3]);
		Product::symToProcess.emplace(std::make_pair(foreignStr, foreignRateProcess));
	}

	model_ = new QuantLib::MCCrossCurrencySwapExposureModel(
		QuantLib::TARGET(), (*curve_)->dayCounter(), referenceDate
		, foreignRateProcess, boost::dynamic_pointer_cast<QuantLib::CrossCurrencySwap>(*instrument_));
}
