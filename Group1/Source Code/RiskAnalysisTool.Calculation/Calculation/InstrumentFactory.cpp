#include "pch.h"

#include <ql/instruments/vanillaswap.hpp>
#include <ql/instruments/swap/equityswap.hpp>
#include <ql/instruments/swap/crosscurrencyswap.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/time/daycounters/thirty360.hpp>
#include <ql/currencies/asia.hpp>

#include <Calculation/Utilities/DateHelper.hpp>

#include "InstrumentFactory.hpp"


using namespace System;
using namespace RiskAnalysisTool::Instruments;
using namespace RiskAnalysisTool::Calculation;
using namespace RiskAnalysisTool::Calculation::Utilities;
using namespace RiskAnalysisTool::Requests;

InstrumentFactory::InstrumentFactory(ComputationRequest^ request, boost::shared_ptr<QuantLib::YieldTermStructure> *curve) {
	curve_ = curve;
	request_ = request;
}

InstrumentFactory::~InstrumentFactory() {
	this->!InstrumentFactory();
}
InstrumentFactory::!InstrumentFactory() {
	if (instrument_)
	{
		delete instrument_;
		instrument_ = nullptr;
	}
}


QuantLib::Instrument *InstrumentFactory::CreateQLInstrument(RiskAnalysisTool::Instruments::Instrument ^instrument) {
	this->Visit(instrument);
	auto r = instrument_;
	instrument_ = nullptr;
	return r;
}


void InstrumentFactory::OnVisit(RiskAnalysisTool::Instruments::InterestRateSwap ^instrument) {
	if (instrument_) {
		delete instrument_;
		instrument_ = nullptr;
	}

	QuantLib::Date referenceDate
		= Utilities::DateHelper::ToQLDate(request_->EvaluationData);

	QuantLib::VanillaSwap::Type irsType;
	QuantLib::Schedule fixedSchedule;
	QuantLib::Schedule floatSchedule;

	QuantLib::Date startDate = Utilities::DateHelper::ToQLDate(instrument->StartDate);
	QuantLib::Date maturityDate = Utilities::DateHelper::ToQLDate(instrument->MaturityDate);

	if (instrument->IsPayFloating){
		irsType = QuantLib::VanillaSwap::Type::Receiver;
		fixedSchedule = QuantLib::Schedule(startDate, maturityDate,
			instrument->ReceiveLegFrequency->Value * QuantLib::TimeUnit((int)instrument->ReceiveLegFrequency->Unit + 1),
			QuantLib::TARGET(), QuantLib::Following, QuantLib::Following, QuantLib::DateGeneration::Forward, false);
		floatSchedule = QuantLib::Schedule(startDate, maturityDate,
			instrument->PayLegFrequency->Value * QuantLib::TimeUnit((int)instrument->PayLegFrequency->Unit + 1),
			QuantLib::TARGET(), QuantLib::Following, QuantLib::Following, QuantLib::DateGeneration::Forward, false);
	}
	else{
		irsType = QuantLib::VanillaSwap::Type::Payer;
		floatSchedule = QuantLib::Schedule(startDate, maturityDate,
			instrument->ReceiveLegFrequency->Value * QuantLib::TimeUnit((int)instrument->ReceiveLegFrequency->Unit + 1),
			QuantLib::TARGET(), QuantLib::Following, QuantLib::Following, QuantLib::DateGeneration::Forward, false);
		fixedSchedule = QuantLib::Schedule(startDate, maturityDate,
			instrument->PayLegFrequency->Value * QuantLib::TimeUnit((int)instrument->PayLegFrequency->Unit + 1),
			QuantLib::TARGET(), QuantLib::Following, QuantLib::Following, QuantLib::DateGeneration::Forward, false);
	}

	boost::shared_ptr<QuantLib::IborIndex> ibor(new QuantLib::Euribor(3 * QuantLib::Months
		, QuantLib::Handle<QuantLib::YieldTermStructure>(*curve_)));
	ibor->addFixing(ibor->fixingDate(referenceDate), 0., true);

	instrument_ = new QuantLib::VanillaSwap(
		irsType, instrument->Nominal, fixedSchedule, instrument->FixedRate, QuantLib::Thirty360(),
		floatSchedule, ibor, instrument->FloatingSpread, QuantLib::Actual360());
}

void InstrumentFactory::OnVisit(RiskAnalysisTool::Instruments::EquitySwap ^instrument) {
	if (instrument_) {
		delete instrument_;
		instrument_ = nullptr;
	}
	//TODO: create it
	QuantLib::Date startDate = Utilities::DateHelper::ToQLDate(instrument->StartDate);
	QuantLib::Date maturityDate = Utilities::DateHelper::ToQLDate(instrument->MaturityDate);
	QuantLib::EquitySwap::Type ersType;

	String^ paySymbol = instrument->PaySymbol;
	QuantLib::Date referenceDate
		= Utilities::DateHelper::ToQLDate(request_->EvaluationData);

	if (paySymbol->Contains("$"))
		ersType = QuantLib::EquitySwap::Type::Payer;
	else
		ersType = QuantLib::EquitySwap::Type::Receiver;

	//boost::shared_ptr<QuantLib::YieldTermStructure> curve;


	instrument_ = new QuantLib::EquitySwap(
		ersType, startDate, referenceDate, instrument->StartPrice, instrument->SpotPrice,
		instrument->Amount, 0., 0., instrument->FixedRate, instrument->DividendYield, (*curve_)->dayCounter().yearFraction(startDate, maturityDate)
		, instrument->Volatility, *curve_);
}

void InstrumentFactory::OnVisit(RiskAnalysisTool::Instruments::CrossCurrencySwap ^instrument) {
	if (instrument_) {
		delete instrument_;
		instrument_ = nullptr;
	}
	//TODO: create it
	QuantLib::Date startDate = Utilities::DateHelper::ToQLDate(instrument->StartDate);
	QuantLib::Date maturityDate = Utilities::DateHelper::ToQLDate(instrument->MaturityDate);
	String^ paySymbol = instrument->PaySymbol;
	QuantLib::Date referenceDate = Utilities::DateHelper::ToQLDate(request_->EvaluationData);

	QuantLib::CrossCurrencySwap::Type ccsType;
	double exchangeRate = instrument->PayNominal / instrument->ReceiveNominal;
	if (paySymbol->Contains(request_->DomesticCurrency)){
		ccsType = QuantLib::CrossCurrencySwap::Type::payDomestic;
	}
	else {
		ccsType = QuantLib::CrossCurrencySwap::Type::payForeign;
		exchangeRate = 1. / exchangeRate;
	}

	QuantLib::Schedule paySchedule(startDate, maturityDate,
		instrument->PayLegFrequency->Value * QuantLib::TimeUnit((int)instrument->PayLegFrequency->Unit + 1),
		QuantLib::TARGET(), QuantLib::Following, QuantLib::Following, QuantLib::DateGeneration::Forward, false);

	QuantLib::Schedule receiveSchedule(startDate, maturityDate,
		instrument->ReceiveLegFrequency->Value * QuantLib::TimeUnit((int)instrument->ReceiveLegFrequency->Unit + 1),
		QuantLib::TARGET(), QuantLib::Following, QuantLib::Following, QuantLib::DateGeneration::Forward, false);

	instrument_ = new QuantLib::CrossCurrencySwap(
		ccsType, exchangeRate, instrument->ExchangeVolatility, referenceDate,
		instrument->PayNominal, instrument->ReceiveNominal, paySchedule, instrument->PayRate, QuantLib::Actual360(),
		receiveSchedule, instrument->ReceiveRate, QuantLib::Actual360(), boost::make_shared<QuantLib::Currency>(QuantLib::CNYCurrency()),
		boost::make_shared<QuantLib::Currency>(QuantLib::CNYCurrency()));
}
