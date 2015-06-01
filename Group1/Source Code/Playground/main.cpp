#include <pch.h>

#pragma managed

#include <Calculation/Utilities/cliext.hpp>
#include <iostream>

void OnDomainUnload(System::Object ^sender, System::EventArgs ^e);
class NativeClass {
public:
	int v;
	NativeClass(int x) {
		std::cout << "[" << x << "] In NativeClass::NativeClass() " << std::endl;
		v = x;
	}
	NativeClass(const NativeClass &) = delete;
	~NativeClass() {
		std::cout << "[" << v << "] In NativeClass::~NativeClass() " << std::endl;
	}
};

public ref class Holder {
public:
	Holder(int x) {
		p_ = boost::make_shared<NativeClass>(x);
	}
public:
	boost::shared_ptr<NativeClass> getT() {
		return p_.native();
	}
private:
	cliext::shared_ptr<NativeClass> p_;
};

ref class Global
{
public:
	static Holder ^h = nullptr;
};


int main3() {

	Holder ^h1 = gcnew Holder(1);
	Holder ^h2 = gcnew Holder(2);
	Holder ^h3 = gcnew Holder(3);
	boost::shared_ptr<NativeClass> p = h1->getT();

	// Should not delete the native instance 1 until function exits, because p is holding it
	delete h1;

	// Now instance 2 should be deleted
	delete h2;

	Global::h = h3;

	std::cout << "main() returning" << std::endl;

	//Instance 3 should be deleted lastly happens lastly, before program exists
	return 0;
}

// irs request
int main(array<System::String ^> ^args) {
	using namespace System;

	using namespace RiskAnalysisTool::Instruments;
	using namespace RiskAnalysisTool::Requests;
	using namespace RiskAnalysisTool::Results;
	using namespace RiskAnalysisTool::Calculation;
	using namespace RiskAnalysisTool::Time;
	using namespace System::Collections::Generic;
	// prepare market data
	array<double, 2>^ correlationMatrix = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 }
	};

	//double riskFreeRate = 0.01;
	double investorRecoveryRate = 0.4;
	double counterpartyRecoveryRate = 0.4;

	// case one
	array<Bond^>^ bonds = gcnew array<Bond^>(22);
	array<DateTime>^ bondTenor = {
		DateTime(2015, 3, 18), DateTime(2015, 6, 18), DateTime(2015, 9, 18), DateTime(2015, 12, 18),
		DateTime(2016, 3, 18), DateTime(2016, 6, 20), DateTime(2016, 9, 19), DateTime(2016, 12, 19),
		DateTime(2017, 3, 20), DateTime(2017, 6, 19), DateTime(2017, 9, 18), DateTime(2017, 12, 18),
		DateTime(2018, 3, 19), DateTime(2018, 6, 18), DateTime(2018, 9, 18), DateTime(2018, 12, 18),
		DateTime(2019, 3, 18), DateTime(2019, 6, 18), DateTime(2019, 9, 18), DateTime(2019, 12, 18),
		DateTime(2020, 3, 18), DateTime(2021, 3, 18)
	};

	array<double>^ dfs = {
		1., 0.99931, 0.998298, 0.996803, 0.994802, 0.99218,
		0.989108, 0.985539, 0.981484, 0.977083, 0.972384,
		0.967452, 0.962292, 0.956992, 0.951799, 0.946447,
		0.940943, 0.935505, 0.929912, 0.924224, 0.918389,
		0.895473
	};

	//// case two
	//array<DateTime>^ bondTenor = {
	//	DateTime(2015, 5, 5), DateTime(2015, 8, 5), DateTime(2015, 11, 5), DateTime(2016, 2, 5),
	//	DateTime(2016, 5, 5), DateTime(2016, 8, 5), DateTime(2016, 11, 7), DateTime(2017, 2, 6),
	//	DateTime(2017, 5, 5), DateTime(2017, 8, 7), DateTime(2017, 11, 6), DateTime(2018, 2, 5),
	//	DateTime(2018, 5, 8), DateTime(2018, 8, 6), DateTime(2018, 11, 5), DateTime(2019, 2, 5),
	//	DateTime(2019, 5, 7), DateTime(2019, 8, 5), DateTime(2019, 11, 5), DateTime(2020, 2, 5),
	//	DateTime(2020, 5, 5), DateTime(2021, 5, 5)
	//};
	//array<double>^ dfs = {
	//	1., 0.999295, 0.998334, 0.996943, 0.995086, 0.992695,
	//	0.989714, 0.986306, 0.982583, 0.978212, 0.973639,
	//	0.968765, 0.963582, 0.958442, 0.95315, 0.947538,
	//	0.941737, 0.936129, 0.9302, 0.924081, 0.91792,
	//	0.892829
	//};

	for (int i = 0; i < 22; ++i) {
		bonds[i] = gcnew Bond;
		bonds[i]->MaturityDate = bondTenor[i];
		bonds[i]->Price = dfs[i];
		bonds[i]->Principle = 1.;
		bonds[i]->Symbol = "$USD";
	}

	array<CreditDefaultSwap^>^ issuerCds = gcnew array<CreditDefaultSwap^>(5);
	array<Period ^>^ issuerCdsTenor = { gcnew Period(PeriodUnit::Year, 1), gcnew Period(PeriodUnit::Year, 2),
		gcnew Period(PeriodUnit::Year, 3), gcnew Period(PeriodUnit::Year, 4), gcnew Period(PeriodUnit::Year, 5) };
	array<double>^ issuerCdsSpread = { 0.031621, 0.050727, 0.090606, 0.117975, 0.126978 };
	/*array<double>^ issuerCdsSpread = { 0.008001, 0.025403, 0.050331, 0.062066, 0.069991 };*/

	for (int i = 0; i < issuerCds->Length; ++i) {
		issuerCds[i] = gcnew CreditDefaultSwap;
		issuerCds[i]->Tenor = issuerCdsTenor[i];
		issuerCds[i]->Spread = issuerCdsSpread[i];
		issuerCds[i]->Symbol = "CDS:ISSUER";
	}

	array<CreditDefaultSwap^>^ investorCds = gcnew array<CreditDefaultSwap^>(5);
	array<Period ^>^ investorCdsTenor = { gcnew Period(PeriodUnit::Year, 1), gcnew Period(PeriodUnit::Year, 2),
		gcnew Period(PeriodUnit::Year, 3), gcnew Period(PeriodUnit::Year, 4), gcnew Period(PeriodUnit::Year, 5) };
	array<double>^ investorCdsSpread = { 0.034021, 0.04743, 0.059706, 0.065624, 0.077995 };
	/*array<double>^ investorCdsSpread = { 0.003303, 0.004424, 0.005720, 0.006850, 0.008279 }*/;

	for (int i = 0; i < investorCds->Length; ++i) {
		investorCds[i] = gcnew CreditDefaultSwap;
		investorCds[i]->Tenor = investorCdsTenor[i];
		investorCds[i]->Spread = investorCdsSpread[i];
		investorCds[i]->Symbol = "CDS:INVESTOR";
	}

	ICollection<Instrument^>^ marketData = gcnew List < Instrument^ > ;
	for each (Bond^ bond in bonds)
		marketData->Add(bond);

	for each (CreditDefaultSwap^ cds in issuerCds)
		marketData->Add(cds);

	for each (CreditDefaultSwap^ cds in investorCds)
		marketData->Add(cds);

	// set up instrument
	InterestRateSwap^ irs = gcnew InterestRateSwap;
	// case one
	irs->Price = 0.;
	irs->StartDate = DateTime(2015, 3, 18);
	irs->EvaluationData = DateTime(2015, 3, 18);
	irs->MaturityDate = DateTime(2020, 3, 18);
	irs->Symbol = "Interest Rate Swap";
	irs->Nominal = 1.0;
	irs->FixedRate = 0.01551176;
	irs->FloatingSpread = 0.0;
	irs->IsPayFloating = false;
	irs->ReceiveLegFrequency = gcnew Period(PeriodUnit::Month, 3);
	irs->PayLegFrequency = gcnew Period(PeriodUnit::Month, 6);
	irs->PaySymbol = "$USD";
	irs->ReceiveSymbol = "$USD";

	//// case two
	//irs->Price = 0.;
	//irs->StartDate = DateTime(2015, 5, 5);
	//irs->EvaluationData = DateTime(2015, 5, 5);
	//irs->MaturityDate = DateTime(2020, 5, 5);
	//irs->Symbol = "Interest Rate Swap";
	//irs->Nominal = 1.0;
	//irs->FixedRate = 0.01672222;
	//irs->FloatingSpread = 0.0;
	//irs->IsPayFloating = true;
	//irs->FixedLegFrequency = gcnew Period(PeriodUnit::Month, 6);
	//irs->FloatingLegFrequency = gcnew Period(PeriodUnit::Month, 3);
	//irs->PaySymbol = "$USD";
	//irs->ReceiveSymbol = "$USD";

	ICollection<Instrument^>^ portfolio = gcnew List < Instrument^ > ;
	portfolio->Add(irs);

	// initiate a Request
	BvaRequest^ request = gcnew BvaRequest;
	request->EvaluationData = DateTime(2015, 3, 18);
	/*request->EvaluationData = DateTime(2015, 5, 5);*/
	//request->RiskFreeRate = riskFreeRate;
	request->InvestorRecoveryRate = investorRecoveryRate;
	request->CounterpartyRecoveryRate = counterpartyRecoveryRate;
	request->Portfolio = portfolio;
	request->MarketData = marketData;
	request->DomesticCurrency = "$USD";

	// set up computation engine
	ComputeEngine^ engine = gcnew ComputeEngine;
	BvaResult^ result = safe_cast<BvaResult^>(engine->calculate(request));
	Console::WriteLine("BVA is: {0}", result->Cva - result->Dva);

	return 0;
}

//int main_ccs(array<System::String ^> ^args) {
//	using namespace System;
//
//	using namespace RiskAnalysisTool::Instruments;
//	using namespace RiskAnalysisTool::Requests;
//	using namespace RiskAnalysisTool::Results;
//	using namespace RiskAnalysisTool::Calculation;
//	using namespace RiskAnalysisTool::Time;
//	using namespace System::Collections::Generic;
//	// prepare market data
//	array<double, 2>^ correlationMatrix = {
//		{ 1, 0, 0, 0 },
//		{ 0, 1, 0, 0 },
//		{ 0, 0, 1, 0 },
//		{ 0, 0, 0, 1 }/*,
//					  { 0, 0, 0, 0, 1 },*/
//	};
//	//double riskFreeRate = 0.04;
//	double investorRecoveryRate = 0.4;
//	double counterpartyRecoveryRate = 0.4;
//
//	array<Bond^>^ bonds = gcnew array<Bond^>(42);
//	array<DateTime>^ bondTenor = {
//		DateTime(2015, 3, 18), DateTime(2015, 9, 18), DateTime(2016, 3, 18),
//		DateTime(2016, 9, 18), DateTime(2017, 3, 18), DateTime(2017, 9, 18),
//		DateTime(2018, 3, 18), DateTime(2018, 9, 18), DateTime(2019, 3, 18),
//		DateTime(2019, 9, 18), DateTime(2020, 3, 18), DateTime(2021, 3, 18),
//		DateTime(2022, 3, 18), DateTime(2023, 3, 18), DateTime(2024, 3, 18),
//		DateTime(2025, 3, 18), DateTime(2026, 3, 18), DateTime(2027, 3, 18),
//		DateTime(2030, 3, 18), DateTime(2035, 3, 18), DateTime(2040, 3, 18)
//	};
//	array<double>^ dfs = {
//		1., 0.998298, 0.99479, 0.989097, 0.981462, 0.972366,
//		0.962267, 0.951817, 0.941055, 0.930079, 0.918616,
//		0.895473, 0.872712, 0.850075, 0.827519, 0.805479,
//		0.783812, 0.762483, 0.70229, 0.613123, 0.5379
//	};
//
//	for (int i = 0; i < 21; ++i) {
//		bonds[i] = gcnew Bond;
//		bonds[i]->MaturityDate = bondTenor[i];
//		bonds[i]->Price = dfs[i];
//		bonds[i]->Principle = 1.;
//		bonds[i]->Symbol = "$CNY";
//	}
//
//	for (int i = 21; i < 42; ++i) {
//		bonds[i] = gcnew Bond;
//		bonds[i]->MaturityDate = bondTenor[i - 21];
//		bonds[i]->Price = dfs[i - 21];
//		bonds[i]->Principle = 1.;
//		bonds[i]->Symbol = "$USD";
//	}
//
//	array<CreditDefaultSwap^>^ issuerCds = gcnew array<CreditDefaultSwap^>(7);
//	array<Period ^>^ issuerCdsTenor = { gcnew Period(PeriodUnit::Year, 1), gcnew Period(PeriodUnit::Year, 2),
//		gcnew Period(PeriodUnit::Year, 3), gcnew Period(PeriodUnit::Year, 4), gcnew Period(PeriodUnit::Year, 5),
//		gcnew Period(PeriodUnit::Year, 7), gcnew Period(PeriodUnit::Year, 10) };
//	/*array<Period ^>^ issuerCdsTenor = { gcnew Period(PeriodUnit::Month, 3), gcnew Period(PeriodUnit::Month, 6),
//	gcnew Period(PeriodUnit::Year, 1), gcnew Period(PeriodUnit::Year, 2)};*/
//	array<double>^ issuerCdsSpread = { 26.766, 35.064, 42.953, 53.080, 68.203, 91.453, 111.109 };
//	const double bps = 10000.;
//	for (int i = 0; i < issuerCdsSpread->Length; ++i)
//		issuerCdsSpread[i] /= bps;
//
//
//	for (int i = 0; i < issuerCds->Length; ++i) {
//		issuerCds[i] = gcnew CreditDefaultSwap;
//		issuerCds[i]->Tenor = issuerCdsTenor[i];
//		issuerCds[i]->Spread = issuerCdsSpread[i];
//		issuerCds[i]->Symbol = "CDS:ISSUER";
//	}
//
//	array<CreditDefaultSwap^>^ investorCds = gcnew array<CreditDefaultSwap^>(7);
//	array<Period ^>^ investorCdsTenor = { gcnew Period(PeriodUnit::Year, 1), gcnew Period(PeriodUnit::Year, 2),
//		gcnew Period(PeriodUnit::Year, 3), gcnew Period(PeriodUnit::Year, 4), gcnew Period(PeriodUnit::Year, 5),
//		gcnew Period(PeriodUnit::Year, 7), gcnew Period(PeriodUnit::Year, 10) };
//	array<double>^ investorCdsSpread = { 30.235, 41.876, 54.027, 65.151, 78.457, 102.7, 122.061 };
//
//	for (int i = 0; i < investorCdsSpread->Length; ++i)
//		investorCdsSpread[i] /= bps;
//
//	for (int i = 0; i < investorCds->Length; ++i) {
//		investorCds[i] = gcnew CreditDefaultSwap;
//		investorCds[i]->Tenor = investorCdsTenor[i];
//		investorCds[i]->Spread = investorCdsSpread[i];
//		investorCds[i]->Symbol = "CDS:INVESTOR";
//	}
//
//	ICollection<Instrument^>^ marketData = gcnew List < Instrument^ >;
//	for each (Bond^ bond in bonds)
//		marketData->Add(bond);
//
//	for each (CreditDefaultSwap^ cds in issuerCds)
//		marketData->Add(cds);
//
//	for each (CreditDefaultSwap^ cds in investorCds)
//		marketData->Add(cds);
//
//	// set up instrument
//	CrossCurrencySwap^ ccs = gcnew CrossCurrencySwap;
//	ccs->Price = 0.;
//	ccs->StartDate = DateTime(2015, 3, 18);
//	ccs->MaturityDate = DateTime(2020, 3, 18);
//	ccs->Symbol = "Cross Currency Swap";
//	ccs->PayLegFrequency = gcnew Period(PeriodUnit::Month, 6);
//	ccs->ReceiveLegFrequency = gcnew Period(PeriodUnit::Month, 6);
//
//	ccs->PayNominal = 62046000;
//	ccs->PayRate = 3.294000 / 100.;
//	ccs->ReceiveNominal = 1e7;
//	ccs->ReceiveRate = 1.567594 / 100.;
//	//ccs->exchangeRate = 0.161171;
//	ccs->ExchangeVolatility = 0.12;
//	ccs->EvaluationData = DateTime(2015, 3, 18);
//	ccs->PaySymbol = "$USD";
//	ccs->ReceiveSymbol = "$CNY";
//
//	ICollection<Instrument^>^ portfolio = gcnew List < Instrument^ >;
//	portfolio->Add(ccs);
//
//	// initiate a TvaRequest
//	TvaRequest^ request = gcnew TvaRequest;
//	request->EvaluationData = DateTime(2015, 3, 18);
//	//request->RiskFreeRate = riskFreeRate;
//	request->InvestorRecoveryRate = investorRecoveryRate;
//	request->CounterpartyRecoveryRate = counterpartyRecoveryRate;
//	request->Correlation = correlationMatrix;
//	request->Portfolio = portfolio;
//	request->MarketData = marketData;
//	request->DomesticCurrency = "$CNY";
//
//	//BvaRequest^ request = gcnew BvaRequest;
//	//request->EvaluationData = DateTime(2015, 3, 18);
//	//request->RiskFreeRate = riskFreeRate;
//	//request->InvestorRecoveryRate = investorRecoveryRate;
//	//request->CounterpartyRecoveryRate = counterpartyRecoveryRate;
//	//request->Portfolio = portfolio;
//	//request->MarketData = marketData;
//	//request->DomesticCurrency = "$CNY";
//
//	// set up computation engine
//	ComputeEngine^ engine = gcnew ComputeEngine;
//	TvaResult^ result = safe_cast<TvaResult^>(engine->calculate(request));
//	Console::WriteLine("FVA is: {0}", result->Fva);
//
//	//Console::WriteLine("Welcome to the demo for BVA computation.");
//	//Console::WriteLine("Here you can calculate BVA for the following product:");
//	//Console::WriteLine("1. Interest Rate Swap");
//	//Console::WriteLine("2. Equity Return Swap");
//	//Console::WriteLine("3. Cross Currency Swap");
//	//Console::WriteLine("Please enter the number of product: ");
//	//int productIndex = int::Parse(Console::ReadLine());
//
//	//switch (productIndex){
//	//case 1:
//	//	irsBvaDemo();
//	//	break;
//	//case 2:
//	//	ersBvaDemo();
//	//	break;
//	//case 3:
//	//	ccsBvaDemo();
//	//	break;
//	//default:
//	//	Console::WriteLine("Invalid product index.");
//	//}
//
//	//Console::ReadLine();
//
//	return 0;
//}

int main_es(array<System::String ^> ^args) {

	return 0;
}