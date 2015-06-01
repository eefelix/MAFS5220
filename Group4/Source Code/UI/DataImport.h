#pragma warning(disable:4819)
#include <windows.h>
#include <ql/quantlib.hpp>
#include <fstream>

using namespace std;
using namespace QuantLib;

//struct EquityForwardData {
//	string Code_;
//	string Underlying_;
//	Real Strike_;
//	int SettlementDay_;
//	int ValueDay_;
//	int MaturityYear_;
//	string PositionType_;
//};

struct EquityMarketData {
	string Code_;
	Real Spot_;
	Real Dividend_;
	Real Drift_;
	Real Volatility_;
	Real Cds6Mths_;
	Real Cds1Yrs_;
	Real Cds2Yrs_;
	Real Cds3Yrs_;
	Real Cds4Yrs_;
	Real Cds5Yrs_;
	Real Cds7Yrs_;
	Real Cds10Yrs_;
};

//class InstrumentData
//{
//private:
//
//public:
//	EquityForwardData LoadEquityForward(string tradeRef);
//};

class MarketData
{
private:
	std::string GetCurrentWorkingDirectory();
	bool IsFileExist(string filePath);

public:
	EquityMarketData LoadEquity(string bbgCode);
};