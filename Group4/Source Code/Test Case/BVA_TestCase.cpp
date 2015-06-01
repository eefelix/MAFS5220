#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <iostream>
#include "BVA_calculator.h"
#include "DataImport.h"

using namespace std;
using namespace QuantLib;



void TestBVA() {
	// ==========================================
	// *** specify date/calendar config below ***
	// ==========================================

	Calendar calendar = TARGET();
	DayCounter dayCounter = Actual365Fixed();   //<== it is agreed that 365 basis to used?
	Date todaysDate(10, July, 2007);

	// ============================================================================================
	// *** BVA Test (Hazard rate calibration): CDS spread of Investor & Counterparty are equal  ***
	// ============================================================================================

	// *** initialize default time related objects below ***

	// CDS premium/schedule user input

	vector<Real> oInvestorCdsSpread;
	vector<Real> oCounterptyCdsSpread;

	Real investor_cds[] = {0.0016, 0.0029, 0.0045, 0.005, 0.0058};
	for (int i = 0; i < sizeof(investor_cds); i++)
		oInvestorCdsSpread.push_back(investor_cds[i]);

	Real counterparty_cds[] = {0.0016, 0.0029, 0.0045, 0.005, 0.0058};
	for (int i = 0; i < sizeof(counterparty_cds); i++)
		oCounterptyCdsSpread.push_back(counterparty_cds[i]);

	vector<Period> tenors;
	tenors.push_back(1 * Years);
	tenors.push_back(3 * Years);
	tenors.push_back(5 * Years);
	tenors.push_back(7 * Years);
	tenors.push_back(10 * Years);

	vector<Date> maturities;
	for (int i = 0; i < tenors.size(); i++) 
		maturities.push_back(calendar.adjust(todaysDate + tenors[i], Following));

	// initiazlize investor hazard rate / default time object
	HazardRateCalibrate oInvestorHazardRate = HazardRateCalibrate(todaysDate);
	oInvestorHazardRate.SetYieldRate(0.000);
	oInvestorHazardRate.SetRecoveryRate(0.40);
	oInvestorHazardRate.SetCdsSpread(oInvestorCdsSpread);
	oInvestorHazardRate.SetCdsTenor(tenors);
	oInvestorHazardRate.PerformCalibration();
	DefaultTimeGenerator oInvestorDefaultTime(todaysDate);
	oInvestorDefaultTime.SetHazardRates(maturities, oInvestorHazardRate.oHazardRates);

	// initiazlize counterparty hazard rate / default time object

	HazardRateCalibrate oCounterptyHazardRate = HazardRateCalibrate(todaysDate);
	oCounterptyHazardRate.SetYieldRate(0.000);
	oCounterptyHazardRate.SetRecoveryRate(0.40);
	oCounterptyHazardRate.SetCdsSpread(oCounterptyCdsSpread);
	oCounterptyHazardRate.SetCdsTenor(tenors);
	oCounterptyHazardRate.PerformCalibration();
	DefaultTimeGenerator oCounterptyDefaultTime(todaysDate);
	oCounterptyDefaultTime.SetHazardRates(maturities, oCounterptyHazardRate.oHazardRates);

	
	// *** details of the configuration of the forward ***

	int no_of_Date_path = 10000;
	int no_of_stock_path = 10;
	Real S0 = 99;
	Real K = 100;
	Natural settle = 2;
	Spread q = 0.0;
	Real spot_dividend = 0;
	Rate r = 0.01;
	Volatility sigma = 0.3287;
	Position::Type type(Position::Long);  //please change it for testing the different direction of the position
	Date valueDate = calendar.adjust(todaysDate + 2 * Days, Following);
	Date maturity = calendar.adjust(todaysDate + 1 * Years, Following);
	Real own_recovery = 0.4;
	Real counterparty_recovery = 0.4;
	Real coll_inv = 0;
	Real coll_cpt = 0;
	Real coll_inv_yield = 0;
	Real coll_cpt_yield = 0;


	cout << "Value date: " << todaysDate << endl;
	cout << "The configuration of the forward: " << endl;
	cout << "Spot: " << S0 << endl;
	cout << "Strike: " << K << endl;
	cout << "risk free rate: " << r << endl;
	cout << "Dividend: " << spot_dividend << endl;
	cout << "Dividend yield: " << q << endl;
	cout << "Long/ Short: " << "Long" << endl;
	cout << "Tenor (year): " << 1 << endl;
	cout << endl;
	/*
	vector<Real> testing = oInvestorHazardRate.oHazardRates;
	
	int i=0;
	std::vector<Real>::iterator at;
	for (at = testing.begin()+1; at != testing.end(); ++at)
	{   
		cout << "Investor's CDS spreads: " << oInvestorCdsSpread[i] << ", hazard rate: " << *(at) << endl;
		i++;
	}
	testing.clear();
	testing = oCounterptyHazardRate.oHazardRates;
	
	i=0;

	for (at = testing.begin()+1; at != testing.end(); ++at)
	{   
		cout << "Counterparty's CDS spreads: " << oCounterptyCdsSpread[i] << ", hazard rate: " << *(at) << endl;
		i++;
	}
	cout << endl;
	*/
	//initial the class for Monte Carlo simulation of the forward and stock with the input above 
	MCgenerator f(todaysDate, maturity, valueDate, settle, calendar, S0, dayCounter, K, q, spot_dividend, r, coll_inv_yield, coll_cpt_yield, sigma, type);

	// *** calculate CVA/DVA ***
	//initial the class BVA_calculator for calculating the CVA and DVA with the class initialized above
	BVA_calculator newbva(&f, &oInvestorHazardRate, &oCounterptyHazardRate, no_of_Date_path, no_of_stock_path,coll_inv,coll_cpt);
	
	cout << "Test Case 1: " << endl;
	cout << "BVA by the Hazard rate calibration with collateral posted (No re-hypothecation)" << endl;
	cout << endl;
	
	vector<map<string, Real>> CVA_DVA;

	std::map<string,Real>::iterator ot;
	std::vector<map<string,Real>>::iterator it;

	//calculated the CBVA
	CVA_DVA = newbva.BVA_cal(maturity,'y',coll_inv,coll_cpt);

	//print the results
	for (it = CVA_DVA.begin(); it != CVA_DVA.end(); ++it)
	{
		for (ot = (*it).begin(); ot != (*it).end(); ++ot)
			cout << (ot->first) << " " << (ot->second) << endl;
	}
	cout << endl << endl;

	CVA_DVA.clear();

	cout << "Test Case 2: " << endl;
	cout << "BVA by the Hazard rate calibration with collateral posted" << endl;
	cout << "(Collateral re-hypothecation both investor and counterparty with 50% collateral recovery rate)" << endl;
	cout << endl;

	coll_inv = 0.6;
	coll_cpt = 0.6;

	CVA_DVA = newbva.BVA_cal(maturity,'y',coll_inv,coll_cpt);

	for (it = CVA_DVA.begin(); it != CVA_DVA.end(); ++it)
	{
		for (ot = (*it).begin(); ot != (*it).end(); ++ot)
			cout << (ot->first) << " " << (ot->second) << endl;
	}
	cout << endl << endl;
	
	CVA_DVA.clear();
	//calculate the BVA without collateral
	cout << "Test Case 3:" << endl;
	cout << "BVA by the Hazard rate calibration without collateral posted" << endl;
	
	CVA_DVA = newbva.BVA_cal(maturity,'n',coll_inv,coll_cpt);
	
	for (it = CVA_DVA.begin(); it != CVA_DVA.end(); ++it)
	{
		for (ot = (*it).begin(); ot != (*it).end(); ++ot)
			cout << (ot->first) << " " << (ot->second) << endl;
	}
	cout << endl << endl;

	// ============================================================================================
	// *** BVA Test (AT1P CDS calibration): CDS spread of Investor & Counterparty are equal  ***
	// ============================================================================================
	vector<Real> discount_rate;
	for(int i=0; i<tenors.size(); i++)
		discount_rate.push_back(5); //set the discount rate as 5%

	//Setup the AT1P model calibration for investor and counterparty
	At1pCalibration investor = At1pCalibration(oInvestorCdsSpread,discount_rate,todaysDate,tenors);
	At1pCalibration counterparty = At1pCalibration(oCounterptyCdsSpread,discount_rate,todaysDate,tenors);
	BVA_calculator at1p_bva(&f, &investor, &counterparty,own_recovery,counterparty_recovery, coll_inv, coll_cpt, no_of_Date_path, no_of_stock_path);
	
	vector<map<string, Real>> AT1P_CVA_DVA;

	//calculated the BVA with collateral
	AT1P_CVA_DVA = at1p_bva.BVA_cal_AT1P(maturity,'y',coll_inv,coll_cpt);


	vector<Real> calibrated_volat = investor.get_cal_vol();
	//print the results
	
	for (int i = 0; i < tenors.size(); ++i)
	{
		cout << "Investor's AT1P calibrated volatilities: " << calibrated_volat[i] << endl;
	}
	
	calibrated_volat.clear();
    calibrated_volat = counterparty.get_cal_vol();

	for (int i = 0; i < tenors.size(); ++i)
	{
		cout << "Counterparty's AT1P calibrated volatilities: " << calibrated_volat[i] << endl;
	}
	cout << endl; 

	cout << "Test Case 4: " << endl;
	cout << "BVA by the AT1P CDS calibration with collateral posted (No re-hypothecation)" << endl;
	cout << endl;

	for (it = AT1P_CVA_DVA.begin(); it != AT1P_CVA_DVA.end(); ++it)
	{
		for (ot = (*it).begin(); ot != (*it).end(); ++ot)
			cout << (ot->first) << " " << (ot->second) << endl;
	}
	cout << endl;
	cout << endl;
	AT1P_CVA_DVA.clear();

	//calculate the BVA without collatreral
	AT1P_CVA_DVA = at1p_bva.BVA_cal_AT1P(maturity,'n',coll_inv,coll_cpt);
	cout << "Test Case 5: " << endl;
	cout << "BVA by the AT1P CDS calibration without collateral posted" << endl;

	for (it = AT1P_CVA_DVA.begin(); it != AT1P_CVA_DVA.end(); ++it)
	{
		for (ot = (*it).begin(); ot != (*it).end(); ++ot)
			cout << (ot->first) << " " << (ot->second) << endl;
	}
	cout << endl;

}

void CalculateBVAWWR() {

	cout << endl;
	cout << " ======================================================= " << endl;
	// ==========================================
	// *** specify date/calendar config below ***
	// ==========================================

	Calendar calendar = TARGET();
	DayCounter dayCounter = Actual365Fixed();   //<== agreed that 365 basis to used
	Date todaysDate(13, May, 2015);
	
	// ========================================================
	// *** details of the configuration of the forward ***
	// ========================================================

	int no_of_Date_path = 1000;		// input parameter to BVACalculatorWWR 
	int no_of_stock_path = 10;		// input parameter to BVACalculatorWWR but actually not used in BVA WWR
	Real S0 = 100;
	Real K = 100;
	Natural settle = 2;
	Spread q = 0.0;
	Real spot_dividend = 0.00;
	Rate r = 0.00;
	Volatility sigma = 0.15;
	Position::Type type(Position::Long);  //please change it for testing the different direction of the position
	Position::Type typeInvestorCalibrate;
	Position::Type typeCounterpartyCalibrate;
	Date valueDate = calendar.adjust(todaysDate + 2 * Days, Following);
	Date maturity = calendar.adjust(todaysDate + 1 * Years, Following);

	//initial the class for Monte Carlo simulation of the forward and stock with the input above 
	MCgenerator oMCgenerator(todaysDate, maturity, valueDate, settle, calendar, S0, dayCounter, K, q, spot_dividend, r, sigma, type);

	// =====================================================================================================================
	// *** BVA WWR Test Case: CDS spread of Investor & Counterparty are equal  ***
	// =====================================================================================================================

	// CDS premium/schedule user input

	Real investor_cds_1[] = { 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125 };
	vector<Real> oInvestorCdsSpread(std::begin(investor_cds_1), std::end(investor_cds_1));		// check if can run in VS2012

	Real counterparty_cds_1[] = { 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125, 0.0125 };
	vector<Real> oCounterptyCdsSpread(std::begin(counterparty_cds_1), std::end(counterparty_cds_1));		// check if can run in VS2012

	vector<Period> tenors;
	tenors.push_back(6 * Months);
	tenors.push_back(1 * Years);
	tenors.push_back(2 * Years);
	tenors.push_back(3 * Years);
	tenors.push_back(4 * Years);
	tenors.push_back(5 * Years);
	tenors.push_back(7 * Years);
	tenors.push_back(10 * Years);

	vector<Date> maturities;
	for (int i = 0; i < tenors.size(); i++) {
		maturities.push_back(calendar.adjust(todaysDate + tenors[i],
			Following));
	}

	// initialize the HazardRateWWR objects for investor & counterparty
	if (type == Position::Long) {
		typeInvestorCalibrate = Position::Short;		// !!! may a bit counter-intuitive but plz remember we need to calibrate investor hazard rate using payoff from counterparty perspective 
		typeCounterpartyCalibrate = Position::Long;		// !!! may a bit counter-intuitive but plz remember we need to calibrate counterparty hazard rate using payoff from investor perspective 
	}
	else if (type == Position::Short) {
		typeInvestorCalibrate = Position::Long;
		typeCounterpartyCalibrate = Position::Short;
	}

	// =====================================================================================================================
	// *** BVA WWR Test Case 1: No WWR (b = 0)
	// =====================================================================================================================

	cout << "Value Adjustment results when there's no WWR. i.e. b = 0" << endl;

	HazardRateWWR oInvestorHazardRateWWR = HazardRateWWR(todaysDate);
	oInvestorHazardRateWWR.SetCdsPayDate(maturities);
	oInvestorHazardRateWWR.SetCdsSpread(oInvestorCdsSpread);
	oInvestorHazardRateWWR.SetRecoveryRate(0.4);
	oInvestorHazardRateWWR.SetForwardType(typeInvestorCalibrate);
	oInvestorHazardRateWWR.SetSensitivityB(0.00);
	//oInvestorHazardRateWWR.SetWWRHazardModel(EXPONENTIALLY);		// Select either exponential or linear for investor hazard rate function
	oInvestorHazardRateWWR.SetWWRHazardModel(LINEARLY);		// Select either exponential or linear for investor hazard rate function
	oInvestorHazardRateWWR.SimulateStockPrices(oMCgenerator);
	cout << "Start calibrating a(ti) of HazardRateWWR of investor " << endl;
	oInvestorHazardRateWWR.PerformCalibration();

	HazardRateWWR oCptyHazardRateWWR = HazardRateWWR(todaysDate);
	oCptyHazardRateWWR.SetCdsPayDate(maturities);
	oCptyHazardRateWWR.SetCdsSpread(oCounterptyCdsSpread);
	oCptyHazardRateWWR.SetRecoveryRate(0.4);
	oCptyHazardRateWWR.SetForwardType(typeCounterpartyCalibrate);
	oCptyHazardRateWWR.SetSensitivityB(0.00);
	//oCptyHazardRateWWR.SetWWRHazardModel(EXPONENTIALLY);		// Select either exponential or linear for counterparty hazard rate function
	oCptyHazardRateWWR.SetWWRHazardModel(LINEARLY);		// Select either exponential or linear for counterparty hazard rate function
	oCptyHazardRateWWR.SimulateStockPrices(oMCgenerator);
	cout << "Start calibrating a(ti) of HazardRateWWR of counterparty " << endl;
	oCptyHazardRateWWR.PerformCalibration();

	// initial the class BVACalculatorWWR for calculating the CVA and DVA with the class initialized above

	BVACalculatorWWR oBVACalculatorWWR(&oMCgenerator, &oInvestorHazardRateWWR, &oCptyHazardRateWWR, no_of_Date_path, no_of_stock_path);

	map<string, Real> valueAdjustment;

	// Calculate the BVA

	valueAdjustment = oBVACalculatorWWR.BVA_cal(maturity);

	Real oUDVA;
	Real oUCVA;
	oUDVA = valueAdjustment.find("UDVA")->second;
	oUCVA = valueAdjustment.find("UCVA")->second;

	// Print the results
	cout << "UDVA to investor perspective = " << oUDVA << endl;
	cout << "UCVA to investor perspective = " << oUCVA << " vs John Hull/Alan White's result CVA = 0.048" << endl;
	cout << "BVA = UDVA - UCVA = " << oUDVA - oUCVA << endl << endl;

	// =====================================================================================================================
	// *** BVA WWR Test Case 2: WWR (b = 0.03)
	// =====================================================================================================================

	cout << "Value Adjustment results when there's WWR. i.e. b = 0.03" << endl;

	oInvestorHazardRateWWR.SetSensitivityB(0.03);
	cout << "Start calibrating a(ti) of HazardRateWWR of investor " << endl;
	oInvestorHazardRateWWR.PerformCalibration();
	oCptyHazardRateWWR.SetSensitivityB(0.03);
	cout << "Start calibrating a(ti) of HazardRateWWR of counterparty " << endl;
	oCptyHazardRateWWR.PerformCalibration();

	valueAdjustment = oBVACalculatorWWR.BVA_cal(maturity);
	oUDVA = valueAdjustment.find("UDVA")->second;
	oUCVA = valueAdjustment.find("UCVA")->second;

	// Print the results
	cout << "UDVA to investor perspective = " << oUDVA << endl;
	cout << "UCVA to investor perspective = " << oUCVA << " vs John Hull/Alan White's result CVA = 0.074" << endl;
	cout << "BVA = UDVA - UCVA = " << oUDVA - oUCVA << endl << endl;

	// =====================================================================================================================
	// *** BVA WWR Test Case 3: RWR (b = -0.03)
	// =====================================================================================================================

	cout << "Value Adjustment results when there's WWR. i.e. b = -0.03" << endl;

	oInvestorHazardRateWWR.SetSensitivityB(-0.03);
	cout << "Start calibrating a(ti) of HazardRateWWR of investor " << endl;
	oInvestorHazardRateWWR.PerformCalibration();
	oCptyHazardRateWWR.SetSensitivityB(-0.03);
	cout << "Start calibrating a(ti) of HazardRateWWR of counterparty " << endl;
	oCptyHazardRateWWR.PerformCalibration();

	valueAdjustment = oBVACalculatorWWR.BVA_cal(maturity);
	oUDVA = valueAdjustment.find("UDVA")->second;
	oUCVA = valueAdjustment.find("UCVA")->second;

	// Print the results
	cout << "UDVA to investor perspective = " << oUDVA << endl;
	cout << "UCVA to investor perspective = " << oUCVA << " vs John Hull/Alan White's result CVA = 0.03" << endl;
	cout << "BVA = UDVA - UCVA = " << oUDVA - oUCVA << endl;

}



int main(int, char*[]) {

	TestBVA();

	CalculateBVAWWR();

	system("pause");
}

