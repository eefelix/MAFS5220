#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <iostream>
#include "BVA_calculator.h"
#include "DataImport.h"

using namespace std;
using namespace QuantLib;


void RunGuiBVA(int iBVAModel) {

	try {
		// ==========================================
		// *** specify date/calendar config below ***
		// ==========================================

		Calendar calendar = TARGET();
		DayCounter dayCounter = Actual365Fixed();   //<== it is agreed that 365 basis to used?
		Date todaysDate = Date::todaysDate();

		// ============================================================================================
		// *** Coding for user input ***
		// ============================================================================================

		Real K;
		Real maturityYear;
		int iType; Position::Type type;
		string bbgCodeUly;
		string bbgCodeInvestor;
		string bbgCodeCounterparty;

		cout << endl << "Please input the BBG code of investor in the Equity Forward: " << endl << "[Ans]: ";
		std::getline(std::cin, bbgCodeInvestor);
		cout << endl << "Please input the BBG code of underlying stock of the Equity Forward: " << endl << "[Ans]: ";
		std::getline(std::cin, bbgCodeUly);
		cout << endl << "Please input the strike of the Equity Forward: " << endl << "[Ans]: ";
		cin >> K;
		cout << endl << "Please input the number of year to matiruty for the Equity Forward: (e.g. 1)" << endl << "[Ans]: ";
		cin >> maturityYear;
		cout << endl << "Please input the position type of the Equity Forward: (1=Long, 2=Short)" << endl << "[Ans]: ";
		cin >> iType; std::cin.ignore(); // skip over end-of-line character
		type = (iType == 1) ? Position::Long : Position::Short;
		cout << endl << "Please input the BBG code of counterparty of the Equity Forward: " << endl << "[Ans]: ";
		std::getline(std::cin, bbgCodeCounterparty);

		// ============================================================================================
		// *** Coding for loading market data from csv accordingly ***
		// ============================================================================================

		MarketData oMarketData;		// initiate MarketData object for loading market data of all equity
		EquityMarketData oInvestorMarketData = oMarketData.LoadEquity(bbgCodeInvestor);
		EquityMarketData oUlyMarketData = oMarketData.LoadEquity(bbgCodeUly);
		EquityMarketData oCounterpartyMarketData = oMarketData.LoadEquity(bbgCodeCounterparty);

		Real S0 = oUlyMarketData.Spot_;
		Real spot_dividend = oUlyMarketData.Dividend_;
		Rate r = oUlyMarketData.Drift_;
		Volatility sigma = oUlyMarketData.Volatility_;

		vector<Real> oInvestorCdsSpread;
		vector<Real> oCounterptyCdsSpread;
		vector<Period> tenors;

		tenors.push_back(1 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds1Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds1Yrs_);
		tenors.push_back(3 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds3Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds3Yrs_);
		tenors.push_back(5 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds5Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds5Yrs_);
		tenors.push_back(7 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds7Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds7Yrs_);
		tenors.push_back(10 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds10Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds10Yrs_);

		vector<Date> maturities;
		for (int i = 0; i < tenors.size(); i++)
			maturities.push_back(calendar.adjust(todaysDate + tenors[i], Following));

		// ====================================================================================

		// initiazlize investor stepwise hazard rate / default time object
		HazardRateCalibrate oInvestorHazardRate = HazardRateCalibrate(todaysDate);
		oInvestorHazardRate.SetYieldRate(0.000);
		oInvestorHazardRate.SetRecoveryRate(0.40);
		oInvestorHazardRate.SetCdsSpread(oInvestorCdsSpread);
		oInvestorHazardRate.SetCdsTenor(tenors);
		oInvestorHazardRate.PerformCalibration();
		DefaultTimeGenerator oInvestorDefaultTime(todaysDate);
		oInvestorDefaultTime.SetHazardRates(maturities, oInvestorHazardRate.oHazardRates);

		// initiazlize counterparty stepwise hazard rate / default time object
		HazardRateCalibrate oCounterptyHazardRate = HazardRateCalibrate(todaysDate);
		oCounterptyHazardRate.SetYieldRate(0.000);
		oCounterptyHazardRate.SetRecoveryRate(0.40);
		oCounterptyHazardRate.SetCdsSpread(oCounterptyCdsSpread);
		oCounterptyHazardRate.SetCdsTenor(tenors);
		oCounterptyHazardRate.PerformCalibration();
		DefaultTimeGenerator oCounterptyDefaultTime(todaysDate);
		oCounterptyDefaultTime.SetHazardRates(maturities, oCounterptyHazardRate.oHazardRates);

		// ==========================================================================
		// *** Hardcoded configuration of the equity forward ***
		// ==========================================================================

		int no_of_Date_path = 5000;
		int no_of_stock_path = 10;
		Natural settle = 2;
		Spread q = 0.0;

		Date valueDate = calendar.adjust(todaysDate + 2 * Days, Following);
		Date maturity = calendar.adjust(todaysDate + maturityYear * Years, Following);
		Real own_recovery = 0.4;
		Real counterparty_recovery = 0.4;
		Real coll_inv = 0;
		Real coll_cpt = 0;
		Real coll_inv_yield = 0;
		Real coll_cpt_yield = 0;

		// ==========================================================================

		//initial the class for Monte Carlo simulation of the forward and stock with the input above 
		MCgenerator f(todaysDate, maturity, valueDate, settle, calendar, S0, dayCounter, K, q, spot_dividend, r, coll_inv_yield, coll_cpt_yield, sigma, type);

		// *** calculate CVA/DVA ***
		//initial the class BVA_calculator for calculating the CVA and DVA with the class initialized above
		BVA_calculator newbva(&f, &oInvestorHazardRate, &oCounterptyHazardRate, no_of_Date_path, no_of_stock_path, coll_inv, coll_cpt);

		vector<map<string, Real>> CVA_DVA;

		std::map<string, Real>::iterator ot;
		std::vector<map<string, Real>>::iterator it;

		// ==================================================================================
		// !!! Swich statement do not allow initialization inside the switch bracket
		Real oBVA = 0.0;
		Real investorLGDC = 0.0;
		Real cptyLGDC = 0.0;

		vector<Real> discount_rate;
		for (int i = 0; i < tenors.size(); i++)
			discount_rate.push_back(5);		//set the discount rate as 5%

		//Setup the AT1P model calibration for investor and counterparty
		At1pCalibration investor = At1pCalibration(oInvestorCdsSpread, discount_rate, todaysDate, tenors);
		At1pCalibration counterparty = At1pCalibration(oCounterptyCdsSpread, discount_rate, todaysDate, tenors);
		BVA_calculator at1p_bva(&f, &investor, &counterparty, own_recovery, counterparty_recovery, coll_inv, coll_cpt, no_of_Date_path, no_of_stock_path);

		vector<map<string, Real>> AT1P_CVA_DVA;
		// ==================================================================================

		switch (iBVAModel) {
		case 1:
			// ==========================================================================
			// *** Stepwise Hazard Rate model without Collateralization***
			// ==========================================================================
			CVA_DVA = newbva.BVA_cal(maturity, 'n',0,0);
			for (it = CVA_DVA.begin(); it != CVA_DVA.end(); ++it) {
				for (ot = (*it).begin(); ot != (*it).end(); ++ot) {
					cout << (ot->first) << " " << (ot->second) << endl;
					oBVA += ot->second;
				}
			}
			cout << "BVA = " << oBVA << endl;
			break;

		case 2:
			// ==========================================================================
			// *** Stepwise Hazard Rate model with Collateralization***
			// ==========================================================================
			
			cout << endl << "Please input the collateral LGD of investor: (e.g. 0.4)" << endl << "[Ans]: ";
			cin >> investorLGDC;
			cout << endl << "Please input the collateral LGD of counterparty: (e.g. 0.4)" << endl << "[Ans]: ";
			cin >> cptyLGDC;
			
			newbva = BVA_calculator(&f, &oInvestorHazardRate, &oCounterptyHazardRate, no_of_Date_path, no_of_stock_path, investorLGDC, cptyLGDC);
			CVA_DVA = newbva.BVA_cal(maturity, 'y', investorLGDC, cptyLGDC);

			for (it = CVA_DVA.begin(); it != CVA_DVA.end(); ++it) {
				for (ot = (*it).begin(); ot != (*it).end(); ++ot) {
					cout << (ot->first) << " " << (ot->second) << endl;
					oBVA += ot->second;
				}
			}
			cout << "BVA = " << oBVA << endl;
			break;

		case 3:
			// ==========================================================================
			// *** AT1P Hazard Rate model without Collateralization***
			// ==========================================================================

			AT1P_CVA_DVA = at1p_bva.BVA_cal_AT1P(maturity, 'n',0,0);

			for (it = AT1P_CVA_DVA.begin(); it != AT1P_CVA_DVA.end(); ++it) {
				for (ot = (*it).begin(); ot != (*it).end(); ++ot) {
					cout << (ot->first) << " " << (ot->second) << endl;
					oBVA += ot->second;
				}
			}
			cout << "BVA = " << oBVA << endl;
			break;

		case 4:
			// ==========================================================================
			// *** AT1P Hazard Rate model with Collateralization***
			// ==========================================================================

			cout << endl << "Please input the collateral LGD of investor: (e.g. 0.4)" << endl << "[Ans]: ";
			cin >> investorLGDC;
			cout << endl << "Please input the collateral LGD of counterparty: (e.g. 0.4)" << endl << "[Ans]: ";
			cin >> cptyLGDC;

			at1p_bva = BVA_calculator(&f, &investor, &counterparty, own_recovery, counterparty_recovery, investorLGDC, cptyLGDC, no_of_Date_path, no_of_stock_path);
			AT1P_CVA_DVA = at1p_bva.BVA_cal_AT1P(maturity, 'y', investorLGDC, cptyLGDC);

			for (it = AT1P_CVA_DVA.begin(); it != AT1P_CVA_DVA.end(); ++it) {
				for (ot = (*it).begin(); ot != (*it).end(); ++ot) {
					cout << (ot->first) << " " << (ot->second) << endl;
					oBVA += ot->second;
				}
			}
			cout << "BVA = " << oBVA << endl;
			break;
		}
	}
	catch (const std::runtime_error& re)
	{
		std::cerr << "Runtime error: " << re.what() << std::endl;
	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error occurred: " << ex.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown failure occured. Possible memory corruption" << std::endl;
	}

}

void RunGuiBVAWWR(int iBVAModel) {

	try {

		// ==========================================
		// *** specify date/calendar config below ***
		// ==========================================

		Calendar calendar = TARGET();
		DayCounter dayCounter = Actual365Fixed();   //<== agreed that 365 basis to used
		Date todaysDate = Date::todaysDate();

		// ============================================================================================
		// *** Coding for user input ***
		// ============================================================================================
				
		Real K;
		Real maturityYear;
		int iType; Position::Type type;
		string bbgCodeUly;
		string bbgCodeInvestor;
		string bbgCodeCounterparty;
		Real iWwrModel; WWRHazardModel eWwrModel;
		Real wwrSensiB;		// Sensitivity b in the WWR model

		cout << endl << "Please input the BBG code of investor in the Equity Forward: " << endl << "[Ans]: ";
		std::getline(std::cin, bbgCodeInvestor);
		cout << endl << "Please input the BBG code of underlying of the Equity Forward: " << endl << "[Ans]: ";
		std::getline(std::cin, bbgCodeUly);
		cout << endl << "Please input the strike of the Equity Forward: " << endl << "[Ans]: ";
		cin >> K;
		cout << endl << "Please input the number of year to matiruty for the Equity Forward: (e.g. 1)" << endl << "[Ans]: ";
		cin >> maturityYear;
		cout << endl << "Please input the position type of the Equity Forward: (1=Long, 2=Short)" << endl << "[Ans]: ";
		cin >> iType; std::cin.ignore(); // skip over end-of-line character
		type = (iType == 1) ? Position::Long : Position::Short;
		cout << endl << "Please input the BBG code of counterparty in the Equity Forward: " << endl << "[Ans]: ";
		std::getline(std::cin, bbgCodeCounterparty);
		cout << endl << "Please select the the WWR model to use: (1=exponential, 2=linear)" << endl << "[Ans]: ";
		cin >> iWwrModel;
		eWwrModel = (iWwrModel == 1) ? EXPONENTIALLY : LINEARLY;
		cout << endl << "Please input the sensitivity b in the WWR model: " << endl << "[Ans]: ";
		cin >> wwrSensiB;

		// ============================================================================================
		// *** Coding for loading market data from csv accordingly ***
		// ============================================================================================

		MarketData oMarketData;		// initiate MarketData object for loading market data of all equity
		EquityMarketData oInvestorMarketData = oMarketData.LoadEquity(bbgCodeInvestor);
		EquityMarketData oUlyMarketData = oMarketData.LoadEquity(bbgCodeUly);
		EquityMarketData oCounterpartyMarketData = oMarketData.LoadEquity(bbgCodeCounterparty);

		Real S0 = oUlyMarketData.Spot_;
		Real spot_dividend = oUlyMarketData.Dividend_;
		Rate r = oUlyMarketData.Drift_;
		Volatility sigma = oUlyMarketData.Volatility_;

		vector<Real> oInvestorCdsSpread;
		vector<Real> oCounterptyCdsSpread;
		vector<Period> tenors;

		tenors.push_back(6 * Months);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds6Mths_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds6Mths_);
		tenors.push_back(1 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds1Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds1Yrs_);
		tenors.push_back(2 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds2Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds2Yrs_);
		tenors.push_back(3 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds3Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds3Yrs_);
		tenors.push_back(4 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds4Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds4Yrs_);
		tenors.push_back(5 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds5Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds5Yrs_);
		tenors.push_back(7 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds7Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds7Yrs_);
		tenors.push_back(10 * Years);
		oInvestorCdsSpread.push_back(oInvestorMarketData.Cds10Yrs_);
		oCounterptyCdsSpread.push_back(oCounterpartyMarketData.Cds10Yrs_);

		vector<Date> maturities;
		for (int i = 0; i < tenors.size(); i++) {
			maturities.push_back(calendar.adjust(todaysDate + tenors[i],
				Following));
		}

		// ============================================================================================
		// *** hard-code parameter specify here ***
		// ============================================================================================

		int no_of_Date_path = 1000;		// input parameter to BVACalculatorWWR 
		int no_of_stock_path = 10;		// input parameter to BVACalculatorWWR but actually not used in BVA WWR
		Natural settle = 2;
		Spread q = 0.0;
		Date valueDate = calendar.adjust(todaysDate + 2 * Days, Following);
		Date maturity = calendar.adjust(todaysDate + maturityYear * Years, Following);
		Real recoveryInvestor = 0.4;
		Real recoveryCpty = 0.4;

		// ============================================================================================

		//initial the class for Monte Carlo simulation of the forward and stock with the input above 
		MCgenerator oMCgenerator(todaysDate, maturity, valueDate, settle, calendar, S0, dayCounter, K, q, spot_dividend, r, sigma, type);

		// Specify forward payoff type for investor & counterparty
		Position::Type typeInvestorCalibrate;
		Position::Type typeCounterpartyCalibrate;

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

		cout << endl;
		HazardRateWWR oInvestorHazardRateWWR = HazardRateWWR(todaysDate);
		oInvestorHazardRateWWR.SetCdsPayDate(maturities);
		oInvestorHazardRateWWR.SetCdsSpread(oInvestorCdsSpread);
		oInvestorHazardRateWWR.SetRecoveryRate(recoveryInvestor);
		oInvestorHazardRateWWR.SetForwardType(typeInvestorCalibrate);
		oInvestorHazardRateWWR.SetSensitivityB(wwrSensiB);
		oInvestorHazardRateWWR.SetWWRHazardModel(eWwrModel);		// Select either exponential or linear for investor hazard rate function
		oInvestorHazardRateWWR.SimulateStockPrices(oMCgenerator);
		cout << "Start calibrating a(ti) of HazardRateWWR of investor " << endl;
		oInvestorHazardRateWWR.PerformCalibration();

		HazardRateWWR oCptyHazardRateWWR = HazardRateWWR(todaysDate);
		oCptyHazardRateWWR.SetCdsPayDate(maturities);
		oCptyHazardRateWWR.SetCdsSpread(oCounterptyCdsSpread);
		oCptyHazardRateWWR.SetRecoveryRate(recoveryCpty);
		oCptyHazardRateWWR.SetForwardType(typeCounterpartyCalibrate);
		oCptyHazardRateWWR.SetSensitivityB(wwrSensiB);
		oCptyHazardRateWWR.SetWWRHazardModel(eWwrModel);		// Select either exponential or linear for counterparty hazard rate function
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
		cout << "UCVA to investor perspective = " << oUCVA << endl;
		cout << "BVA = UDVA - UCVA = " << oUDVA - oUCVA << endl << endl;

	}
	catch (const std::runtime_error& re)
	{
		std::cerr << "Runtime error: " << re.what() << std::endl;
	}
	catch (const std::exception& ex)
	{
		std::cerr << "Error occurred: " << ex.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown failure occured. Possible memory corruption" << std::endl;
	}

}

void RunGui() {

	int iBVAModel;

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "+                                                                +" << endl;
	cout << "+ MAFS5220 Group 4: Bilateral Value Adjustment of Equity Forward +" << endl;
	cout << "+                                                                +" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

	cout << "Please select which model to use in calculating BVA: (No.1-4)" << endl;
	cout << "1) Stepwise Hazard Rate model without Collateralization" << endl;
	cout << "2) Stepwise Hazard Rate model with Collateralization (no re-hypothecation)" << endl;
	cout << "3) AT1P Hazard Rate model without Collateralization" << endl;
	cout << "4) AT1P Hazard Rate model with Collateralization (no re-hypothecation)" << endl;
	cout << "5) Hazard Rate with WWR" << endl << "[Ans]: ";
	cin >> iBVAModel; std::cin.ignore(); // skip over end-of-line character

	if (iBVAModel == 1 || iBVAModel == 2 || iBVAModel == 3 || iBVAModel == 4) {
		RunGuiBVA(iBVAModel);
	}
	else if (iBVAModel == 5) {
		RunGuiBVAWWR(iBVAModel);
	}
}


int main(int, char*[]) {

	//TestBVA();

	//CalculateBVAWWR();

	RunGui();

	system("pause");
}

