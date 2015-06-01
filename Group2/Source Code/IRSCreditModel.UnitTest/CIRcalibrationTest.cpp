#include "stdafx.h"
#include <fstream>
#include "CppUnitTest.h"
#include "CIRcalibration.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(CDScalibrationTest)
	{
	private:
		CIRcalibration* calibrator;

	public:
		TEST_METHOD_INITIALIZE(setup)
		{
			Logger::WriteMessage("CIRcalibrationTest Setup\n");
			calibrator = new CIRcalibration();
		}
		TEST_METHOD_CLEANUP(tearDown)
		{
			Logger::WriteMessage("CIRcalibrationTest tearDown\n");
			delete calibrator;
		}

		TEST_METHOD(CIRCalculationTest_Manmade)
		{
			std::vector<Real> sampleQuote;
			for (double j = 0; j != 6; j++){
				sampleQuote.push_back((j + 1) / 100);
			}

			calibrator->Initialization(sampleQuote);
			calibrator->StartCalibrate();
			calibrator->CalculateParameters();

			std::string msg1 = "[CIRCalculationTest_Manmade: Mean Level] Value is " 
				+ std::to_string(calibrator->getFinalizedCenter()) + "\n";
			Logger::WriteMessage(msg1.c_str());

			std::string msg2 = "[CIRCalculationTest_Manmade: Volatility] Value is " 
				+ std::to_string(calibrator->getFinalizedVolatilty()) + "\n";
			Logger::WriteMessage(msg2.c_str());

			std::string msg3 = "[CIRCalculationTest_Manmade: Drift] Value is "
				+ std::to_string(calibrator->getFinalizedDrift()) + "\n";
			Logger::WriteMessage(msg3.c_str());
		}

		/**
		* This test method use US daily Treasury Yield Curve data from 2-Jan-2009 to 27-Jul-2010.
		* Data source: http://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield
		* We use 1-month yield data as proxy for short rate, which is also 1-month zero rate.
		* We do the calibration in excel in comparison to our result in IRSCreditModel. Calibration Method refers to:
		* http://financetrainingcourse.com/education/2012/06/cox-ingersoll-ross-cir-interest-rate-model-parameter-calibration-short-rates-simulation-and-modeling-of-longer-term-interest-rates-an-example/
		*/
		TEST_METHOD(CIRCalibrationTest_USTreasury)
		{
			Real expectedPhi = 0.9921707;
			Real expectedGamma = 0.0009959;
			Real expectedMinimizedRSS = 0.0205703;
			Real expectedKappa = 0.0078602;
			Real expectedContinuousSigma = 0.0072818;

			std::vector<Real> ShortRate1M;
			Real tmpValue;
			std::ifstream infile("..\\IRSCreditModel.UnitTest\\CIRCalibrationTestData_USTreasury.txt");
			if (!infile.is_open()) {
				Logger::WriteMessage("Open file error!");
				Assert::Fail();
			}
			while (infile >> tmpValue) {
				ShortRate1M.push_back(tmpValue);
			}
			infile.close();
			
			calibrator->Initialization(ShortRate1M);
			calibrator->StartCalibrate();
			calibrator->CalculateParameters();

			std::string msg1 = "[CIRCalibrationTest_USTreasury: Mean level] Value is " 
				+ std::to_string(calibrator->getFinalizedCenter()) + ", expected value is " + std::to_string(expectedGamma) + "\n";
			Logger::WriteMessage(msg1.c_str());
			Assert::AreEqual(calibrator->getFinalizedCenter(), expectedGamma, 0.01);
			std::string msg2 = "[CIRCalibrationTest_USTreasury: Volatility] Value is " 
				+ std::to_string(calibrator->getFinalizedVolatilty()) + ", expected value is " + std::to_string(expectedContinuousSigma) + "\n";
			Logger::WriteMessage(msg2.c_str());
			Assert::AreEqual(calibrator->getFinalizedVolatilty(), expectedContinuousSigma, 0.01);
			std::string msg3 = "[CIRCalibrationTest_USTreasury: Drift] Value is " 
				+ std::to_string(calibrator->getFinalizedDrift()) + ", expected value is " + std::to_string(expectedKappa) + "\n";
			Logger::WriteMessage(msg3.c_str());
			Assert::AreEqual(calibrator->getFinalizedDrift(), expectedKappa, 0.01);
		}

	};
}