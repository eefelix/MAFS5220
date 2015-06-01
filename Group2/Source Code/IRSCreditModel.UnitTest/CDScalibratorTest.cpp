#include "stdafx.h"
#include "CppUnitTest.h"
#include "CDScalibrator.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(CDScalibratorTest)
	{
	private:
		Real RecoveryRate = 0.4;
		CDSCalibrator* calibrator;

	public:
		TEST_METHOD_INITIALIZE(setup)
		{
			Logger::WriteMessage("CDScalibratorTest Setup\n");
			calibrator = new CDSCalibrator(RecoveryRate);
		}
		TEST_METHOD_CLEANUP(tearDown)
		{
			Logger::WriteMessage("CDScalibratorTest tearDown\n");
			delete calibrator;
		}

		TEST_METHOD(CDSLehmanDefaultTestLowRisk)
		{
			Date todayDate(10, June, 2007);
			std::vector<Real> sampleQuote;
			sampleQuote.push_back(0.0016);
			sampleQuote.push_back(0.0029);
			sampleQuote.push_back(0.0045);
			sampleQuote.push_back(0.0050);
			sampleQuote.push_back(0.0058);

			std::vector<Period> sampleTenor;
			sampleTenor.push_back(1 * Years);
			sampleTenor.push_back(3 * Years);
			sampleTenor.push_back(5 * Years);
			sampleTenor.push_back(7 * Years);
			sampleTenor.push_back(10 * Years);

			std::vector<Volatility> expectedVol;
			expectedVol.push_back(0.292);
			expectedVol.push_back(0.140);
			expectedVol.push_back(0.145);
			expectedVol.push_back(0.120);
			expectedVol.push_back(0.127);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(0.997);
			expectedProb.push_back(0.985);
			expectedProb.push_back(0.961);
			expectedProb.push_back(0.941);
			expectedProb.push_back(0.902);

			calibrator->setStartDate(todayDate);
			calibrator->setQuotes(sampleQuote);
			calibrator->setTermStructure(sampleTenor);
			calibrator->extractDefaultProbabilityCurve();
			calibrator->StartCalibrate();

			Logger::WriteMessage("Unit test case for Lehman Default under low risk.");
			for (std::size_t i = 0; i < sampleQuote.size(); i++) {
				Probability resultedProb = calibrator->FinalizedSurvivalProbability(todayDate + sampleTenor[i]);
				Volatility resultedVol = calibrator->FinalizedVolatility(todayDate + sampleTenor[i]);
				std::string msg = "[Year " + std::to_string(sampleTenor[i].length()) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n" + "Expected Vol is " +
					std::to_string(expectedVol[i]) + ", Real Vol is " + std::to_string(resultedVol) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
			
		}

		TEST_METHOD(CDSLehmanDefaultTestMediumRisk)
		{
			Date todayDate(12, June, 2008);
			std::vector<Real> sampleQuote;
			sampleQuote.push_back(0.0397);
			sampleQuote.push_back(0.0315);
			sampleQuote.push_back(0.0277);
			sampleQuote.push_back(0.0258);
			sampleQuote.push_back(0.0240);

			std::vector<Period> sampleTenor;
			sampleTenor.push_back(1 * Years);
			sampleTenor.push_back(3 * Years);
			sampleTenor.push_back(5 * Years);
			sampleTenor.push_back(7 * Years);
			sampleTenor.push_back(10 * Years);

			std::vector<Volatility> expectedVol;
			expectedVol.push_back(0.450);
			expectedVol.push_back(0.219);
			expectedVol.push_back(0.186);
			expectedVol.push_back(0.181);
			expectedVol.push_back(0.175);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(0.935);
			expectedProb.push_back(0.856);
			expectedProb.push_back(0.799);
			expectedProb.push_back(0.750);
			expectedProb.push_back(0.687);

			calibrator->setStartDate(todayDate);
			calibrator->setQuotes(sampleQuote);
			calibrator->setTermStructure(sampleTenor);
			calibrator->extractDefaultProbabilityCurve();
			calibrator->StartCalibrate();

			Logger::WriteMessage("Unit test case for Lehman Default under medium risk.");
			for (std::size_t i = 0; i < sampleQuote.size(); i++) {
				Probability resultedProb = calibrator->FinalizedSurvivalProbability(todayDate + sampleTenor[i]);
				Volatility resultedVol = calibrator->FinalizedVolatility(todayDate + sampleTenor[i]);
				std::string msg = "[Year " + std::to_string(sampleTenor[i].length()) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n" + "Expected Vol is " +
					std::to_string(expectedVol[i]) + ", Real Vol is " + std::to_string(resultedVol) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}

		}

		TEST_METHOD(CDSLehmanDefaultTestHighRisk)
		{
			Date todayDate(12, Sep, 2008);
			std::vector<Real> sampleQuote;
			sampleQuote.push_back(0.1437);
			sampleQuote.push_back(0.0902);
			sampleQuote.push_back(0.0710);
			sampleQuote.push_back(0.0636);
			sampleQuote.push_back(0.0588);

			std::vector<Period> sampleTenor;
			sampleTenor.push_back(1 * Years);
			sampleTenor.push_back(3 * Years);
			sampleTenor.push_back(5 * Years);
			sampleTenor.push_back(7 * Years);
			sampleTenor.push_back(10 * Years);
			
			std::vector<Volatility> expectedVol;
			expectedVol.push_back(0.622);
			expectedVol.push_back(0.308);
			expectedVol.push_back(0.243);
			expectedVol.push_back(0.269);
			expectedVol.push_back(0.295);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(0.784);
			expectedProb.push_back(0.655);
			expectedProb.push_back(0.591);
			expectedProb.push_back(0.525);
			expectedProb.push_back(0.434);

			calibrator->setStartDate(todayDate);
			calibrator->setQuotes(sampleQuote);
			calibrator->setTermStructure(sampleTenor);
			calibrator->extractDefaultProbabilityCurve();
			calibrator->StartCalibrate();
			
			Logger::WriteMessage("Unit test case for Lehman Default under high risk.");
			for (std::size_t i = 0; i < sampleQuote.size(); i++) {
				Probability resultedProb = calibrator->FinalizedSurvivalProbability(todayDate + sampleTenor[i]);
				Volatility resultedVol = calibrator->FinalizedVolatility(todayDate + sampleTenor[i]);
				std::string msg = "[Year " + std::to_string(sampleTenor[i].length()) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n" + "Expected Vol is " +
					std::to_string(expectedVol[i]) + ", Real Vol is " + std::to_string(resultedVol) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}

	};
}