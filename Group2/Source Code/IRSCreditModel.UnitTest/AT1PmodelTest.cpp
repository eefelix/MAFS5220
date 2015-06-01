#include "stdafx.h"
#include "CppUnitTest.h"
#include "AT1Pmodel.hpp"

#include <cstring>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{		
	TEST_CLASS(AT1PmodelTest)
	{
	private:
		AT1Pmodel *at1p;

	public:

		TEST_METHOD_INITIALIZE(setup)
		{
			Logger::WriteMessage("AT1PTest Setup\n");
			at1p = new AT1Pmodel();
		}

		TEST_METHOD_CLEANUP(tearDown)
		{
			Logger::WriteMessage("AT1PTest tearDown\n");
			delete at1p;
		}
		
		TEST_METHOD(AT1PLehmanDefaultTest1)
		{
			std::vector<Date> sampleDate;
			Date tmpDate(10, June, 2007);
			sampleDate.push_back(tmpDate);
			sampleDate.push_back(tmpDate + 1 * Years);
			sampleDate.push_back(tmpDate + 3 * Years);
			sampleDate.push_back(tmpDate + 5 * Years);
			sampleDate.push_back(tmpDate + 7 * Years);
			sampleDate.push_back(tmpDate + 10 * Years);

			std::vector<Volatility> sampleVol;
			sampleVol.push_back(0.292);
			sampleVol.push_back(0.140);
			sampleVol.push_back(0.145);
			sampleVol.push_back(0.120);
			sampleVol.push_back(0.127);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(1.0);
			expectedProb.push_back(0.997);
			expectedProb.push_back(0.985);
			expectedProb.push_back(0.961);
			expectedProb.push_back(0.941);
			expectedProb.push_back(0.902);
			
			at1p->Initialization(sampleDate, sampleVol);
			
			for (std::size_t i = 0; i < expectedProb.size(); i++) {
				Probability resultedProb = at1p->SurvivalProbability(sampleDate[i]);
				std::string msg = "[Test" + std::to_string(i) + "] Expected Prob is " 
					+ std::to_string(expectedProb[i]) + ", Real Prob is " + 
					std::to_string(resultedProb) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}

		TEST_METHOD(AT1PLehmanDefaultTest2)
		{
			std::vector<Date> sampleDate;
			Date tmpDate(10, June, 2008);
			sampleDate.push_back(tmpDate);
			sampleDate.push_back(tmpDate + 1 * Years);
			sampleDate.push_back(tmpDate + 3 * Years);
			sampleDate.push_back(tmpDate + 5 * Years);
			sampleDate.push_back(tmpDate + 7 * Years);
			sampleDate.push_back(tmpDate + 10 * Years);

			std::vector<Volatility> sampleVol;
			sampleVol.push_back(0.450);
			sampleVol.push_back(0.219);
			sampleVol.push_back(0.186);
			sampleVol.push_back(0.181);
			sampleVol.push_back(0.175);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(1.0);
			expectedProb.push_back(0.935);
			expectedProb.push_back(0.856);
			expectedProb.push_back(0.799);
			expectedProb.push_back(0.750);
			expectedProb.push_back(0.687);

			at1p->Initialization(sampleDate, sampleVol);

			for (std::size_t i = 0; i < expectedProb.size(); i++) {
				Probability resultedProb = at1p->SurvivalProbability(sampleDate[i]);
				std::string msg = "[Test" + std::to_string(i) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}

		TEST_METHOD(AT1PLehmanDefaultTest3)
		{
			std::vector<Date> sampleDate;
			Date tmpDate(12, Sep, 2008);
			sampleDate.push_back(tmpDate);
			sampleDate.push_back(tmpDate + 1 * Years);
			sampleDate.push_back(tmpDate + 3 * Years);
			sampleDate.push_back(tmpDate + 5 * Years);
			sampleDate.push_back(tmpDate + 7 * Years);
			sampleDate.push_back(tmpDate + 10 * Years);

			std::vector<Volatility> sampleVol;
			sampleVol.push_back(0.622);
			sampleVol.push_back(0.308);
			sampleVol.push_back(0.243);
			sampleVol.push_back(0.269);
			sampleVol.push_back(0.295);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(1.0);
			expectedProb.push_back(0.784);
			expectedProb.push_back(0.655);
			expectedProb.push_back(0.591);
			expectedProb.push_back(0.525);
			expectedProb.push_back(0.434);

			at1p->Initialization(sampleDate, sampleVol);

			for (std::size_t i = 0; i < expectedProb.size(); i++) {
				Probability resultedProb = at1p->SurvivalProbability(sampleDate[i]);
				std::string msg = "[Test" + std::to_string(i) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}
	};
}