#include "stdafx.h"
#include "CppUnitTest.h"
#include "IntensityModel.hpp"

#include <cstring>
#include <iostream>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(IntensityModelTest)
	{
	private:
		IntensityModel *model;

	public:

		TEST_METHOD_INITIALIZE(setup)
		{
			Logger::WriteMessage("IntensityModelTest Setup\n");
			model = new IntensityModel();
		}

		TEST_METHOD_CLEANUP(tearDown)
		{
			Logger::WriteMessage("IntensityModelTest tearDown\n");
			delete model;
		}

		TEST_METHOD(IntensityLehmanDefaultTest1)
		{
			std::vector<Date> sampleDate;
			Date tmpDate(10, June, 2007);
			sampleDate.push_back(tmpDate);
			sampleDate.push_back(tmpDate + 1 * Years);
			sampleDate.push_back(tmpDate + 3 * Years);
			sampleDate.push_back(tmpDate + 5 * Years);
			sampleDate.push_back(tmpDate + 7 * Years);
			sampleDate.push_back(tmpDate + 10 * Years);

			std::vector<Rate> sampleIntensity;
			sampleIntensity.push_back(0.00267);
			sampleIntensity.push_back(0.00601);
			sampleIntensity.push_back(0.01217);
			sampleIntensity.push_back(0.01096);
			sampleIntensity.push_back(0.01407);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(1.0);
			expectedProb.push_back(0.997);
			expectedProb.push_back(0.985);
			expectedProb.push_back(0.962);
			expectedProb.push_back(0.941);
			expectedProb.push_back(0.902);

			model->Initialization(sampleDate, sampleIntensity);

			for (std::size_t i = 0; i < expectedProb.size(); i++) {
				Probability resultedProb = model->SurvivalProbability(sampleDate[i]);
				std::string msg = "[Test" + std::to_string(i) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}

		TEST_METHOD(IntensityLehmanDefaultTest2)
		{
			std::vector<Date> sampleDate;
			Date tmpDate(10, June, 2008);
			sampleDate.push_back(tmpDate);
			sampleDate.push_back(tmpDate + 1 * Years);
			sampleDate.push_back(tmpDate + 3 * Years);
			sampleDate.push_back(tmpDate + 5 * Years);
			sampleDate.push_back(tmpDate + 7 * Years);
			sampleDate.push_back(tmpDate + 10 * Years);

			std::vector<Rate> sampleIntensity;
			sampleIntensity.push_back(0.06563);
			sampleIntensity.push_back(0.0444);
			sampleIntensity.push_back(0.03411);
			sampleIntensity.push_back(0.03207);
			sampleIntensity.push_back(0.02907);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(1.0);
			expectedProb.push_back(0.936);
			expectedProb.push_back(0.857);
			expectedProb.push_back(0.800);
			expectedProb.push_back(0.751);
			expectedProb.push_back(0.688);

			model->Initialization(sampleDate, sampleIntensity);

			for (std::size_t i = 0; i < expectedProb.size(); i++) {
				Probability resultedProb = model->SurvivalProbability(sampleDate[i]);
				std::string msg = "[Test" + std::to_string(i) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}

		TEST_METHOD(IntensityLehmanDefaultTest3)
		{
			std::vector<Date> sampleDate;
			Date tmpDate(12, Sep, 2008);
			sampleDate.push_back(tmpDate);
			sampleDate.push_back(tmpDate + 1 * Years);
			sampleDate.push_back(tmpDate + 3 * Years);
			sampleDate.push_back(tmpDate + 5 * Years);
			sampleDate.push_back(tmpDate + 7 * Years);
			sampleDate.push_back(tmpDate + 10 * Years);

			std::vector<Rate> sampleIntensity;
			sampleIntensity.push_back(0.2326);
			sampleIntensity.push_back(0.09248);
			sampleIntensity.push_back(0.05245);
			sampleIntensity.push_back(0.05947);
			sampleIntensity.push_back(0.06422);

			std::vector<Probability> expectedProb;
			expectedProb.push_back(1.0);
			expectedProb.push_back(0.792);
			expectedProb.push_back(0.659);
			expectedProb.push_back(0.593);
			expectedProb.push_back(0.527);
			expectedProb.push_back(0.434);

			model->Initialization(sampleDate, sampleIntensity);

			for (std::size_t i = 0; i < expectedProb.size(); i++) {
				Probability resultedProb = model->SurvivalProbability(sampleDate[i]);
				std::string msg = "[Test" + std::to_string(i) + "] Expected Prob is "
					+ std::to_string(expectedProb[i]) + ", Real Prob is " +
					std::to_string(resultedProb) + "\n";
				Logger::WriteMessage(msg.c_str());
				Assert::AreEqual(resultedProb, expectedProb[i], 0.01);
			}
		}
	};
}