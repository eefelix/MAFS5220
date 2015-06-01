#include "stdafx.h"
#include "CppUnitTest.h"
#include "IRSUCVAWithWWR.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Quantlib{
	TEST_CLASS(IRSUCVAWithWWRTest)
	{
	private:

	public:

		TEST_METHOD(wwr_test)
		{
			for (int i = 1; i <= 10; i++){
				boost::shared_ptr<CVAwithWWR> calculator = boost::shared_ptr<CVAwithWWR>(new CVAwithWWR());
				calculator->WithShortRateCIRParam(0.04, 0.0404, 0.1, 0.2);
				calculator->WithIntensityCIRParam(0.02*i, 0.1, 0.05, 0.2);
				calculator->WithRecoveryRate(0.4);
				Date startDate(31, December, 2012);
				Period tenor(5, Years);
				//Real nominal, Real spread, Real fixedRate, Date startDate, Frequency fixedFrequency,
				//Frequency floatingFrequency, Period tenor, Type type
				calculator->WithIRSSwap(1.0, 0.0001, 0.0405, startDate, Annual,
					Semiannual, tenor, CVAwithWWR::Payer);
				calculator->WithReferenceDate(startDate);
				calculator->WithCIRCorrelation(-0.3);
				Real count = calculator->cvaWithWWRCalculation(10000);
				std::string msg = "UCVA with WWR result is " + std::to_string(count) + ".\n";
				Logger::WriteMessage(msg.c_str());
			}
		}
	};
}