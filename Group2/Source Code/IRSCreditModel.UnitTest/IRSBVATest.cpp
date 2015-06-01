#include "stdafx.h"
#include "CppUnitTest.h"
#include "IRSBVA.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Quantlib{
	TEST_CLASS(IRSBVATest)
	{
	private:

	public:

		TEST_METHOD(irsbva_test)
		{
			boost::shared_ptr<IRSBVA> calculator = boost::shared_ptr<IRSBVA>(new IRSBVA());
			calculator->WithShortRateCIRParam(0.05, 0.1, 0.1, 0.1);
			calculator->WithBIntensityCIRParam(0.0165, 0.4, 0.026, 0.14);
			calculator->WithCIntensityCIRParam(0.0165, 0.4, 0.026, 0.14);
			calculator->WithCIRCorrelation(0.1,0.2,0.3);
			calculator->WithRecoveryRate(0.4,0.4);
			Date startDate(31, December, 2012);
			Period tenor(5, Years);
			calculator->WithIRSSwap(1.0, 0.02, 0.04, startDate, Annual,
				Quarterly, tenor, IRSBVA::Payer);
			calculator->WithReferenceDate(startDate);

			Real count = calculator->bvaCalculation(10000);
			std::string msg = "BVA result is " + std::to_string(count) + ".\n";
			Logger::WriteMessage(msg.c_str());
		}
	};
}