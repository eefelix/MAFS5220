#include "stdafx.h"
#include "CppUnitTest.h"
#include "LeastSquareMC.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Quantlib{
	TEST_CLASS(LeastSquareMCTest)
	{
	private:

	public:
		/**
		* This test method examine the correctness of American monte carlo calculation
		* We adopt artificial American option 
		* But compare with our result with the benchmark from
		* http://www.math.columbia.edu/~smirnov/options13s.html
		*/
		TEST_METHOD(LSMC_test)
		{
				LeastSquareMC mc(10000, 20, 100.0, 1.0, 95.0, 0.1, 0.3, LeastSquareMC::Put);
				double price = mc.calcLeastSquareMC();
				std::string msg = "Price of the american put option is " + std::to_string(price) + ".\n";
				Logger::WriteMessage(msg.c_str());

		}
	};
}