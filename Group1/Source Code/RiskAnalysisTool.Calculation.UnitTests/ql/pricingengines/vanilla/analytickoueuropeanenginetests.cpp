#include <pch.h>
#include <limits>
#include <ql/pricingengines/vanilla/analytickoueuropeanengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/exercise.hpp>
#include <ql/time/daycounters/business252.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(AnalyticKouEuropeanEngineTests)
	{
	private:
		boost::shared_ptr<KouProcess> CreateProcess(
			const Calendar &cal, const DayCounter &dc, const Date &d0,
			Real x0, Real q, Real r, Real vol, const Real jumpIntensity,
			const Real posProbability, const Real posJumpMean, const Real negJumpMean)
		{
			return boost::make_shared<KouProcess>(
				Handle<Quote>(boost::make_shared<SimpleQuote>(x0)),
				Handle<YieldTermStructure>(boost::make_shared<FlatForward>(d0, q, dc)),
				Handle<YieldTermStructure>(boost::make_shared<FlatForward>(d0, r, dc)),
				Handle<BlackVolTermStructure>(boost::make_shared<BlackConstantVol>(d0, cal, vol, dc)),
				jumpIntensity,
				posProbability,
				posJumpMean,
				negJumpMean
				);
		}
	public:

		TEST_METHOD(ConstructorTest)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 2;
			const Real negJumpMean = 1;
            const Real k = 80;

			std::vector<Real> seq1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
			std::vector<Real> seq2 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);

			Assert::IsNotNull(engine.get());
		}

		TEST_METHOD(CalculateTest_InTheMoneyCall)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = exp(0.5);
			const Real negJumpMean = exp(1 / 1.5);
			const Real k = 80;

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Call, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);

			Assert::IsNotNull(engine.get());

			option->setPricingEngine(engine);

			Real actual = option->NPV();
			Real expected = 20.31; 

			Assert::AreEqual(expected, actual, expected * 0.01);
		}

		TEST_METHOD(CalculateTest_OutOfTheMoneyPut)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;
			
			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 2;
			const Real negJumpMean = 1;
			const Real k = 80;

			std::vector<Real> seq1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
			std::vector<Real> seq2 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Put, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);

			Assert::IsNotNull(engine.get());

			option->setPricingEngine(engine);

			Real actual = option->NPV();
			Real expected = 3.85;    // Don't know what the value should be, but definitely not +infinity

			Assert::AreEqual(expected, actual, expected * 0.01);
		}

		TEST_METHOD(CalculateTest_OutOfTheMoneyCall)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 2;
			const Real negJumpMean = 1;
			const Real k = 120;

			std::vector<Real> seq1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
			std::vector<Real> seq2 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Call, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);

			Assert::IsNotNull(engine.get());

			option->setPricingEngine(engine);

			Real actual = option->NPV();
			Real expected = 5.96;

			Assert::AreEqual(expected, actual, expected * 0.01);
		}

		TEST_METHOD(CalculateTest_InTheMoneyPut)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 2;
			const Real negJumpMean = 1;
			const Real k = 120;

			std::vector<Real> seq1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
			std::vector<Real> seq2 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Put, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);

			Assert::IsNotNull(engine.get());

			option->setPricingEngine(engine);

			Real actual = option->NPV();
			Real expected = 27.29;

			Assert::AreEqual(expected, actual, expected * 0.01);
		}

		TEST_METHOD(CalculateTest_InTheMoneyCall_ReducedToBS)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.0;
			const Real posProb = 0.5;
			const Real posJumpMean = exp(0.5);
			const Real negJumpMean = exp(1 / 1.5);
			const Real k = 80;

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Call, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);
			option->setPricingEngine(engine);
			Real actual = option->NPV();

			auto engineBS = boost::make_shared<AnalyticEuropeanEngine>(process);
			option->setPricingEngine(engineBS);
			Real expected = option->NPV();

			Assert::AreEqual(expected, actual, expected * 0.01);
		}

		TEST_METHOD(CalculateTest_OutOfTheMoneyPut_ReducedToBS)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.0;
			const Real posProb = 0.5;
			const Real posJumpMean = exp(0.5);
			const Real negJumpMean = exp(1 / 1.5);
			const Real k = 80;

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Put, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);
			option->setPricingEngine(engine);
			Real actual = option->NPV();

			auto engineBS = boost::make_shared<AnalyticEuropeanEngine>(process);
			option->setPricingEngine(engineBS);
			Real expected = option->NPV();

			Assert::AreEqual(expected, actual, expected * 0.02);
		}

		TEST_METHOD(CalculateTest_OutOfTheMoneyCall_ReducedToBS)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.0;
			const Real posProb = 0.5;
			const Real posJumpMean = exp(0.5);
			const Real negJumpMean = exp(1 / 1.5);
			const Real k = 120;

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Call, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);
			option->setPricingEngine(engine);
			Real actual = option->NPV();

			auto engineBS = boost::make_shared<AnalyticEuropeanEngine>(process);
			option->setPricingEngine(engineBS);
			Real expected = option->NPV();

			Assert::AreEqual(expected, actual, expected * 0.02);
		}

		TEST_METHOD(CalculateTest_InTheMoneyPut_ReducedToBS)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.00;
			const Real posProb = 0.5;
			const Real posJumpMean = exp(0.5);
			const Real negJumpMean = exp(1 / 1.5);
			const Real k = 120;

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Put, k),
				boost::make_shared<EuropeanExercise>(Date(2, Jan, 2014))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);
			option->setPricingEngine(engine);
			Real actual = option->NPV();

			auto engineBS = boost::make_shared<AnalyticEuropeanEngine>(process);
			option->setPricingEngine(engineBS);
			Real expected = option->NPV();

			Assert::AreEqual(expected, actual, expected * 0.01);
		}


		TEST_METHOD(CalculateTest_From_Kou_Original_Paper)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			//dc = Actual365Fixed();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			Settings::instance().evaluationDate() = d0;

			const Real x0 = 100.0;
			const Real r = 0.05;
			const Real q = 0;
			const Real vol = 0.16;
			const Real jumpIntensity = 1;
			const Real posProb = 0.4;
			const Real posJumpMean = 10;
			const Real negJumpMean = 5;
			const Real k = 98.0;

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);
				
			auto option = boost::make_shared<VanillaOption>(
				boost::make_shared<PlainVanillaPayoff>(Option::Call, k),
				boost::make_shared<EuropeanExercise>(d0+Period(6,Months))
				);

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);
			option->setPricingEngine(engine);
			Real actual = option->NPV();


			Real expected = 9.14732; // Kou's paper

			Assert::AreEqual(expected, actual, expected * 0.02);
		}
	};
}