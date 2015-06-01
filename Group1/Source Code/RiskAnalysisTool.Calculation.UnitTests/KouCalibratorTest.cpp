#include <pch.h>
#include <limits>
#include <ql/pricingengines/vanilla/analytickoueuropeanengine.hpp>
#include <ql/calibration/kouprocesscalibrator.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/business252.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/exercise.hpp>
#include <ql/option.hpp>
#include <ql/instruments/vanillaoption.hpp>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(KouCalibratorTests)
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

		TEST_METHOD(Calibrate)
		{
			Calendar cal;
			DayCounter dc;
			Date d0,d1;


			cal = TARGET();
			dc = Business252();
			d0 = Date(22, Feb, 2015);
			d1 = Date(22, Jan, 2016);

			// calibrate Kou Process
			std::shared_ptr<std::vector<double>> strike(new std::vector<double>);
			std::shared_ptr<std::vector<Date>> maturity(new std::vector<Date>);
			std::shared_ptr<std::vector<double>> optionPrice(new std::vector<double>);
			std::shared_ptr<QuantLib::EndCriteria> endCriteria(new QuantLib::EndCriteria(10000, 100, 1e-5, 1e-5, 1e-5));
			double price[5] = { 34.75, 31.5, 30.0, 25.85, 21.00 };
			double s[5] = { 95.71, 100.00, 104.29, 110.0, 114.29 };

			for (int i = 0; i < 5; ++i){
				strike->push_back(s[i]);
				maturity->push_back(d1);
				optionPrice->push_back(price[i]);
			}

			auto calibrator
				= boost::make_shared<KouProcessCalibrator>(d0, cal, dc, 0.05, 129.49, 0.01, strike, maturity, optionPrice, endCriteria);

			calibrator->calibrate();

			// recover option price, here maturity is not exactly the same
			Calendar cal1;
			DayCounter dc1;
			Date d2;


			cal1 = TARGET();
			dc1 = Business252();
			d2 = Date(22, Feb, 2015);

			Settings::instance().evaluationDate() = d2;

			const Real x0 = 129.49;
			const Real r = 0.05;
			const Real q = 0.01;
			const Real vol = calibrator->getVol().blackVol(0,0,true);
			const Real jumpIntensity = calibrator->getJumpIntensity();
			const Real posProb = calibrator->getPosProb();
			const Real posJumpMean = calibrator->getPosJumpMean();
			const Real negJumpMean = calibrator->getNegJumpMean();

			auto process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean
				);

			std::vector<boost::shared_ptr<QuantLib::VanillaOption>> option;

			for (int i = 0; i < 5; ++i)
				option.push_back(boost::make_shared<VanillaOption>(
					boost::make_shared<PlainVanillaPayoff>(Option::Call, s[i]),
					boost::make_shared<EuropeanExercise>(Date(22, Jan, 2016))
					));

			auto engine = boost::make_shared<AnalyticKouEuropeanEngine>(process);

			Assert::IsNotNull(engine.get());

			for (int i = 0; i < option.size(); ++i)
				option[i]->setPricingEngine(engine);

			// price with calib parameters and price difference(relative)
			std::vector<double> calibPrice(5), priceDiff(5);
			for (int i = 0; i < calibPrice.size(); ++i){
				calibPrice[i] = option[i]->NPV();
				priceDiff[i] = (calibPrice[i] - price[i]) / price[i];
			}

			for (int i = 0; i < priceDiff.size(); ++i)
				Assert::IsTrue(priceDiff[i] < 0.1);
		}
	};
}
