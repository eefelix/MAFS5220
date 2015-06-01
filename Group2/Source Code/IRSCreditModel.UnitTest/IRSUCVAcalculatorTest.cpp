#include "stdafx.h"
#include "CppUnitTest.h"
#include "IRSUCVAcalculator.hpp"
#include "DefaultCurve.hpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace Quantlib{

	TEST_CLASS(IRSUCVAcalculatorTest)
	{
	private:
		boost::shared_ptr<IRSUCVAcalculator> calculator;
		std::vector<Date> sampleDate;
		std::vector<Period> contractMaturity;
		std::vector<Rate> contractSwapRate;

	public:

		TEST_METHOD_INITIALIZE(setup)
		{
			Logger::WriteMessage("IRSUCVAcalculatorTest Setup\n");
			calculator = boost::shared_ptr<IRSUCVAcalculator>(new IRSUCVAcalculator());
			Date Date1(12, Mar, 2004);
			sampleDate.push_back(Date1);
			Date Date2(12, Mar, 2005);
			sampleDate.push_back(Date2);
			sampleDate.push_back(Date2 + 2 * Years);
			sampleDate.push_back(Date2 + 4 * Years);
			sampleDate.push_back(Date2 + 6 * Years);
			sampleDate.push_back(Date2 + 9 * Years);
			sampleDate.push_back(Date2 + 14 * Years);
			sampleDate.push_back(Date2 + 19 * Years);
			sampleDate.push_back(Date2 + 24 * Years);
			sampleDate.push_back(Date2 + 29 * Years);

			contractMaturity.push_back(Period(5, Years));
			contractMaturity.push_back(Period(10, Years));
			contractMaturity.push_back(Period(15, Years));
			contractMaturity.push_back(Period(20, Years));
			contractMaturity.push_back(Period(25, Years));
			contractMaturity.push_back(Period(30, Years));

			contractSwapRate.push_back(0.03248);
			contractSwapRate.push_back(0.04075);
			contractSwapRate.push_back(0.04462);
			contractSwapRate.push_back(0.04676);
			contractSwapRate.push_back(0.04775);
			contractSwapRate.push_back(0.04810);
		}

		TEST_METHOD_CLEANUP(tearDown)
		{
			Logger::WriteMessage("IRSUCVAcalculatorTest tearDown\n");
			sampleDate.clear();
			contractMaturity.clear();
			contractSwapRate.clear();
		}

		TEST_METHOD(SWAP_EXAMPLE)
		{
			std::vector<Date> dates;
			std::vector<DiscountFactor> discountFactor;

			Date valuationDate(31, December, 2012);
			Settings::instance().evaluationDate() = valuationDate;

			dates.push_back(valuationDate); discountFactor.push_back(1.0);
			dates.push_back(Date(31, December, 2013));  discountFactor.push_back(0.99);
			dates.push_back(Date(31, December, 2024));  discountFactor.push_back(0.80);
			boost::shared_ptr<YieldTermStructure> forwardCurve(
				new InterpolatedDiscountCurve<LogLinear>(dates, discountFactor, Actual360()));

			discountFactor.pop_back(); discountFactor.pop_back();

			discountFactor.push_back(0.999);
			discountFactor.push_back(0.89);

			boost::shared_ptr<YieldTermStructure> oisCurve(
				new InterpolatedDiscountCurve<LogLinear>(dates, discountFactor, Actual360()));

			Handle<YieldTermStructure> discountingTermStructure(oisCurve);
			Handle<YieldTermStructure> forwardingTermStructure(forwardCurve);

			Real nominal = 1000000.0;
			Date previousResetDate(20, November, 2012);
			Date maturity(20, November, 2022);
			double spread = 0.02;
			double fixedRate = 0.04;

			boost::shared_ptr<IborIndex> euribor(new Euribor(3 * Months, forwardingTermStructure));
			euribor->addFixing(euribor->fixingDate(previousResetDate), 0.01, true);

			VanillaSwap::Type swapType = VanillaSwap::Payer;

			Schedule fixedSchedule(previousResetDate, maturity, 1 * Years,
				TARGET(), ModifiedFollowing, ModifiedFollowing,
				DateGeneration::Backward, false);

			Schedule floatSchedule(previousResetDate, maturity, 3 * Months,
				TARGET(), ModifiedFollowing, ModifiedFollowing,
				DateGeneration::Backward, false);

			VanillaSwap swap(VanillaSwap::Payer, nominal, fixedSchedule, fixedRate, Thirty360(),
				floatSchedule, euribor, spread, Actual360());

			boost::shared_ptr<PricingEngine> swapEngine(new DiscountingSwapEngine(discountingTermStructure));

			swap.setPricingEngine(swapEngine);

			double res = swap.NPV();
			std::string msg = "[swap value is " + std::to_string(res) + "\n";
			Logger::WriteMessage(msg.c_str());
		}

		TEST_METHOD(UCVA_Paper)
		{
			std::vector<Date> zeroCouponRateDate;
			std::vector<Rate> zeroCouponRate;
			zeroCouponRateDate.push_back(Date(23, Jun, 2006)); zeroCouponRate.push_back(0.0001);
			zeroCouponRateDate.push_back(Date(26, Jun, 2006)); zeroCouponRate.push_back(0.0283);
			zeroCouponRateDate.push_back(Date(27, Jun, 2006)); zeroCouponRate.push_back(0.0283);
			zeroCouponRateDate.push_back(Date(28, Jun, 2006)); zeroCouponRate.push_back(0.0283);
			zeroCouponRateDate.push_back(Date(4, Jul, 2006)); zeroCouponRate.push_back(0.0287);
			zeroCouponRateDate.push_back(Date(11, Jul, 2006)); zeroCouponRate.push_back(0.0287);
			zeroCouponRateDate.push_back(Date(18, Jul, 2006)); zeroCouponRate.push_back(0.0287);
			zeroCouponRateDate.push_back(Date(27, Jul, 2006)); zeroCouponRate.push_back(0.0288);
			zeroCouponRateDate.push_back(Date(28, Aug, 2006)); zeroCouponRate.push_back(0.0292);
			zeroCouponRateDate.push_back(Date(20, Sep, 2006)); zeroCouponRate.push_back(0.0296);
			zeroCouponRateDate.push_back(Date(20, Dec, 2006)); zeroCouponRate.push_back(0.0314);
			zeroCouponRateDate.push_back(Date(20, Mar, 2007)); zeroCouponRate.push_back(0.0327);
			zeroCouponRateDate.push_back(Date(21, Jun, 2007)); zeroCouponRate.push_back(0.0338);
			zeroCouponRateDate.push_back(Date(20, Sep, 2007)); zeroCouponRate.push_back(0.0346);
			zeroCouponRateDate.push_back(Date(19, Dec, 2007)); zeroCouponRate.push_back(0.0352);
			zeroCouponRateDate.push_back(Date(19, Mar, 2008)); zeroCouponRate.push_back(0.0357);
			zeroCouponRateDate.push_back(Date(19, Jun, 2008)); zeroCouponRate.push_back(0.0361);
			zeroCouponRateDate.push_back(Date(18, Sep, 2008)); zeroCouponRate.push_back(0.0365);
			zeroCouponRateDate.push_back(Date(29, Jun, 2009)); zeroCouponRate.push_back(0.0375);
			zeroCouponRateDate.push_back(Date(28, Jun, 2010)); zeroCouponRate.push_back(0.0384);
			zeroCouponRateDate.push_back(Date(27, Jun, 2011)); zeroCouponRate.push_back(0.0391);
			zeroCouponRateDate.push_back(Date(27, Jun, 2012)); zeroCouponRate.push_back(0.0398);
			zeroCouponRateDate.push_back(Date(27, Jun, 2013)); zeroCouponRate.push_back(0.0403);
			zeroCouponRateDate.push_back(Date(27, Jun, 2014)); zeroCouponRate.push_back(0.0409);
			zeroCouponRateDate.push_back(Date(29, Jun, 2015)); zeroCouponRate.push_back(0.0414);
			zeroCouponRateDate.push_back(Date(27, Jun, 2016)); zeroCouponRate.push_back(0.0419);
			zeroCouponRateDate.push_back(Date(27, Jun, 2017)); zeroCouponRate.push_back(0.0423);
			zeroCouponRateDate.push_back(Date(27, Jun, 2018)); zeroCouponRate.push_back(0.0427);
			zeroCouponRateDate.push_back(Date(27, Jun, 2019)); zeroCouponRate.push_back(0.0431);
			zeroCouponRateDate.push_back(Date(29, Jun, 2020)); zeroCouponRate.push_back(0.0435);
			zeroCouponRateDate.push_back(Date(28, Jun, 2021)); zeroCouponRate.push_back(0.0438);
			zeroCouponRateDate.push_back(Date(27, Jun, 2022)); zeroCouponRate.push_back(0.0441);
			zeroCouponRateDate.push_back(Date(27, Jun, 2023)); zeroCouponRate.push_back(0.0443);
			zeroCouponRateDate.push_back(Date(27, Jun, 2024)); zeroCouponRate.push_back(0.0445);
			zeroCouponRateDate.push_back(Date(27, Jun, 2025)); zeroCouponRate.push_back(0.0447);
			zeroCouponRateDate.push_back(Date(29, Jun, 2026)); zeroCouponRate.push_back(0.0448);
			zeroCouponRateDate.push_back(Date(28, Jun, 2027)); zeroCouponRate.push_back(0.0450);

			EURLibor6M eulibor;
			Calendar cal = TARGET();
			DayCounter dc = eulibor.dayCounter();
			boost::shared_ptr<InterpolatedZeroCurve<LogLinear>> zeroCouponCurve
				(new InterpolatedZeroCurve<LogLinear>(zeroCouponRateDate, zeroCouponRate, dc, cal));

			boost::shared_ptr<DefaultCurve> defaultCurve(
				new DefaultCurve(DefaultCurve::ConstantIntensityGiven, Date(26, Jun, 2007), 0.05));
			calculator->InitializeDefaultCurve(defaultCurve);
			Real UCVA = calculator->getAnticUCVA(0.0405, 10 * Years, zeroCouponCurve);
			std::string msg = "[UCVA value is " + std::to_string(UCVA) + "\n";
			Logger::WriteMessage(msg.c_str());
			/*boost::shared_ptr<VanillaSwap> swap = calculator->getSwap();
			Leg floatinglegs = swap->floatingLeg();
			for (boost::shared_ptr<CashFlow> singleLeg : floatinglegs) {
				Logger::WriteMessage(std::to_string(singleLeg->date().serialNumber()).c_str());
				Logger::WriteMessage(std::to_string(singleLeg->amount()).c_str());
			}*/
		}
		/*
		TEST_METHOD(IRS_UCVA_PricingLowRisk)
		{
			std::vector<Probability> survivalProb;
			survivalProb.push_back(1.0000);
			survivalProb.push_back(0.9964);
			survivalProb.push_back(0.9834);
			survivalProb.push_back(0.9638);
			survivalProb.push_back(0.9424);
			survivalProb.push_back(0.8931);
			survivalProb.push_back(0.8164);
			survivalProb.push_back(0.7463);
			survivalProb.push_back(0.6822);
			survivalProb.push_back(0.6236);

			calculator->InitializationDefaultCurve(sampleDate, survivalProb);
			std::vector<Real> expectedAnticUCVA;
			expectedAnticUCVA.push_back(0.64);
			expectedAnticUCVA.push_back(2.52);
			expectedAnticUCVA.push_back(4.92);
			expectedAnticUCVA.push_back(7.24);
			expectedAnticUCVA.push_back(9.1);
			expectedAnticUCVA.push_back(10.51);

			std::vector<Real> expectedPostpUCVA;
			expectedPostpUCVA.push_back(0.50);
			expectedPostpUCVA.push_back(2.16);
			expectedPostpUCVA.push_back(4.47);
			expectedPostpUCVA.push_back(6.78);
			expectedPostpUCVA.push_back(8.63);
			expectedPostpUCVA.push_back(10.07);

			for (size_t i = 0; i < contractMaturity.size(); i++)
			{
				Assert::AreEqual(calculator->getAnticUCVA(contractSwapRate[i], contractMaturity[i]),
					expectedAnticUCVA[i], 0.01);
				Assert::AreEqual(calculator->getPostpUCVA(contractSwapRate[i], contractMaturity[i]),
					expectedPostpUCVA[i], 0.01);
			}
		}

		TEST_METHOD(IRS_UCVA_PricingMediaRisk)
		{
			std::vector<Probability> survivalProb;
			survivalProb.push_back(1.0000);
			survivalProb.push_back(0.9796);
			survivalProb.push_back(0.9348);
			survivalProb.push_back(0.8857);
			survivalProb.push_back(0.8371);
			survivalProb.push_back(0.7527);
			survivalProb.push_back(0.6305);
			survivalProb.push_back(0.5280);
			survivalProb.push_back(0.4423);
			survivalProb.push_back(0.3705);

			calculator->InitializationDefaultCurve(sampleDate, survivalProb);
			std::vector<Real> expectedAnticUCVA;
			expectedAnticUCVA.push_back(1.91);
			expectedAnticUCVA.push_back(6.09);
			expectedAnticUCVA.push_back(10.52);
			expectedAnticUCVA.push_back(14.51);
			expectedAnticUCVA.push_back(17.53);
			expectedAnticUCVA.push_back(19.66);

			std::vector<Real> expectedPostpUCVA;
			expectedPostpUCVA.push_back(1.80);
			expectedPostpUCVA.push_back(5.8);
			expectedPostpUCVA.push_back(10.2);
			expectedPostpUCVA.push_back(14.22);
			expectedPostpUCVA.push_back(17.28);
			expectedPostpUCVA.push_back(19.45);

			for (size_t i = 0; i < contractMaturity.size(); i++)
			{
				Assert::AreEqual(calculator->getAnticUCVA(contractSwapRate[i], contractMaturity[i]),
					expectedAnticUCVA[i], 0.01);
				Assert::AreEqual(calculator->getPostpUCVA(contractSwapRate[i], contractMaturity[i]),
					expectedPostpUCVA[i], 0.01);
			}
		}

		TEST_METHOD(IRS_UCVA_PricingHighRisk)
		{
			std::vector<Probability> survivalProb;
			survivalProb.push_back(1.0000);
			survivalProb.push_back(0.9470);
			survivalProb.push_back(0.8447);
			survivalProb.push_back(0.7478);
			survivalProb.push_back(0.6603);
			survivalProb.push_back(0.5342);
			survivalProb.push_back(0.3753);
			survivalProb.push_back(0.2636);
			survivalProb.push_back(0.1851);
			survivalProb.push_back(0.1301);

			calculator->InitializationDefaultCurve(sampleDate, survivalProb);
			std::vector<Real> expectedAnticUCVA;
			expectedAnticUCVA.push_back(4.27);
			expectedAnticUCVA.push_back(12.28);
			expectedAnticUCVA.push_back(19.55);
			expectedAnticUCVA.push_back(25.44);
			expectedAnticUCVA.push_back(29.46);
			expectedAnticUCVA.push_back(31.97);

			std::vector<Real> expectedPostpUCVA;
			expectedPostpUCVA.push_back(4.25);
			expectedPostpUCVA.push_back(12.26);
			expectedPostpUCVA.push_back(19.68);
			expectedPostpUCVA.push_back(25.77);
			expectedPostpUCVA.push_back(29.93);
			expectedPostpUCVA.push_back(32.54);

			for (size_t i = 0; i < contractMaturity.size(); i++)
			{
				Assert::AreEqual(calculator->getAnticUCVA(contractSwapRate[i], contractMaturity[i]),
					expectedAnticUCVA[i], 0.01);
				Assert::AreEqual(calculator->getPostpUCVA(contractSwapRate[i], contractMaturity[i]),
					expectedPostpUCVA[i], 0.01);
			}
		}*/
	};
}