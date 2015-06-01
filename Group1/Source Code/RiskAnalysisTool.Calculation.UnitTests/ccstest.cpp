#include <pch.h>
#include <limits.h>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/models/default/cirdefaultmodel.hpp>
#include <ql/models/default/deterministicdefaultmodel.hpp>
#include <ql/pricingengines/credit/mccrosscurrencyswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mcexposuremodel.hpp>
#include <ql/pricingengines/credit/mctvaengine.hpp>
#include <ql/pricingengines/credit/mcbvaengine.hpp>
#include <ql/pricingengines/credit/mccreditvarengine.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <ql/models/shortrate/calibrationhelpers/zerocouponbondhelper.hpp>
#include <ql/pricingengines/cirbondengine.hpp>
#include <ql/instruments/credit/tva.hpp>
#include <ql/instruments/credit/bva.hpp>
#include <ql/instruments/credit/var.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/currencies/america.hpp>
#include <ql/currencies/asia.hpp>
#include <ql/math/optimization/simplex.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(CCSTest)
	{
	public:

		TEST_METHOD(CCSTVAMCEngineTest)
		{
			Date d0(18, Mar, 2015);
			Settings::instance().evaluationDate() = d0;
			Real fxRate = 0.161171;
			Volatility fxVol = 0.12;
			Volatility domesticVol = 0.001;
			Volatility foreignVol = 0.001;
			CrossCurrencySwap::Type payType = CrossCurrencySwap::Type::payForeign;

			// receive information
			Real receiveNomial = 1e7;
			boost::shared_ptr<Currency>	receiveCurr(new USDCurrency());
			Rate receiveRate = 1.567594 / 100.;
			DayCounter receiveDC = Thirty360();
			Schedule paySche(d0, d0 + 5 * Years, 6 * Months, TARGET(), ModifiedFollowing, ModifiedFollowing, DateGeneration::Forward, false);


			// pay information
			Real payNomial = 62046000;
			Rate payRate = 3.294000 / 100.;
			boost::shared_ptr<Currency>	payCurr(new CNYCurrency());
			DayCounter payDC = Thirty360();
			Schedule receiveSche(d0, d0 + 5 * Years, 3 * Months, TARGET(), ModifiedFollowing, ModifiedFollowing, DateGeneration::Forward, false);

			// discount curve 
			vector<Date> discDates{
				Date(18, Mar, 2015), Date(18, Sep, 2015), Date(18, Mar, 2016),
				Date(18, Sep, 2016), Date(18, Mar, 2017), Date(18, Sep, 2017),
				Date(18, Mar, 2018), Date(18, Sep, 2018), Date(18, Mar, 2019),
				Date(18, Sep, 2019), Date(18, Mar, 2020), Date(18, Mar, 2021),
				Date(18, Mar, 2022), Date(18, Mar, 2023), Date(18, Mar, 2024),
				Date(18, Mar, 2025), Date(18, Mar, 2026), Date(18, Mar, 2027),
				Date(18, Mar, 2030), Date(18, Mar, 2035), Date(18, Mar, 2040)
			};
			vector<DiscountFactor> dfs{
				1., 0.998298, 0.99479, 0.989097, 0.981462, 0.972366,
				0.962267, 0.951817, 0.941055, 0.930079, 0.918616,
				0.895473, 0.872712, 0.850075, 0.827519, 0.805479,
				0.783812, 0.762483, 0.70229, 0.613123, 0.5379
			};

			boost::shared_ptr<YieldTermStructure> disCurve(
				new DiscountCurve(
				discDates, dfs, payDC, TARGET()));

			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss());
			//instruments zero-coupon bonds
			std::vector<boost::shared_ptr<CalibrationHelper>> zcbonds;
			boost::shared_ptr<PricingEngine> engine(new CIRBondEngine(CIRModel, Handle<YieldTermStructure>(disCurve)));
			for (int i = 1; i != discDates.size(); ++i) {
				if (discDates[i] <= d0 + 5 * Years) {
					zcbonds.push_back(boost::shared_ptr<ZerocouponbondHelper>(new ZerocouponbondHelper(
						0,
						TARGET(),
						1,
						discDates[i],
						Following,
						100.0,
						d0,
						Handle<Quote>(new SimpleQuote(0.1)),
						Handle<YieldTermStructure>(disCurve))));
					zcbonds.back()->setPricingEngine(engine);
				}
			}
			//cir model
			Simplex solver(0.001);
			const QuantLib::Size maxIteration = 10000;
			const QuantLib::Size minStatIteration = 50;
			const QuantLib::Real rootEpsilon = 1e-8;
			const QuantLib::Real FunctionEpsilon = 1e-8;
			const QuantLib::Real gradientNormEpsilon = 1e-8;
			const QuantLib::EndCriteria endcriteria = QuantLib::EndCriteria(maxIteration, minStatIteration, rootEpsilon, FunctionEpsilon, gradientNormEpsilon);
			//calibration
			CIRModel->calibrate(zcbonds, solver, endcriteria, *(CIRModel->constraint()));

			boost::shared_ptr<QuantLib::StochasticProcess1D> foreignRateProcess
				= boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]);

			double bps = 1. / 10000.;
			// investor(myself)
			boost::shared_ptr<vector <double> > investorCDSSpreads(
				new vector<double>(
			{ /*23.362,*/ 30.235, 41.876, 54.027, 65.151, 78.457, 102.7, 122.061 }
			));

			for (auto& i = investorCDSSpreads->begin(); i != investorCDSSpreads->end(); ++i)
				*i *= bps;

			boost::shared_ptr<vector < Period >> investorCDSTenor(
				new vector<Period>(
			{ /*Period(6, Months),*/ Period(1, Years), Period(2, Years), Period(3, Years),
			Period(4, Years), Period(5, Years), Period(7, Years), Period(10, Years) }
			));

			boost::shared_ptr<DefaultModel>	investorModel
				= boost::make_shared<AT1Pmodel>(investorCDSTenor, disCurve);

			/*boost::shared_ptr<DefaultModel>	investorModel
				= boost::make_shared<CirDefaultModel>();*/
			/*boost::shared_ptr<DefaultModel>	investorModel
				= boost::make_shared<DeterministicDefaultModel>();*/


			boost::shared_ptr<Counterparty> investor(
				new Counterparty(0, 0.4, investorCDSSpreads, investorCDSTenor, disCurve, investorModel,
				Quarterly, Following, DateGeneration::TwentiethIMM, d0, receiveDC, TARGET(), disCurve));

			// issurer(counterparty)
			boost::shared_ptr<vector <double> > issuerCDSSpreads(
				new vector<double>(
			{/* 19.128,*/ 26.766, 35.064, 42.953, 53.080, 68.203, 91.453, 111.109 }
			));

			for (auto& i = issuerCDSSpreads->begin(); i != issuerCDSSpreads->end(); ++i)
				*i *= bps;

			boost::shared_ptr<vector < Period >> issuerCDSTenor(
				new vector<Period>(
			{ /*Period(6, Months),*/ Period(1, Years), Period(2, Years), Period(3, Years),
			Period(4, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<DefaultModel>	issuerModel
				= boost::make_shared<AT1Pmodel>(issuerCDSTenor, disCurve);
			/*boost::shared_ptr<DefaultModel>	issuerModel
				= boost::make_shared<CirDefaultModel>();*/
			/*boost::shared_ptr<DefaultModel>	issuerModel
				= boost::make_shared<DeterministicDefaultModel>();*/

			boost::shared_ptr<Counterparty> issuer(
				new Counterparty(0, 0.4, issuerCDSSpreads, issuerCDSTenor, disCurve, issuerModel,
				Quarterly, Following, DateGeneration::TwentiethIMM, d0, receiveDC, TARGET(), disCurve));


			boost::shared_ptr<CrossCurrencySwap>
				ccs(new CrossCurrencySwap(payType, fxRate, fxVol, d0, payNomial, receiveNomial, paySche,
				payRate, payDC, receiveSche, receiveRate,
				receiveDC, payCurr, receiveCurr));

			Simplex optMethod(0.05);
			EndCriteria endCriteria(1000, 10, 1e-3, 1e-3, 1e-3);

			issuer->modelCalibrate(optMethod, endCriteria);
			investor->modelCalibrate(optMethod, endCriteria);

			Matrix correlation(4, 4, 0.);
			for (int i = 0; i < 4; ++i)
				correlation[i][i] = 1.;

			boost::shared_ptr<MCCrossCurrencySwapExposureModel> ccsModel = boost::make_shared<MCCrossCurrencySwapExposureModel>(
				TARGET(), payDC, d0, foreignRateProcess, ccs);

			std::vector < boost::shared_ptr<QuantLib::StochasticProcess1D> >
				process;
			process.push_back(issuer->createDefaultProcess());
			process.push_back(investor->createDefaultProcess());
			process.push_back(boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]));
			process.push_back(foreignRateProcess);

			boost::shared_ptr<MCTVAEngine<>> mcTVAEngine = boost::make_shared<MCTVAEngine<>>(
				std::vector<boost::shared_ptr<const Counterparty>>({ issuer, investor }), std::vector<boost::shared_ptr<const MCExposureModel>>({ ccsModel }),
				correlation, process, Handle<YieldTermStructure>(disCurve), CIRModel, 1800, 360,
				true, 500, Null<Real>(), 50000, 5);

			TVA tva;
			tva.setPricingEngine(mcTVAEngine);

			//double npv = ccs->NPV();
			double fva = tva.FVA();
			//double expectedBVA = 7169.72;
			/*double expectedCVA = 0.;
			double expectedDVA = 0.;*/

			Assert::AreEqual(0, fva, 0.);
			/*Assert::AreEqual(expectedcVA, cva, 0.0);
			Assert::AreEqual(expectedDVA, dva, 0.0);*/
		}

		TEST_METHOD(CCSBVAMCEngineTest)
		{
			Date d0(18, Mar, 2015);
			Settings::instance().evaluationDate() = d0;
			Real fxRate = 0.161171;
			Volatility fxVol = 0.4;

			CrossCurrencySwap::Type payType = CrossCurrencySwap::Type::payForeign;

			// receive information
			Real receiveNomial = 1e7;
			boost::shared_ptr<Currency>	receiveCurr(new USDCurrency());
			Rate receiveRate = 1.567594 / 100.;
			DayCounter receiveDC = Thirty360();
			Schedule paySche(d0, d0 + 5 * Years, 6 * Months, TARGET(), ModifiedFollowing, ModifiedFollowing, DateGeneration::Forward, false);


			// pay information
			Real payNomial = 62046000;
			Rate payRate = 3.294000 / 100.;
			boost::shared_ptr<Currency>	payCurr(new CNYCurrency());
			DayCounter payDC = Thirty360();
			Schedule receiveSche(d0, d0 + 5 * Years, 3 * Months, TARGET(), ModifiedFollowing, ModifiedFollowing, DateGeneration::Forward, false);

			// discount curve 
			vector<Date> discDates{
				Date(18, Mar, 2015), Date(18, Sep, 2015), Date(18, Mar, 2016),
				Date(18, Sep, 2016), Date(18, Mar, 2017), Date(18, Sep, 2017),
				Date(18, Mar, 2018), Date(18, Sep, 2018), Date(18, Mar, 2019),
				Date(18, Sep, 2019), Date(18, Mar, 2020), Date(18, Mar, 2021),
				Date(18, Mar, 2022), Date(18, Mar, 2023), Date(18, Mar, 2024),
				Date(18, Mar, 2025), Date(18, Mar, 2026), Date(18, Mar, 2027),
				Date(18, Mar, 2030), Date(18, Mar, 2035), Date(18, Mar, 2040)
			};
			vector<DiscountFactor> dfs{
				1., 0.998298, 0.99479, 0.989097, 0.981462, 0.972366,
				0.962267, 0.951817, 0.941055, 0.930079, 0.918616,
				0.895473, 0.872712, 0.850075, 0.827519, 0.805479,
				0.783812, 0.762483, 0.70229, 0.613123, 0.5379
			};

			boost::shared_ptr<YieldTermStructure> disCurve(
				new DiscountCurve(
				discDates, dfs, payDC, TARGET()));

			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss());
			//instruments zero-coupon bonds
			std::vector<boost::shared_ptr<CalibrationHelper>> zcbonds;
			boost::shared_ptr<PricingEngine> engine(new CIRBondEngine(CIRModel, Handle<YieldTermStructure>(disCurve)));
			for (int i = 1; i != discDates.size(); ++i) {
				if (discDates[i] <= d0 + 5 * Years) {
					zcbonds.push_back(boost::shared_ptr<ZerocouponbondHelper>(new ZerocouponbondHelper(
						0,
						TARGET(),
						1,
						discDates[i],
						Following,
						100.0,
						d0,
						Handle<Quote>(new SimpleQuote(0.1)),
						Handle<YieldTermStructure>(disCurve))));
					zcbonds.back()->setPricingEngine(engine);
				}
			}
			//cir model
			Simplex solver(0.001);
			const QuantLib::Size maxIteration = 10000;
			const QuantLib::Size minStatIteration = 50;
			const QuantLib::Real rootEpsilon = 1e-8;
			const QuantLib::Real FunctionEpsilon = 1e-8;
			const QuantLib::Real gradientNormEpsilon = 1e-8;
			const QuantLib::EndCriteria endcriteria = QuantLib::EndCriteria(maxIteration, minStatIteration, rootEpsilon, FunctionEpsilon, gradientNormEpsilon);
			//calibration
			CIRModel->calibrate(zcbonds, solver, endcriteria, *(CIRModel->constraint()));

			boost::shared_ptr<QuantLib::StochasticProcess1D> foreignRateProcess
				= boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]);

			double bps = 1. / 10000.;

			// investor(myself)
			boost::shared_ptr<vector <double> > investorCDSSpreads(
				new vector<double>(
			{ /*23.362,*/ 30.235, 41.876, 54.027, 65.151, 78.457, 102.7, 122.061 }
			));

			for (auto& i = investorCDSSpreads->begin(); i != investorCDSSpreads->end(); ++i)
				*i *= bps;

			boost::shared_ptr<vector < Period >> investorCDSTenor(
				new vector<Period>(
			{ /*Period(6, Months),*/ Period(1, Years), Period(2, Years), Period(3, Years),
			Period(4, Years), Period(5, Years), Period(7, Years), Period(10, Years) }
			));

			boost::shared_ptr<DefaultModel>	investorModel
				= boost::make_shared<AT1Pmodel>(investorCDSTenor, disCurve);

			boost::shared_ptr<Counterparty> investor(
				new Counterparty(0, 0.4, investorCDSSpreads, investorCDSTenor, disCurve, investorModel,
				Quarterly, Following, DateGeneration::TwentiethIMM, d0, receiveDC, TARGET(), disCurve));

			// issurer(counterparty)
			boost::shared_ptr<vector <double> > issuerCDSSpreads(
				new vector<double>(
			{/* 19.128,*/ 26.766, 35.064, 42.953, 53.080, 68.203, 91.453, 111.109 }
			));

			for (auto& i = issuerCDSSpreads->begin(); i != issuerCDSSpreads->end(); ++i)
				*i *= bps;

			boost::shared_ptr<vector < Period >> issuerCDSTenor(
				new vector<Period>(
			{ /*Period(6, Months),*/ Period(1, Years), Period(2, Years), Period(3, Years),
			Period(4, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<DefaultModel>	issuerModel
				= boost::make_shared<AT1Pmodel>(issuerCDSTenor, disCurve);

			boost::shared_ptr<Counterparty> issuer(
				new Counterparty(0, 0.4, issuerCDSSpreads, issuerCDSTenor, disCurve, issuerModel,
				Quarterly, Following, DateGeneration::TwentiethIMM, d0, receiveDC, TARGET(), disCurve));


			boost::shared_ptr<CrossCurrencySwap>
				ccs(new CrossCurrencySwap(payType, fxRate, fxVol, d0, payNomial, receiveNomial, paySche,
				payRate, payDC, receiveSche, receiveRate,
				receiveDC, payCurr, receiveCurr));

			Simplex optMethod(0.05);
			EndCriteria endCriteria(1000, 10, 1e-3, 1e-3, 1e-3);

			issuer->modelCalibrate(optMethod, endCriteria);
			investor->modelCalibrate(optMethod, endCriteria);

			Matrix correlation(4, 4, 0.);
			for (int i = 0; i < 4; ++i)
				correlation[i][i] = 1.;

			boost::shared_ptr<MCCrossCurrencySwapExposureModel> ccsModel = boost::make_shared<MCCrossCurrencySwapExposureModel>(
				TARGET(), payDC, d0, foreignRateProcess, ccs);

			std::vector < boost::shared_ptr<QuantLib::StochasticProcess1D> >
				process;
			process.push_back(issuer->createDefaultProcess());
			process.push_back(investor->createDefaultProcess());
			process.push_back(boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]));
			process.push_back(foreignRateProcess);

			boost::shared_ptr<MCBVAEngine<>> mcBVAEngine = boost::make_shared<MCBVAEngine<>>(
				std::vector<boost::shared_ptr<const Counterparty>>({ issuer, investor }), std::vector<boost::shared_ptr<const MCExposureModel>>({ ccsModel }),
				correlation, process, Handle<YieldTermStructure>(disCurve), CIRModel, 1800, 360,
				true, 10000, Null<Real>(), 50000, 5);

			BVA bva;
			bva.setPricingEngine(mcBVAEngine);

			double cva = bva.CVA();

			double dva = bva.DVA();

			Assert::AreEqual(-2200, cva - dva, -2000 * 0.2);

		}

	};
}