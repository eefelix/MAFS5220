#include <pch.h>
#include <limits.h>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/pricingengines/credit/mcvanillaswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mcexposuremodel.hpp>
#include <ql/pricingengines/credit/mccreditvarengine.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <ql/models/shortrate/calibrationhelpers/zerocouponbondhelper.hpp>
#include <ql/pricingengines/cirbondengine.hpp>
#include <ql/instruments/credit/var.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/indexes/ibor/euribor.hpp>
#include <ql/processes/cirprocess.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	void TestCVarCalculation(
		Date referenceDate, Date maturityDate,
		vector<Date> disTsDates, vector<double> disFactors,
		boost::shared_ptr<vector <double> > cdsspreads1,
		boost::shared_ptr<vector < Period >> tenor1,
		boost::shared_ptr<vector <double> > cdsspreads2,
		boost::shared_ptr<vector < Period >> tenor2,
		Schedule& fixedSchedule,
		Schedule& floatSchedule,
		double nominal,
		double fixedRate,
		double spread,
		VanillaSwap::Type swapType,
		double expectedValue){

		Calendar calendar = TARGET();
		DayCounter dayCounter = Actual360();
		Settings::instance().evaluationDate() = referenceDate;
		// discount curve
		boost::shared_ptr<YieldTermStructure> disTS(new DiscountCurve(disTsDates, disFactors, dayCounter));

		//counterparty
		boost::shared_ptr<DefaultModel>	issuerModel
			= boost::make_shared<AT1Pmodel>(tenor1, disTS);
		boost::shared_ptr<Counterparty> issuer(
			new Counterparty(0, 0.4, cdsspreads1, tenor1, disTS, issuerModel,
			Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, dayCounter, TARGET(), disTS));
		issuer->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

		//investor
		boost::shared_ptr<DefaultModel>	investorModel
			= boost::make_shared<AT1Pmodel>(tenor2, disTS);
		boost::shared_ptr<Counterparty> investor(
			new Counterparty(0, 0.4, cdsspreads2, tenor2, disTS, investorModel,
			Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, dayCounter, TARGET(), disTS));
		investor->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

		Matrix correlation(2, 2, 0.);
		for (int i = 0; i < 2; ++i)
			correlation[i][i] = 1.;

		//swap information
		boost::shared_ptr<QuantLib::IborIndex> ibor(new QuantLib::Euribor(3 * Months, Handle<YieldTermStructure>(disTS)));
		ibor->addFixing(ibor->fixingDate(referenceDate), 0., true);

		boost::shared_ptr<VanillaSwap> swap(new VanillaSwap(swapType, nominal, fixedSchedule, fixedRate, Thirty360(),
			floatSchedule, ibor, spread, Actual360()));

		boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss());
		//instruments zero-coupon bonds
		std::vector<boost::shared_ptr<CalibrationHelper>> zcbonds;
		boost::shared_ptr<PricingEngine> engine(new CIRBondEngine(CIRModel, Handle<YieldTermStructure>(disTS)));
		for (int i = 1; i != disTsDates.size(); ++i) {
			if (disTsDates[i] <= maturityDate) {
				zcbonds.push_back(boost::shared_ptr<ZerocouponbondHelper>(new ZerocouponbondHelper(
					0,
					calendar,
					1,
					disTsDates[i],
					Following,
					100.0,
					referenceDate,
					Handle<Quote>(new SimpleQuote(0.1)),
					Handle<YieldTermStructure>(disTS))));
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

		//Setting engine
		boost::shared_ptr<MCVanillaSwapExposureModel> irsModel(new MCVanillaSwapExposureModel(
			calendar, dayCounter, referenceDate, swap));

		std::vector < boost::shared_ptr<QuantLib::StochasticProcess1D> >
			process;
		process.push_back(issuer->createDefaultProcess());
		process.push_back(boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]));

		boost::shared_ptr<MCCreditVaREngine<>> mcCVarEngine = boost::make_shared<MCCreditVaREngine<>>(
			std::vector<boost::shared_ptr<const Counterparty>>({ issuer, investor }),
			std::vector<boost::shared_ptr<const MCExposureModel>>({ irsModel }),
			correlation, process, Handle<YieldTermStructure>(disTS), CIRModel, Null<Size>(), 360,
			false, 100000, Null<Real>(), 100000, 1);

		VAR var;
		var.setPricingEngine(mcCVarEngine);

		double creditvar = var.CreditVAR();

		Assert::AreEqual(expectedValue, creditvar, 0.1 * expectedValue);
	}

	TEST_CLASS(IRSCVarTestCase)
	{
		TEST_METHOD(IRSCVarTestSpreadSheet)
		{
			Calendar calendar = TARGET();
			DayCounter dayCounter = Actual360();
			Date referenceDate(18, Mar, 2015);
			Date maturityDate(18, Mar, 2020);
			vector<Date> dates({
				Date(18, Mar, 2015), Date(18, Jun, 2015), Date(18, Sep, 2015), Date(18, Dec, 2015),
				Date(18, Mar, 2016), Date(20, Jun, 2016), Date(19, Sep, 2016), Date(19, Dec, 2016),
				Date(20, Mar, 2017), Date(19, Jun, 2017), Date(18, Sep, 2017), Date(18, Dec, 2017),
				Date(19, Mar, 2018), Date(18, Jun, 2018), Date(18, Sep, 2018), Date(18, Dec, 2018),
				Date(18, Mar, 2019), Date(18, Jun, 2019), Date(18, Sep, 2019), Date(18, Dec, 2019),
				Date(18, Mar, 2020), Date(18, Mar, 2021)
			});
			vector<double> discount({
				1., 0.99931, 0.998298, 0.996803, 0.994802, 0.99218,
				0.989108, 0.985539, 0.981484, 0.977083, 0.972384,
				0.967452, 0.962292, 0.956992, 0.951799, 0.946447,
				0.940943, 0.935505, 0.929912, 0.924224, 0.918389,
				0.895473
			});
			//issuer
			boost::shared_ptr<vector <double> > cdsspreads1(new vector<double>({ 0.031621, 0.050727, 0.090606, 0.117975, 0.126978 }));
			boost::shared_ptr<vector < Period >> tenor1(new
				vector<Period>({ Period(1, Years), Period(2, Years), Period(3, Years), Period(4, Years), Period(5, Years) }));

			//investor
			boost::shared_ptr<vector <double> > cdsspreads2(new vector<double>({ 0.034021, 0.04743, 0.059706, 0.065624, 0.077995 }));
			boost::shared_ptr<vector < Period >> tenor2(new
				vector<Period>({ Period(1, Years), Period(2, Years), Period(3, Years), Period(4, Years), Period(5, Years) }));

			Schedule fixedSchedule(referenceDate, Date(18, Mar, 2020), 6 * Months, calendar, Following, Following, DateGeneration::Forward, false);
			Schedule floatSchedule(referenceDate, Date(18, Mar, 2020), 3 * Months, calendar, Following, Following, DateGeneration::Forward, false);
			double nominal = 1.0;
			double fixedRate = 0.01551176;
			double spread = 0.;
			VanillaSwap::Type swapType = VanillaSwap::Payer;

			double expectedValue = 0.0048;

			TestCVarCalculation(
				referenceDate,
				maturityDate,
				dates, discount,
				cdsspreads1, tenor1,
				cdsspreads2, tenor2,
				fixedSchedule, floatSchedule,
				nominal,
				fixedRate,
				spread,
				swapType,
				expectedValue);
		}
	};
}