#include <pch.h>

#include <limits>
#include <ql/instruments/swap/equityswap.hpp>
#include <ql/pricingengines/swap/analyticequityswapengine.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/pricingengines/credit/mcequityswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mctvaengine.hpp>
#include <ql/instruments/credit/tva.hpp>
#include <ql/pricingengines/credit/mcbvaengine.hpp>
#include <ql/instruments/credit/bva.hpp>
#include <ql/pricingengines/credit/mccreditvarengine.hpp>
#include <ql/instruments/credit/var.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/pricingengines/ucva/mcequityswapucvaengine.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace QuantLib;

namespace QuantLib
{
	struct EquitySwapData {
		EquitySwap::Type type;
		Real startprice;	// start price
		Real spotprice;     // spot price
		Integer amount;     // number of stock
		Real maturity;      // time to maturity
		Real fixedrate;		// fixed rate
		Real dividend;		// equity dividend yield
		Volatility v;		// volatility of underlying equity
		Real interest;		// risk free interest
		Date referencedate;	// reference date
		Date startdate;		// start date of equity swap
		// Daycounter
	};

	EquitySwapData data[] = {
		//
		{ EquitySwap::Type::Payer, 35.12, 35.12, 100, 5, 0.04, 0.03, 0.2, 0.04, Date(5, Mar, 2015), Date(5, Mar, 2015) },
	};

	TEST_CLASS(EquitySwapTest)
	{
	public:

		TEST_METHOD(AnalyticEquitySwapEngineTest){
			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			Settings::instance().evaluationDate() = referenceDate;
			boost::shared_ptr<EquitySwap> ers1 = boost::shared_ptr<EquitySwap>(new EquitySwap(
				data[0].type, data[0].startdate, data[0].referencedate, data[0].startprice, data[0].spotprice, data[0].amount,
				0, 0, data[0].fixedrate, data[0].dividend, data[0].maturity, data[0].v,
				boost::shared_ptr<YieldTermStructure>(new FlatForward(data[0].referencedate, data[0].interest, Actual360()))));
			double notional = 3512;
			boost::shared_ptr<EquitySwap::engine> engine(new AnalyticESEngine());

			ers1->setPricingEngine(engine);
			double actual = ers1->NPV();
			Assert::AreEqual(actual, 0, notional*0.01);
		}

		TEST_METHOD(EquitySwapUCVAengineTest){
			// from Brigo's book page 171


			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			DayCounter daycounter = Actual360();
			Settings::instance().evaluationDate() = referenceDate;
			
			double interestrate = 0.04;
			double sigma = 0.2;
			boost::shared_ptr<YieldTermStructure> discountCurve(
				new FlatForward(referenceDate, interestrate, daycounter));


			// construct counterparty
			boost::shared_ptr<vector <double> > cdsspreads_(new vector<double>({ 0.0019, 0.0032, 0.0042, 0.0045, 0.0056 }));
			boost::shared_ptr<vector < Period >> tenor(new
				vector<Period>({ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<AT1Pmodel> model = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<Counterparty> C(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
			C->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			// construct equity return swap
			double StockPrice = 20.0;
			int amount = 1;
			double stockdividend = 0.008;

			double spread = 0.0;
			double fixedrate = 0.04;
			double maturity = 5.0;

			boost::shared_ptr<EquitySwap> ers = boost::make_shared<EquitySwap>(
				EquitySwap::Type::Payer,
				referenceDate,
				referenceDate,
				StockPrice,
				StockPrice,
				amount,
				0,
				0,
				fixedrate,
				stockdividend,
				maturity,
				sigma,
				discountCurve);

			boost::shared_ptr<EquitySwap> ersfordva = boost::make_shared<EquitySwap>(
				EquitySwap::Type::Receiver,
				referenceDate,
				referenceDate,
				StockPrice,
				StockPrice,
				amount,
				0,
				0,
				fixedrate,
				stockdividend,
				maturity,
				sigma,
				discountCurve);

			boost::shared_ptr<EquitySwap::engine> engine(new AnalyticESEngine());
			ers->setPricingEngine(engine);
			ersfordva->setPricingEngine(engine);
			double riskfreeprice1 = ers->NPV();
			double riskfreeprice2 = ersfordva->NPV();

			Matrix correlation(2, 2, 0.);
			correlation[0][0] = 1;
			correlation[1][1] = 1;

			boost::shared_ptr<EquitySwap::engine> engineucva(new MCEquitySwapUCVAEngine(
				Handle<YieldTermStructure>(discountCurve), 
				ers, 
				C, 
				correlation,
				360,
				true, 10000, Null<Real>(), 50000));

			ers->setPricingEngine(engineucva);
			ersfordva->setPricingEngine(engineucva);
			double payerucvaengineresult = ers->result<Real>("UCVA");
			double receiverucvaengineresult = ersfordva->result<Real>("UCVA");

			Assert::AreEqual(0.04, payerucvaengineresult, 0.01);
			Assert::AreEqual(0.02, receiverucvaengineresult, 0.01);

		}

		TEST_METHOD(EquitySwapBVAengineTest){
			// from Brigo's book page 171
			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			Settings::instance().evaluationDate() = referenceDate;
			DayCounter daycounter = Actual360();
			double interestrate = 0.04;
			double sigma = 0.2;
			boost::shared_ptr<YieldTermStructure> discountCurve(
				new FlatForward(referenceDate, interestrate, daycounter));

			// construct counterparty
			boost::shared_ptr<vector <double> > cdsspreads_(new vector<double>({ 0.0019, 0.0032, 0.0042, 0.0045, 0.0056 }));
			boost::shared_ptr<vector < Period >> tenor(new
				vector<Period>({ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<AT1Pmodel> model1 = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<AT1Pmodel> model2 = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<Counterparty> C1(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model1,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
			boost::shared_ptr<Counterparty> C2(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model2,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
			C1->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));
			C2->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			// do not need calibrate (reduced situation)
			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss(interestrate, interestrate, 0.1, QL_EPSILON));

			// construct equity return swap
			double StockPrice = 20.0;
			int amount = 1;
			double stockdividend = 0.008;

			double fixedrate = 0.04;
			double maturity = 5.0;

			boost::shared_ptr<EquitySwap> ers = boost::make_shared<EquitySwap>(
				EquitySwap::Type::Payer,
				referenceDate,
				referenceDate,
				StockPrice,
				StockPrice,
				amount,
				0,
				0,
				fixedrate,
				stockdividend,
				maturity,
				sigma,
				discountCurve);

			std::vector<boost::shared_ptr<StochasticProcess1D>> process;
			process.push_back(C1->createDefaultProcess());
			process.push_back(C2->createDefaultProcess());
			process.push_back(boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]));
			process.push_back(
				boost::shared_ptr<BlackScholesMertonProcess>(new BlackScholesMertonProcess(
				Handle<Quote>(new SimpleQuote(ers->getSpotPrice())),
				Handle<YieldTermStructure>(new FlatForward(
				referenceDate, ers->getDividendYield(), daycounter)),
				Handle<YieldTermStructure>(discountCurve),
				Handle<BlackVolTermStructure>(new BlackConstantVol(
				referenceDate, ers->getCalendar(),
				ers->getVolatility(), daycounter)))));

			boost::shared_ptr<MCEquitySwapExposureModel> ersModel = boost::make_shared<MCEquitySwapExposureModel>(
				calendar, daycounter, referenceDate, process.back(), ers);

			Matrix correlation(4, 4, 0.);
			for (int i = 0; i < correlation.rows(); ++i)
				correlation[i][i] = 1.;

			boost::shared_ptr<MCBVAEngine<>> mcBVAEngine = boost::make_shared<MCBVAEngine<>>(
				std::vector<boost::shared_ptr<const Counterparty>>({ C1, C2 }),
				std::vector<boost::shared_ptr<const MCExposureModel>>({ ersModel }),
				correlation, process, Handle<YieldTermStructure>(discountCurve), CIRModel, 1800, 360,
				true, 20000, Null<Real>(), 50000, 5);

			BVA bva;

			bva.setPricingEngine(mcBVAEngine);

			double CVA = bva.CVA();
			double DVA = bva.DVA();

			Assert::AreEqual(0.03, CVA, 0.01);
			Assert::AreEqual(0.02, DVA, 0.01);
		}

		TEST_METHOD(EquitySwapTVAengineTest){
			// from Brigo's book page 171
			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			Settings::instance().evaluationDate() = referenceDate;
			DayCounter daycounter = Actual360();
			double interestrate = 0.04;
			double sigma = 0.2;
			boost::shared_ptr<YieldTermStructure> discountCurve(
				new FlatForward(referenceDate, interestrate, daycounter));

			// construct counterparty
			boost::shared_ptr<vector <double> > cdsspreads_(new vector<double>({ 0.0019, 0.0032, 0.0042, 0.0045, 0.0056 }));
			boost::shared_ptr<vector < Period >> tenor(new
				vector<Period>({ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<AT1Pmodel> model1 = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<AT1Pmodel> model2 = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<Counterparty> C1(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model1,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar, discountCurve));
			boost::shared_ptr<Counterparty> C2(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model2,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar, discountCurve));
			C1->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));
			C2->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			// do not need calibrate (reduced situation)
			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss(interestrate, interestrate, 0.1, QL_EPSILON));

			// construct equity return swap
			double StockPrice = 20.0;
			int amount = 1;
			double stockdividend = 0.008;

			double spread = 0.0;
			double fixedrate = 0.04;
			double maturity = 5.0;

			boost::shared_ptr<EquitySwap> ers = boost::make_shared<EquitySwap>(
				EquitySwap::Type::Payer,
				referenceDate,
				referenceDate,
				StockPrice,
				StockPrice,
				amount,
				0,
				0,
				fixedrate,
				stockdividend,
				maturity,
				sigma,
				discountCurve);

			std::vector<boost::shared_ptr<StochasticProcess1D>> process;
			process.push_back(C1->createDefaultProcess());
			process.push_back(C2->createDefaultProcess());
			process.push_back(boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]));
			process.push_back(
				boost::shared_ptr<BlackScholesMertonProcess>(new BlackScholesMertonProcess(
				Handle<Quote>(new SimpleQuote(ers->getSpotPrice())),
				Handle<YieldTermStructure>(new FlatForward(
				referenceDate, ers->getDividendYield(), daycounter)),
				Handle<YieldTermStructure>(discountCurve),
				Handle<BlackVolTermStructure>(new BlackConstantVol(
				referenceDate, ers->getCalendar(),
				ers->getVolatility(), daycounter)))));

			boost::shared_ptr<MCEquitySwapExposureModel> ersModel = boost::make_shared<MCEquitySwapExposureModel>(
				calendar, daycounter, referenceDate, process.back(), ers);

			Matrix correlation(4, 4, 0.);
			for (int i = 0; i < correlation.rows(); ++i)
				correlation[i][i] = 1.;

			boost::shared_ptr<MCTVAEngine<>> mcTVAEngine = boost::make_shared<MCTVAEngine<>>(
				std::vector<boost::shared_ptr<const Counterparty>>({ C1, C2 }),
				std::vector<boost::shared_ptr<const MCExposureModel>>({ ersModel }),
				correlation, process, Handle<YieldTermStructure>(discountCurve), CIRModel, 1800, 360,
				false, 500, Null<Real>(), 50000, 5);

			TVA tva;
			tva.setPricingEngine(mcTVAEngine);

			double FVA = tva.FVA();

			Assert::AreEqual(-0.01, FVA, 0.1 * FVA);

		}


		TEST_METHOD(EquitySwapCreditVaRTest){
			// from Brigo's book page 171
			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			Settings::instance().evaluationDate() = referenceDate;
			DayCounter daycounter = Actual360();
			double interestrate = 0.04;
			double sigma = 0.2;
			boost::shared_ptr<YieldTermStructure> discountCurve(
				new FlatForward(referenceDate, interestrate, daycounter));

			// construct counterparty
			boost::shared_ptr<vector <double> > cdsspreads_(new vector<double>({ 0.0019, 0.0032, 0.0042, 0.0045, 0.0056 }));
			boost::shared_ptr<vector < Period >> tenor(new
				vector<Period>({ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<AT1Pmodel> model1 = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<AT1Pmodel> model2 = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));
			boost::shared_ptr<Counterparty> C1(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model1,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
			boost::shared_ptr<Counterparty> C2(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, model2,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
			C1->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));
			C2->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			// do not need calibrate (reduced situation)
			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss(interestrate, interestrate, 0.1, QL_EPSILON));

			// construct equity return swap
			double StockPrice = 20.0;
			int amount = 1;
			double stockdividend = 0.008;
			double spread = 0.0;
			double fixedrate = 0.04;
			double maturity = 5.0;

			boost::shared_ptr<EquitySwap> ers = boost::make_shared<EquitySwap>(
				EquitySwap::Type::Payer,
				referenceDate,
				referenceDate,
				StockPrice,
				StockPrice,
				amount,
				0,
				0,
				fixedrate,
				stockdividend,
				maturity,
				sigma,
				discountCurve);

			Matrix correlation(3, 3, 0.);
			for (int i = 0; i < 3; ++i)
				correlation[i][i] = 1.;

			std::vector < boost::shared_ptr<QuantLib::StochasticProcess1D> >
				process;
			process.push_back(C1->createDefaultProcess());
			process.push_back(boost::make_shared<QuantLib::CIRprocess>(CIRModel->params()[0], CIRModel->params()[1], CIRModel->params()[2], CIRModel->params()[3]));
			process.push_back(
				boost::shared_ptr<BlackScholesMertonProcess>(new BlackScholesMertonProcess(
				Handle<Quote>(new SimpleQuote(ers->getSpotPrice())),
				Handle<YieldTermStructure>(new FlatForward(
				referenceDate, ers->getDividendYield(), daycounter)),
				Handle<YieldTermStructure>(discountCurve),
				Handle<BlackVolTermStructure>(new BlackConstantVol(
				referenceDate, ers->getCalendar(),
				ers->getVolatility(), daycounter)))));

			boost::shared_ptr<MCEquitySwapExposureModel> ersModel = boost::make_shared<MCEquitySwapExposureModel>(
				calendar, daycounter, referenceDate, process.back(), ers);

			boost::shared_ptr<MCCreditVaREngine<>> mcCVarEngine = boost::make_shared<MCCreditVaREngine<>>(
				std::vector<boost::shared_ptr<const Counterparty>>({ C1, C2 }),
				std::vector<boost::shared_ptr<const MCExposureModel>>({ ersModel }),
				correlation, process, Handle<YieldTermStructure>(discountCurve), CIRModel, Null<Size>(), 360,
				false, 100000, Null<Real>(), 100000, 1);

			VAR var;
			var.setPricingEngine(mcCVarEngine);

			double creditvar = var.CreditVAR();

			Assert::AreEqual(0, creditvar, 0.);

		}

	};
}