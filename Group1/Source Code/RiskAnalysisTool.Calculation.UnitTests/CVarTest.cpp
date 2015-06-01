#include <pch.h>
//#include <ql/credit/creditvar/creditvarengine.hpp>
//
//using namespace std;
//using namespace QuantLib;
//using namespace Microsoft::VisualStudio::CppUnitTestFramework;
//
//namespace QuantLib
//{
//	TEST_CLASS(CreditVaRComputationEngineTests)
//	{
//	public:
//		TEST_METHOD(ConstructorTest)
//		{
//			shared_ptr<CreditVaRComputationEngine> CVar1(new CreditVaRComputationEngine());
//			Assert::IsNotNull(CVar1.get());
//		}
//
//		TEST_METHOD(CalculateTest)
//		{
//			Handle<Quote> underlying_(
//				boost::shared_ptr<Quote>(
//				new SimpleQuote(23)
//				));
//
//			Handle<Quote> default_(
//				boost::shared_ptr<Quote>(
//				new SimpleQuote(1)
//				));
//
//			Calendar calendar = TARGET();
//			DayCounter dayCounter = Actual365Fixed();
//			Date t0(29, Dec, 2011);
//			Date VarDate(29, Dec, 2012);
//			Date MaturityDate(29, Dec, 2013);
//
//			std::vector<QuantLib::Time> cdsMaturity;
//			cdsMaturity.push_back(dayCounter.yearFraction(t0, MaturityDate));
//
//			Handle<YieldTermStructure> flatTermStructureQ(
//				boost::shared_ptr<YieldTermStructure>(
//				new FlatForward(t0, 0.02, dayCounter)
//				));
//			Handle<YieldTermStructure> flatTermStructureP(
//				boost::shared_ptr<YieldTermStructure>(
//				new FlatForward(t0, 0.03, dayCounter)
//				));
//			Handle<BlackVolTermStructure> flatVolTS(
//				boost::shared_ptr<BlackVolTermStructure>(
//				new BlackConstantVol(t0, calendar, 0.2, dayCounter)
//				));
//
//
//			shared_ptr<QuantLib::BlackScholesProcess> defaultProcess_(
//				new BlackScholesProcess(default_, flatTermStructureP, flatVolTS));
//			shared_ptr<QuantLib::BlackScholesProcess> defaultQProcess_(
//				new BlackScholesProcess(default_, flatTermStructureQ, flatVolTS));
//
//			shared_ptr<QuantLib::BlackScholesProcess> underlyingProcess_(
//				new BlackScholesProcess(underlying_, flatTermStructureP, flatVolTS));
//			shared_ptr<QuantLib::BlackScholesProcess> underlyingQProcess_(
//				new BlackScholesProcess(underlying_, flatTermStructureQ, flatVolTS));
//
//			//set payoff and need to set the ExerciseType as global variable
//			Option::Type type(Option::Call);
//			Real Strike = 25;
//			shared_ptr<StrikedTypePayoff> payoff_(
//				new PlainVanillaPayoff(type, Strike)
//				);
//			
//			shared_ptr<CreditVaRComputationEngine> CVar(new CreditVaRComputationEngine());
//
//			CVar->setCorrelation(0.2);
//			CVar->setLGD(0.6);
//			CVar->setConfidence(0.99);
//			CVar->setDefaultBarrier(0.4);
//
//			CVar->setReferenceDate(t0);
//			CVar->setVarDate(VarDate);
//			CVar->setMaturityDate(MaturityDate);
//			CVar->setDayCounter(dayCounter);
//
//			CVar->setDefaultProcess(defaultProcess_);
//			CVar->setDefaultQProcess(defaultQProcess_);
//			CVar->setUnderlyingProcess(underlyingProcess_);
//			CVar->setUnderlyingQProcess(underlyingQProcess_);
//			CVar->setPayoff(payoff_);
//			CVar->setNumberOfPath(1000);
//			CVar->setSteps(2000);
//
//			CVar->calculate(cdsMaturity,"AT1P",Strike,type);
//
//			vector<double> dist;
//			dist = *(CVar->lossDistribution_);
//			double CVarValue = CVar->getVaR();
//			double ESValue = CVar->getES();
//			double ExpectedVaR = 0.0;
//
//			Assert::AreEqual(ExpectedVaR, CVarValue, 0.);
//		}
//	};
//}