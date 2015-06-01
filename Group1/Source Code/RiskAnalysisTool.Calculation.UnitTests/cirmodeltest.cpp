#include <pch.h>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <ql/models/shortrate/calibrationhelpers/zerocouponbondhelper.hpp>
#include <ql/pricingengines/cirbondengine.hpp>
#include <ql/termstructures/yield/discountcurve.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/processes/cirprocess.hpp>
#include <ql/math/randomnumbers/all.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib
{
	TEST_CLASS(cirmodeltest)
	{
	public:
		TEST_METHOD(CIRModelCalibration)
		{
			Calendar calendar = TARGET();
			DayCounter dayCounter = Actual360();
			Date referenceDate(18, Mar, 2015);
			Settings::instance().evaluationDate() = referenceDate;

			// discount curve
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
			boost::shared_ptr<YieldTermStructure> disTS(new DiscountCurve(dates, discount, dayCounter));
			
			//cir model
			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss());
			//instruments zero-coupon bonds
			std::vector<boost::shared_ptr<CalibrationHelper>> zcbonds;
			vector<Real> weights;
			boost::shared_ptr<PricingEngine> engine(new CIRBondEngine(CIRModel, Handle<YieldTermStructure>(disTS)));
			for (int i = 0; i != 10; ++i){
				zcbonds.push_back(boost::shared_ptr<ZerocouponbondHelper>(new ZerocouponbondHelper(
					0,
					calendar,
					1,
					dates[2*i + 1],
					Following,
					100.0,
					referenceDate,
					Handle<Quote>(new SimpleQuote(0.1)),
					Handle<YieldTermStructure>(disTS))));
				zcbonds[i]->setPricingEngine(engine);
				weights.push_back(0.1);
			}
			//cir model
			Simplex solver(0.001);
			const QuantLib::Size maxIteration = 1000;
			const QuantLib::Size minStatIteration = 50;
			const QuantLib::Real rootEpsilon = 1e-8;
			const QuantLib::Real FunctionEpsilon = 1e-8;
			const QuantLib::Real gradientNormEpsilon = 1e-8;
			const QuantLib::EndCriteria endcriteria = QuantLib::EndCriteria(maxIteration, minStatIteration, rootEpsilon, FunctionEpsilon, gradientNormEpsilon);

			CIRModel->calibrate(zcbonds, solver, endcriteria, *(CIRModel->constraint()), weights);

			Real theta = CIRModel->params()[0];
			Real k = CIRModel->params()[1];
			Real sigma = CIRModel->params()[2];
			Real r0 = CIRModel->params()[3];

			vector<Real> actualdiscount;
			for (int i = 0; i != 6; ++i){
				actualdiscount.push_back(CIRModel->discountBond(0., dayCounter.yearFraction(referenceDate, dates[2*i]), r0));
			}

			Assert::AreEqual(discount[0], actualdiscount[0], discount[0] * 0.01);
			Assert::AreEqual(discount[2], actualdiscount[1], discount[2] * 0.01);
			Assert::AreEqual(discount[4], actualdiscount[2], discount[4] * 0.01);
			Assert::AreEqual(discount[6], actualdiscount[3], discount[6] * 0.01);
			Assert::AreEqual(discount[8], actualdiscount[4], discount[8] * 0.01);
		}

		TEST_METHOD(CIRProcessEvolveErrorTest)
		{
			// Implicit Milstein Discretization Scheme
			Calendar calendar = TARGET();
			DayCounter dayCounter = Actual360();
			Date referenceDate(18, Mar, 2015);
			Settings::instance().evaluationDate() = referenceDate;

			// discount curve
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
			boost::shared_ptr<YieldTermStructure> disTS(new DiscountCurve(dates, discount, dayCounter));

			//cir model
			boost::shared_ptr<CoxIngersollRoss> CIRModel(new CoxIngersollRoss());
			//instruments zero-coupon bonds
			std::vector<boost::shared_ptr<CalibrationHelper>> zcbonds;
			vector<Real> weights;
			boost::shared_ptr<PricingEngine> engine(new CIRBondEngine(CIRModel, Handle<YieldTermStructure>(disTS)));
			for (int i = 0; i != 10; ++i){
				zcbonds.push_back(boost::shared_ptr<ZerocouponbondHelper>(new ZerocouponbondHelper(
					0,
					calendar,
					1,
					dates[2 * i + 1],
					Following,
					100.0,
					referenceDate,
					Handle<Quote>(new SimpleQuote(0.1)),
					Handle<YieldTermStructure>(disTS))));
				zcbonds[i]->setPricingEngine(engine);
				weights.push_back(0.1);
			}
			//cir model
			Simplex solver(0.001);
			const QuantLib::Size maxIteration = 1000;
			const QuantLib::Size minStatIteration = 50;
			const QuantLib::Real rootEpsilon = 1e-8;
			const QuantLib::Real FunctionEpsilon = 1e-8;
			const QuantLib::Real gradientNormEpsilon = 1e-8;
			const QuantLib::EndCriteria endcriteria = QuantLib::EndCriteria(maxIteration, minStatIteration, rootEpsilon, FunctionEpsilon, gradientNormEpsilon);

			CIRModel->calibrate(zcbonds, solver, endcriteria, *(CIRModel->constraint()), weights);

			Real theta = CIRModel->params()[0];
			Real k = CIRModel->params()[1];
			Real sigma = CIRModel->params()[2];
			Real r0 = CIRModel->params()[3];

			boost::shared_ptr<CIRprocess>cirprocess(new CIRprocess(theta, k, sigma, r0, CIRprocess::ImplicitMilstein));

			Time dt = dayCounter.yearFraction(dates[0], dates[1]) / 60.0;

			QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
			QuantLib::MersenneTwisterUniformRng URng(seed);
			QuantLib::BoxMullerGaussianRng < QuantLib::MersenneTwisterUniformRng > dffusionRng(URng);

			Real Rate = 0.0;
			Real means = 0.0;
			Real msq = 0.0;
			for (int n = 0; n != 100000; ++n){
				Real tmp = r0;
				for (int i = 0; i != 60; ++i){
					Real x = dffusionRng.next().value;
					tmp = cirprocess->evolve(i*dt, tmp, dt, x);
				}
				Rate += tmp;
			}
			Rate = Rate / 100000.0;

			double expected = r0 * std::exp(-k * dt * 60) + theta * (1 - std::exp(-k * dt * 60));

			Assert::AreEqual(expected, Rate, expected * 0.05);

		}
	};
}