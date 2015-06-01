#include <pch.h>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <ql/handle.hpp>
#include <ql/time/date.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/math/optimization/simplex.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace QuantLib;

namespace QuantLib
{
	TEST_CLASS(FundingSpreadTest)
	{
	public:
		TEST_METHOD(FundingSpread){
			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			Settings::instance().evaluationDate() = referenceDate;
			DayCounter daycounter = Actual360();

			double interestrate = 0.04;
			boost::shared_ptr<YieldTermStructure> discountCurve(
				new FlatForward(referenceDate, interestrate, daycounter));

			boost::shared_ptr<vector <double> > cdsspreads_(new vector<double>({ 0.0019, 0.0032, 0.0042, 0.0045, 0.0056 }));
			boost::shared_ptr<vector < Period >> tenor(new
				vector<Period>({ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<AT1Pmodel> model = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor, discountCurve));

			boost::shared_ptr<Counterparty> C(
				new Counterparty(
				0, 0.4, cdsspreads_, tenor, discountCurve, model,
				Quarterly, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
			C->modelCalibrate(Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));
		}
	};
}