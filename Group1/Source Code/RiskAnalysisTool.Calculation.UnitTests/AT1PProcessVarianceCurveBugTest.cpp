#include <pch.h>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/settings.hpp>
#include <ql/time/date.hpp>
#include <ql/time/calendar.hpp>
#include <ql/time/period.hpp>
#include <ql/processes/blackscholesprocess.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib{
	TEST_CLASS(AT1PProcessVarianceBugTest)
	{
	public:
		TEST_METHOD(AT1PBugTest)
		{
			Date referenceDate = Date(10, Mar, 2004);
			Calendar calendar = TARGET();
			Settings::instance().evaluationDate() = referenceDate;
			DayCounter daycounter = Actual365Fixed();
			double interestrate = 0.0000001;
			boost::shared_ptr<YieldTermStructure> discountCurve(
				new FlatForward(referenceDate, interestrate, daycounter));


			// construct counterparty
			boost::shared_ptr<vector <double> > cdsspreads_(new vector<double>({ 0.0019, 0.0032, 0.0042, 0.0045, 0.0056 }));
			boost::shared_ptr<vector < Period >> tenor(new
				vector<Period>({ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			boost::shared_ptr<AT1Pmodel> at1pmodel = boost::shared_ptr<AT1Pmodel>(new AT1Pmodel(tenor,discountCurve));
			boost::shared_ptr<Counterparty> C(
				new Counterparty(0, 0.4, cdsspreads_, tenor, discountCurve, at1pmodel,
				Semiannual, Following, DateGeneration::TwentiethIMM, referenceDate, daycounter, calendar));
		
			double T = daycounter.yearFraction(referenceDate, referenceDate + tenor->back());
			long nSteps = 50;
			double t = 0;
			double dt = T / double(nSteps);
			boost::shared_ptr<BlackScholesMertonProcess> process_ = boost::dynamic_pointer_cast<BlackScholesMertonProcess>(
				C->createDefaultProcess());
			double state =process_->x0();

			while (t < T){
				state = process_->evolve(t,state,dt,0.0);
				t += dt;
				Assert::AreEqual(state, state);
				
			}
		
		}
	};
}