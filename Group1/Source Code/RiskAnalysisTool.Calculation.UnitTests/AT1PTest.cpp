#include <pch.h>
#include <limits>
#include <boost/shared_ptr.hpp>
#include <ql/handle.hpp>
#include <ql/time/date.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/termstructures/credit/defaultprobabilityhelpers.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/dategenerationrule.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/quote.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/models/defaultcdshelper.hpp>

using namespace std;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace QuantLib;

namespace QuantLib
{
	TEST_CLASS(AT1PTest)
	{

	public:
		/*
		* Three test case from chapter three of Brigo's book
		* page 59-60 table 3.1 3.2 3.3
		* CDS maturities are 1y,3y,5y,7y,10y
		*/

		TEST_METHOD(AT1PTest1){
			boost::shared_ptr<vector<double>> CDSspreads(new vector<double>({ 0.0016, 0.0029, 0.0045, 0.0050, 0.0058 }));
			Date referenceDate = Date(10, July, 2007);
			boost::shared_ptr<vector<QuantLib::Period>> tenors(new vector<QuantLib::Period>(
			{ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			Settings::instance().evaluationDate() = referenceDate;
			Calendar calendar = TARGET();
			DayCounter daycounter = Actual360();
			double interestrate = 0.05;
			boost::shared_ptr<YieldTermStructure> discountCurve_(new
				FlatForward(referenceDate, interestrate, daycounter));
			boost::shared_ptr<AT1Pmodel> model = boost::shared_ptr<AT1Pmodel>(
				new AT1Pmodel(tenors, discountCurve_));

			std::vector<boost::shared_ptr<QuantLib::DefaultProbabilityHelper> > instruments;
			for (int i = 0; i < 5; i++) {
				instruments.push_back(
					boost::shared_ptr<QuantLib::DefaultProbabilityHelper>(
					new QuantLib::SpreadCdsHelper(
					QuantLib::Handle<QuantLib::Quote>(
					boost::shared_ptr<QuantLib::Quote>(
					new QuantLib::SimpleQuote(CDSspreads->at(i)))),
					tenors->at(i),
					0,
					TARGET(),
					Quarterly,
					ModifiedFollowing,
					DateGeneration::TwentiethIMM,
					Actual360(),
					0.4,
					Handle<YieldTermStructure>(discountCurve_))));
			}

			boost::shared_ptr < QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> >
				hazardRateStructure_ =
				boost::make_shared<QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat>>(
				referenceDate,
				instruments,
				Actual360());

			std::vector<Time> time(tenors->size());

			for (auto itr = tenors->begin(); itr != tenors->end(); ++itr)
			{
				time[itr - tenors->begin()] = Actual360().yearFraction(referenceDate, referenceDate + *itr);
			}

			std::vector<boost::shared_ptr<CalibrationHelper>> cdsHelper_;
			cdsHelper_.reserve(time.size());

			for (auto itr = time.begin(); itr != time.end(); ++itr){
				cdsHelper_.push_back(
					boost::make_shared<DefaultCdsHelper>(*itr,
					Handle<Quote>(new SimpleQuote(hazardRateStructure_->defaultProbability(*itr, true))),
					Handle<YieldTermStructure>(discountCurve_), model));
			}

			model->calibrate(cdsHelper_, Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			boost::shared_ptr<BlackScholesMertonProcess> BSMprocess = boost::dynamic_pointer_cast<BlackScholesMertonProcess>(
				model->process());

			vector<double> vol(5);
			for (int i = 0; i < 5; ++i){
				vol[i] = BSMprocess->blackVolatility()->blackVol(Actual360().yearFraction(referenceDate, referenceDate + tenors->at(i)), 1);
			}
			boost::shared_ptr<QuantLib::Array> actual(new Array(vol.begin(), vol.end()));
			int n = actual->size();

			std::vector<double> expected = { 0.292, 0.140, 0.145, 0.12, 0.127 };

			std::vector<double> survivalProbabilityexpected = { 0.997, 0.985, 0.961, 0.941, 0.902 };
			std::vector<double> survivalProbabilityactual;
			for (int i = 0; i < n; i++){
				survivalProbabilityactual.push_back(1 - model->defaultProbability(
					(Actual360().yearFraction(referenceDate, referenceDate + tenors->at(i)))));
			}
			for (int i = 0; i < n; i++){
				Assert::AreEqual(survivalProbabilityexpected[i], survivalProbabilityactual[i], survivalProbabilityexpected[i] * 0.05);
			}
		}

		TEST_METHOD(AT1PTest2){
			boost::shared_ptr<vector<double>> CDSspreads(new vector<double>({ 0.0397, 0.0315, 0.0277, 0.0258, 0.0240 }));
			Date referenceDate = Date(12, June, 2008);
			boost::shared_ptr<vector<QuantLib::Period>> tenors(new vector<QuantLib::Period>(
			{ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			Settings::instance().evaluationDate() = referenceDate;
			Calendar calendar = TARGET();
			DayCounter daycounter = Actual360();
			double interestrate = 0.05;
			boost::shared_ptr<YieldTermStructure> discountCurve_(new
				FlatForward(referenceDate, interestrate, daycounter));
			boost::shared_ptr<AT1Pmodel> model = boost::shared_ptr<AT1Pmodel>(
				new AT1Pmodel(tenors, discountCurve_));

			std::vector<boost::shared_ptr<QuantLib::DefaultProbabilityHelper> > instruments;
			for (int i = 0; i < 5; i++) {
				instruments.push_back(
					boost::shared_ptr<QuantLib::DefaultProbabilityHelper>(
					new QuantLib::SpreadCdsHelper(
					QuantLib::Handle<QuantLib::Quote>(
					boost::shared_ptr<QuantLib::Quote>(
					new QuantLib::SimpleQuote(CDSspreads->at(i)))),
					tenors->at(i),
					0,
					TARGET(),
					Quarterly,
					ModifiedFollowing,
					DateGeneration::TwentiethIMM,
					Actual360(),
					0.4,
					Handle<YieldTermStructure>(discountCurve_))));
			}

			boost::shared_ptr < QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> >
				hazardRateStructure_ =
				boost::make_shared<QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat>>(
				referenceDate,
				instruments,
				Actual360());

			std::vector<Time> time(tenors->size());

			for (auto itr = tenors->begin(); itr != tenors->end(); ++itr)
			{
				time[itr - tenors->begin()] = Actual360().yearFraction(referenceDate, referenceDate + *itr);
			}

			std::vector<boost::shared_ptr<CalibrationHelper>> cdsHelper_;
			cdsHelper_.reserve(time.size());

			for (auto itr = time.begin(); itr != time.end(); ++itr){
				cdsHelper_.push_back(
					boost::make_shared<DefaultCdsHelper>(*itr,
					Handle<Quote>(new SimpleQuote(hazardRateStructure_->defaultProbability(*itr, true))),
					Handle<YieldTermStructure>(discountCurve_), model));
			}

			model->calibrate(cdsHelper_, Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			boost::shared_ptr<BlackScholesMertonProcess> BSMprocess = boost::dynamic_pointer_cast<BlackScholesMertonProcess>(
				model->process());

			vector<double> vol(5);
			for (int i = 0; i < 5; ++i){
				vol[i] = BSMprocess->blackVolatility()->blackVol(Actual360().yearFraction(referenceDate, referenceDate + tenors->at(i)), 1);
			}
			boost::shared_ptr<QuantLib::Array> actual(new Array(vol.begin(), vol.end()));
			int n = actual->size();

			std::vector<double> expected = { 0.45, 0.219, 0.186, 0.181, 0.175 };

			std::vector<double> survivalProbabilityexpected = { 0.936, 0.857, 0.800, 0.751, 0.688 };
			std::vector<double> survivalProbabilityactual;
			for (int i = 0; i < n; i++){
				survivalProbabilityactual.push_back(1 - model->defaultProbability(
					(Actual360().yearFraction(referenceDate, referenceDate + tenors->at(i)))));
			}
			for (int i = 0; i < n; i++){
				Assert::AreEqual(survivalProbabilityexpected[i], survivalProbabilityactual[i], survivalProbabilityexpected[i] * 0.05);
			}
		}

		TEST_METHOD(AT1PTest3){
			boost::shared_ptr<vector<double>> CDSspreads(new vector<double>({ 0.1437, 0.0902, 0.0710, 0.0636, 0.0588 }));
			Date referenceDate = Date(12, Sep, 2008);
			boost::shared_ptr<vector<QuantLib::Period>> tenors(new vector<QuantLib::Period>(
			{ Period(1, Years), Period(3, Years), Period(5, Years), Period(7, Years), Period(10, Years) }));

			Settings::instance().evaluationDate() = referenceDate;
			Calendar calendar = TARGET();
			DayCounter daycounter = Actual360();
			double interestrate = 0.05;
			boost::shared_ptr<YieldTermStructure> discountCurve_(new
				FlatForward(referenceDate, interestrate, daycounter));
			boost::shared_ptr<AT1Pmodel> model = boost::shared_ptr<AT1Pmodel>(
				new AT1Pmodel(tenors, discountCurve_));

			std::vector<boost::shared_ptr<QuantLib::DefaultProbabilityHelper> > instruments;
			for (int i = 0; i < 5; i++) {
				instruments.push_back(
					boost::shared_ptr<QuantLib::DefaultProbabilityHelper>(
					new QuantLib::SpreadCdsHelper(
					QuantLib::Handle<QuantLib::Quote>(
					boost::shared_ptr<QuantLib::Quote>(
					new QuantLib::SimpleQuote(CDSspreads->at(i)))),
					tenors->at(i),
					0,
					TARGET(),
					Quarterly,
					ModifiedFollowing,
					DateGeneration::TwentiethIMM,
					Actual360(),
					0.4,
					Handle<YieldTermStructure>(discountCurve_))));
			}

			boost::shared_ptr < QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat> >
				hazardRateStructure_ =
				boost::make_shared<QuantLib::PiecewiseDefaultCurve<QuantLib::HazardRate, QuantLib::BackwardFlat>>(
				referenceDate,
				instruments,
				Actual360());

			std::vector<Time> time(tenors->size());

			for (auto itr = tenors->begin(); itr != tenors->end(); ++itr)
			{
				time[itr - tenors->begin()] = Actual360().yearFraction(referenceDate, referenceDate + *itr);
			}

			std::vector<boost::shared_ptr<CalibrationHelper>> cdsHelper_;
			cdsHelper_.reserve(time.size());

			for (auto itr = time.begin(); itr != time.end(); ++itr){
				cdsHelper_.push_back(
					boost::make_shared<DefaultCdsHelper>(*itr,
					Handle<Quote>(new SimpleQuote(hazardRateStructure_->defaultProbability(*itr, true))),
					Handle<YieldTermStructure>(discountCurve_), model));
			}

			model->calibrate(cdsHelper_, Simplex(0.001), EndCriteria(10000, 50, 1e-8, 1e-8, 1e-8));

			boost::shared_ptr<BlackScholesMertonProcess> BSMprocess = boost::dynamic_pointer_cast<BlackScholesMertonProcess>(
				model->process());

			vector<double> vol(5);
			for (int i = 0; i < 5; ++i){
				vol[i] = BSMprocess->blackVolatility()->blackVol(Actual360().yearFraction(referenceDate, referenceDate + tenors->at(i)), 1);
			}
			boost::shared_ptr<QuantLib::Array> actual(new Array(vol.begin(), vol.end()));
			int n = actual->size();

			std::vector<double> expected = { 0.622, 0.308, 0.243, 0.269, 0.295 };
			std::vector<double> survivalProbabilityexpected = { 0.784, 0.655, 0.591, 0.525, 0.434 };
			std::vector<double> survivalProbabilityactual;
			for (int i = 0; i < n; i++){
				survivalProbabilityactual.push_back(1 - model->defaultProbability(
					(Actual360().yearFraction(referenceDate, referenceDate + tenors->at(i)))));
			}
			for (int i = 0; i < n; i++){
				Assert::AreEqual(survivalProbabilityexpected[i], survivalProbabilityactual[i], survivalProbabilityexpected[i] * 0.05);
			}
		}
	};
}