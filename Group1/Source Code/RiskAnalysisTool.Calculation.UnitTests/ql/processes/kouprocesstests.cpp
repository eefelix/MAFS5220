#include <pch.h>
#include <ql/quantlib.hpp>
#include <boost/shared_ptr.hpp>
#include <ql/processes/kouprocess.hpp>
#include <ql/math/randomnumbers/fakeuniformrng.hpp>

using namespace std;
using namespace QuantLib;
using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace QuantLib 
{
	TEST_CLASS(KouProcessTests)
	{
	private:
		typedef GeneralizedKouProcess<FakeUniformRng<Real>> FakeKouProcess;

		boost::shared_ptr<FakeKouProcess> CreateProcess(
			const Calendar &cal, const DayCounter &dc, const Date &d0,
			Real x0, Real q, Real r, Real vol, const Real jumpIntensity,
			const Real posProbability, const Real posJumpMean, const Real negJumpMean,
			const vector<Real> &seq1, const vector<Real> &seq2)
		{
			return boost::make_shared<FakeKouProcess>(
				Handle<Quote>(boost::make_shared<SimpleQuote>(x0)),
				Handle<YieldTermStructure>(boost::make_shared<FlatForward>(d0, q, dc)),
				Handle<YieldTermStructure>(boost::make_shared<FlatForward>(d0, r, dc)),
				Handle<BlackVolTermStructure>(boost::make_shared<BlackConstantVol>(d0, cal, vol, dc)),
				jumpIntensity,
				posProbability,
				posJumpMean,
				negJumpMean,
				FakeUniformRng<Real>(seq1),
				FakeUniformRng<Real>(seq2)
				);
		}
	public:
		TEST_METHOD(ConstructorTest)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 0.05;
			const Real negJumpMean = -0.06;

			std::vector<Real> seq1 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };
			std::vector<Real> seq2 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

			boost::shared_ptr<FakeKouProcess> process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean,
				seq1, seq2
				);

			Assert::IsNotNull(process.get());

			Assert::AreEqual(100.0, process->x0());
			Assert::AreEqual(jumpIntensity, process->jumpIntensity());
			Assert::AreEqual(posProb, process->posProbability());
			Assert::AreEqual(posJumpMean, process->posJumpMean());
			Assert::AreEqual(negJumpMean, process->negJumpMean());
		}

		TEST_METHOD(SamplePathTest1)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 0.05;
			const Real negJumpMean = -0.06;

			std::vector<Real> seq1 = { 0.1, 0.2, 0.9, 0.4, 0.5, 0.6 };
			std::vector<Real> seq2 = { 0.1, 0.8, 0.3, 0.4, 0.5, 0.6 };
			std::vector<Real> dw = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 };

			boost::shared_ptr<FakeKouProcess> process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean,
				seq1, seq2
				);

			Assert::IsNotNull(process.get());

			vector<Real> expected = { 97.391, 98.370, 99.359, 100.357, 101.366, 102.385 };

			for (int i = 0; i < seq1.size(); ++i) {
				Real dx = process->evolve(0.0, x0, 1.0, dw[i]);
				Assert::AreEqual(expected[i], dx, 0.001);
			}

		}

		TEST_METHOD(SamplePathTest2)
		{
			Calendar cal;
			DayCounter dc;
			Date d0;

			cal = TARGET();
			dc = Business252();
			d0 = Date(2, Jan, 2013);

			const Real x0 = 100.0;
			const Real r = 0.03;
			const Real q = 0.05;
			const Real vol = 0.1;
			const Real jumpIntensity = 0.2;
			const Real posProb = 0.4;
			const Real posJumpMean = 0.05;
			const Real negJumpMean = -0.06;

			std::vector<Real> seq1 = { 0.3, 0.2, 0.3, 0.4, 0.9, 0.5 };
			std::vector<Real> seq2 = { 0.3, 0.2, 0.8, 0.4, 0.5, 0.4 };
			std::vector<Real> dw = { 0.2, 0.2, 0.6, 0.9, 0.1, 0.2 };

			boost::shared_ptr<FakeKouProcess> process = CreateProcess(
				cal, dc, d0, x0, q, r, vol,
				jumpIntensity, posProb, posJumpMean, negJumpMean,
				seq1, seq2
				);

			Assert::IsNotNull(process.get());

			vector<Real> expected = { 98.370, 98.370, 102.385, 105.503, 97.391, 98.370 };

			for (int i = 0; i < seq1.size(); ++i) {
				Real dx = process->evolve(0.0, x0, 1.0, dw[i]);
				Assert::AreEqual(expected[i], dx, 0.001);
			}
		}
	};
}