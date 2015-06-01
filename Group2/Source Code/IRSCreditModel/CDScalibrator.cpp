#include "pch.h"
#include "CDScalibrator.hpp"

using namespace QuantLib;

void 
CDSCalibrator::extractDefaultProbabilityCurve()
{
	Calendar calendar = TARGET();
	startDate = calendar.adjust(startDate);
	Settings::instance().evaluationDate() = startDate;

	boost::shared_ptr<Quote> flatInterestRate(new SimpleQuote(interestRate));
	Handle<YieldTermStructure> tsCurve(boost::shared_ptr<FlatForward>(
		new FlatForward(startDate, Handle<Quote>(flatInterestRate), Actual360())));

	std::vector<boost::shared_ptr<DefaultProbabilityHelper>> instruments;
	for (size_t i = 0; i < Quotes.size(); i++) {
		instruments.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
			new SpreadCdsHelper(
			Handle<Quote>(boost::shared_ptr<Quote>(new SimpleQuote(Quotes[i]))),
			Tenors[i],
			0,
			calendar,
			Quarterly,
			Following,
			DateGeneration::TwentiethIMM,
			Actual360(),
			recoveryRate,
			tsCurve)));
	}

	boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat>> 
		impliedSurvivalProbability = boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat>>
			(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
			startDate, instruments, Actual360()
			));
	impliedSurvivalProbability->nodes();
	for (size_t i = 0; i < Tenors.size(); i++) {
		impliedProbs.push_back(impliedSurvivalProbability->survivalProbability(startDate + Tenors[i]));
	}
}


Probability 
CDSCalibrator::FinalizedSurvivalProbability(Date testDate)
{
	return at1p->SurvivalProbability(testDate);
}

Volatility 
CDSCalibrator::FinalizedVolatility(Date testDate)
{
	QL_REQUIRE(testDate >= startDate && testDate <= startDate + Tenors.back(), "Invalid date!");
	for (size_t ind = 0; ind != Tenors.size(); ind++)
	{
		if (testDate <= startDate + Tenors[ind])
			return segmentedVol[ind];
	}
	return 0.0;
}

void 
CDSCalibrator::StartCalibrate()
{
	Size maxIterations = 1000;
	Size minStatIterations = 100;
	Real rootEpsilon = 1e-12;
	Real functionEpsilon = 1e-12;
	Real gradientNormEpsilon = 1e-12;
	EndCriteria myEndCrit(maxIterations, minStatIterations, rootEpsilon,
		functionEpsilon, gradientNormEpsilon);

	int termNum = Tenors.size();
	Array startVol(termNum, 0.3);
	BoundaryConstraint constraint(0.0, 0.8);

	std::vector<Date> dateVector;
	dateVector.push_back(startDate);
	for (size_t i = 0; i < Tenors.size(); i++) {
		dateVector.push_back(startDate + Tenors[i]);
	}
	optimizeFunc.setTermStructure(dateVector);
	optimizeFunc.setCdsQuotes(Quotes);
	optimizeFunc.setImpliedSurvivalProb(impliedProbs);

	Problem CDScalibrationProblem(optimizeFunc, constraint, startVol);
	Simplex solver(0.1);
	EndCriteria::Type solvedCrit = solver.minimize(CDScalibrationProblem, myEndCrit);

	Array resultedVol = CDScalibrationProblem.currentValue();
	segmentedVol.clear();
	for (Size i = 0; i < resultedVol.size(); i++){
		segmentedVol.push_back(resultedVol[i]);
	}
	at1p->Initialization(dateVector, segmentedVol);
}


Real 
CDSCalibrateProblemFunction::value(const Array& x) const
{
	std::vector<Volatility> segmentedVol;
	for (Size i = 0; i < x.size(); i++){
		segmentedVol.push_back(x[i]);
	}
	at1p->Initialization(CDSTermStructure, segmentedVol);
	Real sum = 0.0;
	for (int i = 0; i < termNum; i++){
		Real tmp = at1p->SurvivalProbability(CDSTermStructure[i + 1]);
		sum += (tmp - impliedProbability[i]) * (tmp - impliedProbability[i]);
	}
	return sum;
}

Disposable<Array> 
CDSCalibrateProblemFunction::values(const Array& x) const
{
	Array res(1);
	res[0] = value(x);
	return res;
}


//Real CDSCalibrateProblemFunction::CDSValueFormula(Spread premium, Date startDate, Date endDate) const
//{
//	Real partResult = 0.0;
//	Time stepWidth = dc.yearFraction(startDate, startDate + 1 * Days);
//	for (Date t = startDate; t < endDate; t++) {
//		if (!calendar.isBusinessDay(t))
//			continue;
//		partResult += stepWidth * IntergratedCallFunction(premium, t);
//	}
//
//	for (std::vector<Date>::const_iterator iterDate = ++CDSPremiumPaymentLegDate.cbegin();
//		iterDate != CDSPremiumPaymentLegDate.cend(); iterDate++)
//	{
//		if (*iterDate <= endDate) {
//			Probability survivalProb = at1p->SurvivalProbability(*iterDate);
//			BigInteger businessDay = calendar.businessDaysBetween(*(iterDate - 1), *iterDate);
//			partResult -= premium * dc.yearFraction(*(iterDate - 1), *iterDate)
//			partResult -= premium * stepWidth * businessDay
//				* exp(-RiskFreeRate * dc.yearFraction(startDate, *iterDate)) * survivalProb;
//		}
//		else {
//			break;
//		}
//	}
//
//	return partResult;
//}
//
//Real CDSCalibrateProblemFunction::differentialization(Date testDate) const
//{
//	Real differential;
//	Time stepWidth = dc.yearFraction(startDate, startDate + 1 * Days);
//	if (testDate == CDSTermDate.front()) {
//		differential = -(at1p->SurvivalProbability(testDate + 1 * Days) - at1p->SurvivalProbability(testDate)) / stepWidth;
//	}
//	else if (testDate == CDSTermDate.back()) {
//		differential = -(at1p->SurvivalProbability(testDate) - at1p->SurvivalProbability(testDate - 1 * Days)) / stepWidth;
//	}
//	else {
//		differential = -(at1p->SurvivalProbability(testDate + 1 * Days) - at1p->SurvivalProbability(testDate - 1 * Days)) / (stepWidth * 2);
//	}
//	return differential;
//}
//
//Real CDSCalibrateProblemFunction::IntergratedCallFunction(Spread premium, Date currDate) const
//{
//	Time timeInterval = dc.yearFraction(startDate, currDate);
//	DiscountFactor df = exp(-RiskFreeRate *  timeInterval);
//	Real differential = differentialization(currDate);
//	Date latestPrevPaymentDate;
//	for (std::vector<Date>::const_iterator iterDate = ++CDSPremiumPaymentLegDate.cbegin();
//		iterDate != CDSPremiumPaymentLegDate.cend(); iterDate++)
//	{
//
//		if (currDate <= *iterDate) {
//			latestPrevPaymentDate = *(iterDate - 1);
//			break;
//		}
//	}
//	BigInteger accrual = calendar.businessDaysBetween(latestPrevPaymentDate, currDate, true, false);
//	return LossGivenDefaultRate * df * differential -
//		premium * df * dc.yearFraction(latestPrevPaymentDate, currDate) * differential;
//		premium * df * accrual * dc.yearFraction(startDate + 1*Days, startDate) * differential;
//}
