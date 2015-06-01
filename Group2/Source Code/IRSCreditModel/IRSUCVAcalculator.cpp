#include "pch.h"
#include "IRSUCVAcalculator.hpp"

using namespace QuantLib;

namespace{
	const Volatility constantSwapRateVol = 0.15;
	const Real LossGivenDefault = 1.0;
	const Period fixLegInterval = 1 * Years;
}

IRSUCVAcalculator::IRSInfo::IRSInfo(Rate swapRate, Period maturity, 
		Date startDate, const boost::shared_ptr<YieldTermStructure>& p){
	this->swapRate = swapRate;
	this->maturity = maturity;
	this->startDate = startDate;
	Settings::instance().evaluationDate() = startDate;
	notional = 1.0;
	swapType = VanillaSwap::Receiver;
	
	fixedConvention = Unadjusted;
	fixedFrequency = Annual;
	fixedDayCounter = Thirty360();

	index = boost::shared_ptr<IborIndex>(new Euribor6M(termStructure));
	floatingConvention = ModifiedFollowing;
	floatingFrequency = Semiannual;
	floatingTenor = index->tenor();
	calendar = TARGET();

	termStructure.linkTo(p);

	//swap = MakeVanillaSwap(maturity, index, swapRate)
	//	.withEffectiveDate(startDate)
	//	.withFixedLegTenor(1 * Years)
	//	.withFixedLegDayCount(fixedDayCounter)
	//	.withFloatingLegSpread(0.0)
	//	.withType(swapType);
	Date maturityDate = calendar.advance(startDate, maturity);
	Schedule fixedSchedule(startDate, maturityDate,
		Period(fixedFrequency), calendar, fixedConvention,
		fixedConvention, DateGeneration::Backward, false);
	Schedule floatingSchedule(startDate, maturityDate, 
		Period(floatingFrequency), calendar, floatingConvention,
		floatingConvention, DateGeneration::Backward, false);
	swap = boost::shared_ptr<VanillaSwap>(new VanillaSwap(swapType, notional,
		fixedSchedule, swapRate, fixedDayCounter, 
		floatingSchedule, index, 0.0, index->dayCounter()));
	swap->setPricingEngine(boost::shared_ptr<PricingEngine>(
		new DiscountingSwapEngine(termStructure)));
}

boost::shared_ptr<Swaption> IRSUCVAcalculator::IRSInfo::makeSwaption(const Date& exerciseDate,
	Volatility volatility, Settlement::Type delivery)
{
	Handle<Quote> vol(boost::shared_ptr<Quote>(new SimpleQuote(volatility)));
	boost::shared_ptr<PricingEngine> engine(new BlackSwaptionEngine(termStructure, vol));
	
	boost::shared_ptr<Swaption> swaption(new Swaption(
		swap,
		boost::shared_ptr<Exercise>(new EuropeanExercise(exerciseDate)),
		delivery));
	swaption->setPricingEngine(engine);
	return swaption;
}

boost::shared_ptr<VanillaSwap> IRSUCVAcalculator::getSwap()
{
	return irsInfo->swap;
}

IRSUCVAcalculator& IRSUCVAcalculator::InitializeDefaultCurve(const boost::shared_ptr<DefaultCurve>& defaultcurve)
{
	defaultCurveHandler.linkTo(defaultcurve);
	return *this;
}

Real IRSUCVAcalculator::getAnticUCVA(Rate swapRate, Period maturity, const boost::shared_ptr<YieldTermStructure>& floatingLegYieldCurve)
{
	Real UCVAValue = 0.0;
	defaultCurveHandler->buildCurve(maturity, fixLegInterval);
	irsInfo = boost::shared_ptr<IRSInfo>(new IRSInfo(swapRate, maturity, 
											defaultCurveHandler->getStartDate(), floatingLegYieldCurve));

	Date endDate = defaultCurveHandler->getStartDate() + maturity;
	std::vector<Date> cpCreditTermDate = defaultCurveHandler->getDefaultTerm();
	Settings::instance().evaluationDate() = cpCreditTermDate.front();
	for (size_t i = 0; cpCreditTermDate[i] < endDate; i++)
	{
		Probability defaultProb = defaultCurveHandler->getSurvivalProb(i) - defaultCurveHandler->getSurvivalProb(i+1);
		boost::shared_ptr<Swaption> swaption = irsInfo->makeSwaption(cpCreditTermDate[i],
			constantSwapRateVol);
		Real swaptionValue = swaption->NPV();
		UCVAValue += LossGivenDefault * defaultProb * swaptionValue;
	}
	return UCVAValue;
}

Real IRSUCVAcalculator::getPostpUCVA(Rate swapRate, Period maturity)
{
	//Real endDate = dc.yearFraction(startDate, startDate + maturity);
	return 0.0;
}