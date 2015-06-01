
#ifndef irsucvacalculator_hpp
#define irsucvacalculator_hpp

#include <ql/types.hpp>
#include <ql/time/date.hpp>
#include "DefaultCurve.hpp"
//#include <ql/time/businessdayconvention.hpp>
//#include <ql/time/calendar.hpp>

#include <vector>
using namespace QuantLib;

/**
* Class Name: IRSUCVAcalculator.
* Usage:
*	1. Unilateral-CVA calculator for interest rate swap
*	2. Without considering wrong way risk
*	3. Using the Vanilla swap class from QuantLib to get the NPV
*/

class IRSUCVAcalculator{
public:
	class IRSInfo;

	IRSUCVAcalculator(){}
	~IRSUCVAcalculator(){}

	IRSUCVAcalculator& InitializeDefaultCurve(const boost::shared_ptr<DefaultCurve>& defaultcurve);
	/**
	* Function Name: getAnticUCVA, getPostpUCVA
	* Usage: calculate the UCVA of anticipated and postponed approximation
	*/
	Real getAnticUCVA(Rate swapRate, Period maturity, const boost::shared_ptr<YieldTermStructure>& floatingLegYieldCurve);
	Real getPostpUCVA(Rate swapRate, Period maturity);
	boost::shared_ptr<VanillaSwap> getSwap();

private:
	RelinkableHandle<DefaultCurve> defaultCurveHandler;
	boost::shared_ptr<IRSInfo> irsInfo;
};


/**
* Class Name: IRSInfo.
* Usage:
*	1. Inner class of IRSUCVAcalculator
*	2. Contain the information of vanilla swap of different tenors
*	3. Useful for generating corresponding swaptions
*/

class IRSUCVAcalculator::IRSInfo{
public:
	IRSInfo(Rate swapRate, Period maturity, Date startDate, const boost::shared_ptr<YieldTermStructure>& p);
	boost::shared_ptr<Swaption> makeSwaption(const Date& exerciseDate,
		Volatility volatility, Settlement::Type delivery = Settlement::Physical);
	boost::shared_ptr<VanillaSwap> swap;
private:
	Date startDate;
	Period maturity;
	Real notional;
	Real swapRate;
	VanillaSwap::Type swapType;

	Calendar calendar;
	BusinessDayConvention fixedConvention;
	Frequency fixedFrequency;
	Frequency floatingFrequency;
	DayCounter fixedDayCounter;

	BusinessDayConvention floatingConvention;
	Period floatingTenor;
	boost::shared_ptr<IborIndex> index;
	RelinkableHandle<YieldTermStructure> termStructure;
};

#endif