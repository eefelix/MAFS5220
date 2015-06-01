#pragma warning(disable:4819)
#include <ql/quantlib.hpp>

using namespace std;
using namespace QuantLib;

//Inherited from base class Forward
class equity_forward: public Forward
{
public:
	equity_forward(const DayCounter& dayCounter_eq,			//constructor for the class
                const Calendar& calendar_eq,
                BusinessDayConvention businessDayConvention_eq,
                Natural settlementDays_eq,
                const boost::shared_ptr<Payoff>& payoff_eq,
                const Date& valueDate_eq,
                const Date& maturityDate_eq,
                const Handle<YieldTermStructure>& discountCurve_eq,
				const Handle<YieldTermStructure>& dividendCurve_eq,
				Real spot, Real spot_dividend)
				: Forward(dayCounter_eq,calendar_eq,businessDayConvention_eq,settlementDays_eq,
				payoff_eq,valueDate_eq,maturityDate_eq,discountCurve_eq){
	
				underlyingSpotValue_=spot;
				underlyingIncome_= spot_dividend;  
				discount_factor = discountCurve_eq->discount(maturityDate_eq);
				dividend_factor = dividendCurve_eq->discount(maturityDate_eq);
				spotprice = spot/discount_factor;
				spotIncome_= spot_dividend/dividend_factor; 
	} 
	
	//member functions
	//return the data member spot value input
	Real spotValue() const;
	
	//return the present value of the dividend
	Real spotIncome(const Handle<YieldTermStructure> &incomeDiscountCurve) const;

	//return the present value of the forward
	Real forwardvalue() const;


private:
	Real spotprice;
	Real spotIncome_;
	Real discount_factor;
	Real dividend_factor;
	
};