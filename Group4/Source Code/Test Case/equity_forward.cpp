#include "equity_forward.h"

//return the underlying spot price of stock
Real equity_forward::spotValue() const {
		return underlyingSpotValue_; 
	}

//return the present value of the dividend of the stock
Real equity_forward::spotIncome(const Handle<YieldTermStructure> &incomeDiscountCurve) const {	
		return underlyingIncome_/(incomeDiscountCurve->discount(maturityDate_));
	}

//perform the calculation of the forward value and return the present value of the forward
Real equity_forward::forwardvalue() const {
		return ((payoff_->operator()(spotprice-spotIncome_))* discountCurve_->discount(maturityDate_));
	}

