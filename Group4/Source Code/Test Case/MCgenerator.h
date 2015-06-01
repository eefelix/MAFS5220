#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <chrono>

#ifndef MCGENERATOR_H		// to avoid error on multiple declaration
#define MCGENERATOR_H		// to avoid error on multiple declaration

#include "equity_forward.h"

//This class is used for simulating the stock price path by monte carlo simulation and calculate its corresponding forward NDV 
class MCgenerator
{
public:
	MCgenerator(Date start, Date mat, Date val,	Natural settle,	Calendar cal, Real init_S0, DayCounter dayCount, Real strike,
		Spread div_q, Real spt_div,	Rate risk_free,	Volatility si,
	Position::Type type): startdate(start), matdate(mat), valdate(val), selday(settle), 
	S0(init_S0), k(strike), div_y(div_q),spt_dividend(spt_div),riskfree(risk_free),calen(cal),daycount(dayCount), vol(si),buysell(type)
	{
	}		//Constructor
	
	MCgenerator(Date start, Date mat, Date val,	Natural settle,	Calendar cal, Real init_S0, DayCounter dayCount, Real strike,
		Spread div_q, Real spt_div,	Rate risk_free,	Real Inv_coll_yield, Real Inv_cpty_yield, Volatility si,
	Position::Type type): startdate(start), matdate(mat), valdate(val), selday(settle), 
	S0(init_S0), k(strike), div_y(div_q),spt_dividend(spt_div),riskfree(risk_free),calen(cal),daycount(dayCount), vol(si), buysell(type),
	coll_I_yield(Inv_coll_yield), coll_c_yield(Inv_cpty_yield)
	{
	}	

	//member functions
	//reset the stock price initialized in the constructor
	void set_stock_price(Real new_price);

	//invoke the class "equity_forward" to calculate the NDV of the forward
	Real stage_NPV(Date maturity_date);

	//perform monte carlo simulation on the stock price and its corresponding forward NDV
	//vector<Real> MCgeneration(int samplepath, Date mature_date, BigInteger seed);
	Real MCgeneration(int samplepath, Date mature_date, BigInteger seed);

	// To calculate the NDV of the forward at future time (Date forwardStart)
	// used for BVA WWR calibration & implementation as its path dependent model
	Real stage_NPV_V2(Date forwardStart, Date maturity_date);

	// Perform monte carlo simulation to simulate payoff of equity forward at future time
	// Return map storing the payoff value of equity forward at specified vSpotDates
	// used for BVA WWR calibration & implementation as its path dependent model
	map<Date, Real> SimulatePayoffPath_V1(vector<Date> vSpotDates);

	// used for BVA WWR calibration & implementation as its path dependent model
	void SetForwardType(Position::Type forwardType);

	//This is the monte carlo simulation with the inclusion of collateral
	pair<Real, Real> MCgeneration_collateral(int samplepath, Date mature_date, BigInteger seed);

	//Reset the collateral accrued interest posted by investor
	void set_inv_coll_yield(Real new_inv_coll_yield);

	//Reset the collateral accrued interest posted by counterparty
	void set_cpty_coll_yield(Real new_cpty_coll_yield);
private:
	Date startdate, valdate, matdate;
	Natural selday;
	Calendar calen;
	DayCounter daycount;
	Real S0, k, spt_dividend, riskfree;
	Real coll_I_yield, coll_c_yield;
	Spread div_y;
	Volatility vol;
	Position::Type buysell;
};

#endif /*MCGENERATOR_H*/
