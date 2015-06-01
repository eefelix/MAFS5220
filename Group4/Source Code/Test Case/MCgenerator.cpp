#include "MCgenerator.h"

using namespace std;
using namespace QuantLib;

//reset the stock price initialized in the constructor
void MCgenerator::set_stock_price(Real new_price)
	{
		S0 = new_price;
	}

//invoke the class "equity_forward" to calculate the NDV of the forward
Real MCgenerator::stage_NPV(Date maturity_date)
	{
		//Construct the dividend TS object, the pointer and the handle of the object
		Handle<YieldTermStructure> DividendTS(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate,div_y,daycount)));

		//Construct the interest rate TS object, the pointer and the handle of the object
		Handle<YieldTermStructure> discountfactor(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate,riskfree,daycount)));
		
		//Construct the option payoff object and the pointer of the object
		boost::shared_ptr<Payoff> forwardpayoff(new ForwardTypePayoff(buysell, k));
		
		//initial the object equity_farward and calculate its NDV
		equity_forward f0(daycount,calen,ModifiedFollowing,selday,forwardpayoff,valdate,maturity_date,discountfactor,DividendTS,S0,spt_dividend);
		return f0.forwardvalue();
	}

//This function has been amended to output a single expected value of the forward NPV for the input maturity date and number of paths
//perform monte carlo simulation on the stock price and invoke the member function "stage_NPV" for calculating its corresponding forward NDV
//vector<Real> MCgenerator::MCgeneration(int samplepath, Date mature_date, BigInteger seed)
Real MCgenerator::MCgeneration(int samplepath, Date mature_date, BigInteger seed)
{
	//Construct the option payoff object
	boost::shared_ptr<Payoff> forwardpayoff(new ForwardTypePayoff(buysell, k));
	
	//construct the stock quote object
	Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(S0)));
	
	//Construct the volatility Term Structure object
	Handle<BlackVolTermStructure> flatVolTS(boost::shared_ptr<BlackVolTermStructure>(new BlackConstantVol(startdate, calen, vol, daycount)));
	
	//Construct the dividend yield term structure object
	Handle<YieldTermStructure> DividendTS(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate,div_y,daycount)));

	//Construct the risk free interest rate term structure object
	Handle<YieldTermStructure> discountfactor(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate,riskfree,daycount)));
	
	//Construct the stock movement process object by using generalized black Scholes Merton Process
	boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(new BlackScholesMertonProcess(underlyingH, DividendTS, discountfactor, flatVolTS));
	
	//generate the ranodm number dW by using the normal distributed Gaussian deviate with average 0.0 andstandard deviation of 1.0
	MersenneTwisterUniformRng unifMt ( seed );
	BoxMullerGaussianRng < MersenneTwisterUniformRng > bmGauss ( unifMt );

	//simulate the stock with time step (1 day)
	Time dt = daycount.yearFraction(startdate,startdate+1) , t = 0;
	Real y = S0; //dummy variable for keeping the initial stock price 
	Real x = S0;
	Real dW;
	Size numVals = (int)(daycount.yearFraction(startdate,mature_date)*365);

	cout << setprecision(10);

	Real accumulator=0;
	//simulation with user input on number of paths  
   	for (int i = 1; i <= samplepath; ++i)
	{
		for (Size j = 1; j <= numVals ; ++j)
		{
			dW = bmGauss.next().value;	     	//generate the random variable
		    x = bsmProcess->evolve(t, x, dt, dW); //E(S(t+dt)) + Sd(t+dt)*dW
			t += dt;
		}
		
		set_stock_price(x);			//reset the stock for calculating the new path of the stock
		accumulator = accumulator + stage_NPV(mature_date); //accumulate the NPVs from the simulated path
		t = 0;
		x = y;		//reset the stock price to initial stock price in the next path's generation
	}

	set_stock_price(y);	//resume S0 to the initial price
	return accumulator/samplepath;    //return the expected value of the forward NDV

}

Real MCgenerator::stage_NPV_V2(Date forwardStart, Date maturity_date)
{
	// =======================================================================================================
	// To calculate the NDV of the forward at future time (Date forwardStart)
	// used for BVA WWR calibration & implementation as its path dependent model
	// =======================================================================================================

	//Construct the dividend TS object, the pointer and the handle of the object
	Handle<YieldTermStructure> DividendTS(boost::shared_ptr<YieldTermStructure>(new FlatForward(forwardStart, div_y, daycount)));

	//Construct the interest rate TS object, the pointer and the handle of the object
	Handle<YieldTermStructure> discountfactor(boost::shared_ptr<YieldTermStructure>(new FlatForward(forwardStart, riskfree, daycount)));

	//Construct the option payoff object and the pointer of the object
	boost::shared_ptr<Payoff> forwardpayoff(new ForwardTypePayoff(buysell, k));

	//initial the object equity_farward and calculate its NDV
	equity_forward f0(daycount, calen, ModifiedFollowing, selday, forwardpayoff, forwardStart, maturity_date, discountfactor, DividendTS, S0, spt_dividend);
	return f0.forwardvalue();
}

map<Date, Real> MCgenerator::SimulatePayoffPath_V1(vector<Date> vSpotDates)
{
	// =======================================================================================================
	// Purpose: To add a new functionality to simulate payoff of equity forward at specified future date
	//			used for BVA WWR calibration & implementation as its path dependent model
	// Input: a vector of specified dates where we want to know the spot price S(ti)
	// Output: a map containing the spot prices results at specified dates 
	// =======================================================================================================

	//Construct the forward payoff object
	boost::shared_ptr<Payoff> forwardpayoff(new ForwardTypePayoff(buysell, k));

	//construct the stock quote object
	Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(S0)));

	//Construct the volatility Term Structure object
	Handle<BlackVolTermStructure> flatVolTS(boost::shared_ptr<BlackVolTermStructure>(new BlackConstantVol(startdate, calen, vol, daycount)));

	//Construct the dividend yield term structure object
	Handle<YieldTermStructure> DividendTS(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate, div_y, daycount)));

	//Construct the risk free interest rate term structure object
	Handle<YieldTermStructure> discountfactor(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate, riskfree, daycount)));

	//Construct the stock movement process object by using generalized black Scholes Merton Process
	boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(new BlackScholesMertonProcess(underlyingH, DividendTS, discountfactor, flatVolTS));

	//generate the ranodm number dW by using the normal distributed Gaussian deviate with average 0.0 and standard deviation of 1.0
	// a time dependent seed was used. 
	unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	MersenneTwisterUniformRng unifMt(seed);
	BoxMullerGaussianRng < MersenneTwisterUniformRng > bmGauss(unifMt);

	Calendar oCalendar = TARGET();
	Time dt, t = 0;
	Real y = S0;	//dummy variable for keeping the initial stock price 
	Real x = S0;
	Real dW;

	cout << setprecision(10);

	map<Date, Real> payoffPath;
	Date currentDate = startdate;
	Date nextDate;

	for (Size j = 0; j < vSpotDates.size(); ++j)
	{
		nextDate = vSpotDates[j];
		dt = Actual365Fixed().yearFraction(currentDate, nextDate);		// Act/365 was used as day count basis
		dW = bmGauss.next().value;		//generate the random variable
		x = bsmProcess->evolve(t, x, dt, dW);		//Browninan drift + diffusion
		set_stock_price(x);		// update the spot price with the new one at nextDate
		payoffPath[nextDate] = stage_NPV_V2(nextDate, vSpotDates.back());		// Find the exposure of EqForward at nextDate. 
		// Assume the EqForward mature matures at 10 years, vSpotDates.back() was the maturity date  
		currentDate = nextDate;
	}

	t = 0;		//reset start time t
	x = y;		//reset the stock price to initial stock price in the next path's generation
	set_stock_price(y);		//resume S0 to the initial price

	return payoffPath;

}

void MCgenerator::SetForwardType(Position::Type forwardType) {
	buysell = forwardType;
}

pair<Real, Real> MCgenerator::MCgeneration_collateral(int samplepath, Date mature_date, BigInteger seed)
{
	//Construct the option payoff object
	boost::shared_ptr<Payoff> forwardpayoff(new ForwardTypePayoff(buysell, k));
	
	//construct the stock quote object
	Handle<Quote> underlyingH(boost::shared_ptr<Quote>(new SimpleQuote(S0)));
	
	//Construct the volatility Term Structure object
	Handle<BlackVolTermStructure> flatVolTS(boost::shared_ptr<BlackVolTermStructure>(new BlackConstantVol(startdate, calen, vol, daycount)));
	
	//Construct the dividend yield term structure object
	Handle<YieldTermStructure> DividendTS(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate,div_y,daycount)));

	//Construct the risk free interest rate term structure object
	Handle<YieldTermStructure> discountfactor(boost::shared_ptr<YieldTermStructure>(new FlatForward(startdate,riskfree,daycount)));
	
	//Construct the stock movement process object by using generalized black Scholes Merton Process
	boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(new BlackScholesMertonProcess(underlyingH, DividendTS, discountfactor, flatVolTS));
	
	//generate the ranodm number dW by using the normal distributed Gaussian deviate with average 0.0 andstandard deviation of 1.0
	MersenneTwisterUniformRng unifMt ( seed );
	BoxMullerGaussianRng < MersenneTwisterUniformRng > bmGauss ( unifMt );

	//simulate the stock with time step (1 day)
	Time dt = daycount.yearFraction(startdate,startdate+1) , t = 0;
	Real y = S0; //dummy variable for keeping the initial stock price 
	Real x = S0;
	Real dW;
	Size numVals = (int)(daycount.yearFraction(startdate,mature_date)*365);
	Real forward_curr_value = 0;
	Real forward_prev_value = stage_NPV(startdate);
	Real collateral_val = 0;
	cout << setprecision(10);
	Date date_in_path = startdate;
	
	//vector<Real> forward_NDV;	//set up the vector for holding the termainal forward NDV value
	Real accumulator=0;
	
	//simulation with user input on number of paths  
   	for (int i = 1; i <= samplepath; ++i)
	{  
		for (Size j = 1; j <= numVals ; ++j)
		{
			dW = bmGauss.next().value;	     	//generate the random variable
		    x = bsmProcess->evolve(t, x, dt, dW); //E(S(t+dt)) + Sd(t+dt)*dW
			set_stock_price(x);
			date_in_path = date_in_path + 1;
			forward_curr_value = stage_NPV(date_in_path);
			
			//Calulate the collateral needed to be posted with interest accrued by either investor or counterparty before one of them had defaulted
			if(j != numVals)
			{   
				if(collateral_val>=0)
					collateral_val = collateral_val*(1+(coll_c_yield/100*(1.0/365))) + (forward_curr_value-forward_prev_value);
				else
					collateral_val = collateral_val*(1+(coll_I_yield/100*(1.0/365))) + (forward_curr_value-forward_prev_value);
			}
			forward_prev_value = forward_curr_value; 
			t += dt;
		} 
		
		set_stock_price(x);			//reset the stock for calculating the new path of the stock
		accumulator = accumulator + stage_NPV(mature_date); //accumulate the NPVs from the simulated path
		t = 0;
		x = y;		//reset the stock price to initial stock price in the next path's generation
		date_in_path = startdate;
	}

	set_stock_price(y);	//resume S0 to the initial price

	//return forward_NDV as well as its corresponding collateral amount in present value
	pair<Real,Real> collateral_forwardpv(accumulator/samplepath,collateral_val*(discountfactor->discount(mature_date))/samplepath);
	return collateral_forwardpv;
}

//Reset the collateral accrued interest posted by investor
void MCgenerator::set_inv_coll_yield(Real new_inv_coll_yield)
{
	coll_I_yield = new_inv_coll_yield;
}

//Reset the collateral accrued interest posted by counterparty
void MCgenerator::set_cpty_coll_yield(Real new_cpty_coll_yield)
{
	coll_c_yield = new_cpty_coll_yield;
}