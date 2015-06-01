#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <fstream>
#include <map>
#include <chrono>

using namespace std;
using namespace QuantLib;


//This is for storing the volatility term structures from the calibration
class VolatilityTS{
private:
	int i;
	Schedule oSchedule_;
	DayCounter dc;
	Date valueDate_;

public:
		
	int oSize;
	std::vector<Date> oDates_;
	map<Date, Real> tsVolatility;
	
	//constructor
	VolatilityTS();
	//constructor 2
	VolatilityTS(Date valueDate, std::vector<Date> oDates);

	//update the volatility term structure each time after the calibration
	void UpdateStepwiseVolatility(const Array& x, int iTenor);

	//calculate the integration of the square of calibrated volatility
	Real GetSumOfVolSquare(Date endDate);
};

//This is for the calculating the survivual probabilities from the calibrated volatility
class AT1P {

private:
	Real V0_;
	Real H0_;
	Real B0_;
	Real d1;
	Real d2;
	Real oQT;
	Date valueDate_;
	vector<Date> survDates_;
	VolatilityTS sigma_;

public:
	AT1P();	//constructor

	AT1P(Real V0, Real H0, Real B0, vector<Date> survDates, Date valueDate); //constructor 2

	//calculate the survival probability
	Real GetSurvivalProbability(Date endDate);

	VolatilityTS get_sigma();
	
	//update the volatility term structure for calculating the survivual probability
	void UpdateVolatilityTS(VolatilityTS oVolatilityTS);
};


class CreditDefaultSwapV2{

private:
	Calendar calendar;
	Protection::Side side_;
	Real notional_;
	Rate spread_;
	Schedule schedule_;
	const Date& protectionStart_;
	AT1P oAT1P_;
	Real lossGivenDefault;
	

public:
	//constructor 1
	CreditDefaultSwapV2(Protection::Side side,
		Real notional,
		const Date& protectionStart);
	
	//constructor 1
	CreditDefaultSwapV2(Protection::Side side,
		Real notional,
		Rate spread,
		Schedule schedule,
		const Date& protectionStart);
	Real risk_free;
	
	//reset the CDS premium
	void SetPremium(Rate premium);
	
	//reset the schedule
	void SetSchedule(Schedule schedule);

	//set the survivual probabilty for calculating the CDS NPV
	void SetAT1PModel(AT1P oAT1P);
	
	//update the risk free rate during the calculation of the CDS NPV
	void SetDis_fac(Real risk_freerate);

	Real pow(Real x, int power);
	
	//return the term's discount factor
	Real discount_factor(Date discount_date);

	//calculating the CDS default leg
	Real DefaultLegPV();

	//calculating the CDS coupon leg
	Real CouponLegPV();

	//calculate the CDS NPV
	Real NPV();

};

class FunctionCdsPv : public CostFunction{

private:
	VolatilityTS oVolTS_1;
	CreditDefaultSwapV2 oCDS_1;
	AT1P oAT1P_1;
	vector<Period> survTenor_;
	vector<Rate> cdsPremiums_;
	Date valueDate_;
	int iTenor_;

public:
	//constructor 1
	FunctionCdsPv(VolatilityTS oVolTS, CreditDefaultSwapV2 oCDS, AT1P oAT1P);
	//constructor 2
	FunctionCdsPv(VolatilityTS oVolTS, CreditDefaultSwapV2 oCDS, AT1P oAT1P, vector<Period> survTenor, vector<Rate> cdsPremiums, Date valueDate);

	//reset the calibration tenor
	void SetTenorToCalibrate(int iTenor);

	//reset the premium for calibration
	void setcds_premium(Rate premium);

	//update the volatility term structure for calculating the CDS NPV
	void UpdateVolatilityTS(VolatilityTS oVolTS);
	
	// Override the cost function for the simplex optimization 
	Real value(const Array& x) const;

	Disposable<Array> values(const Array& x) const;

};

class At1pCalibration
{
private:
	Date valueDate;
	vector<Rate> cdsPremiums;
	vector<Date> survDate;
	vector<Period> survTenor;
	vector<Real> discount_rate;
	Calendar calendar;
	vector<Real> at1p_prob;
	vector<Real> calibrated_vol;
	VolatilityTS oVolTS;
	AT1P oAT1P;
	Real UniformDistributionGenerator(double range_from, double range_to);

public:
	//constructor
	At1pCalibration(vector<Rate> cdsPremiums1,vector<Real> discount_rate1, Date valueDate1, vector<Period> survTenor1):
		valueDate(valueDate1),cdsPremiums(cdsPremiums1),discount_rate(discount_rate1),survTenor(survTenor1)
	{
		calendar = TARGET();
	}
	
	//get the vector for the tenor of survival probabilities
	vector<Period> get_survTenor();

	//get the vector of the discount factors
	vector<Real> get_discount_rate();

	//get the vector of CDS premiums
	vector<Rate> get_cdsPremiums();

	//get the vector of calibrated vol
	vector<Real> get_cal_vol();

	//reset the vector of the tenor of survival probabilities
	void set_survTenor(vector<Period> new_survTenor);

	//reset the vector of the discount factors
	void set_discount_rate(vector<Real> new_discount_rate);

	//reset the CDS premium input
	void set_cdsPremiums(vector<Rate> new_cdsPremiums);

	//reset the value date
	void set_valueDate(Date new_date);

	//The implementation of the AT1P model to calibrate the CDS spreads
	vector<Real> AT1P_CDS_calibration();

	//bisection method
	Date bisection(Date startrange,Date endrange,Real parameters_input);

	//generated the default time
	Date Get_default_time_AT1P();

};


