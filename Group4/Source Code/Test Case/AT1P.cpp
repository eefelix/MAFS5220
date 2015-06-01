#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <fstream>
#include <map>
#include "AT1P.h"



using namespace std;
using namespace QuantLib;


#pragma region VolatilityTS_
	//Code in Region
	VolatilityTS::VolatilityTS() {
		dc = Actual365Fixed();
		oSize = 0;
	}

	VolatilityTS::VolatilityTS(Date valueDate, std::vector<Date> oDates): valueDate_(valueDate), oDates_(oDates) {
		dc = Actual365Fixed();
		oSize = oDates_.size();
		oSchedule_ = Schedule(oDates_);
	}
	

	void VolatilityTS::UpdateStepwiseVolatility(const Array& x, int iTenor) {
		tsVolatility.insert(std::pair<Date, Real>(oDates_[iTenor], x[0]));    // insert the pair<volatility tenor, vol> into the term structure vector
	}

	//calculate the sum of square of the calibrated volatility for the AT1P's usage
	Real VolatilityTS::GetSumOfVolSquare(Date endDate) {
		Real SumOfVolSqr_ = 0.0;

		for (i = 0; i < oSize; ++i)
		{
			if (endDate <= oDates_[i]) {
				Date LboundDate = (i == 0) ? valueDate_ : oDates_[i - 1];
				SumOfVolSqr_ += pow(tsVolatility[oDates_[i]], 2) * dc.yearFraction(LboundDate, endDate);
				return SumOfVolSqr_;
			}
			else if (endDate > oDates_[i]) {
				if (oDates_[i] != oSchedule_.endDate()) {
					Date LboundDate = (i == 0) ? valueDate_ : oDates_[i - 1];
					SumOfVolSqr_ += pow(tsVolatility[oDates_[i]], 2) * dc.yearFraction(LboundDate, oDates_[i]);
				}
				else if (oDates_[i] == oSchedule_.endDate()) {
					Date LboundDate = (i == 0) ? valueDate_ : oDates_[i - 1];
					SumOfVolSqr_ += pow(tsVolatility[oDates_[i]], 2) * dc.yearFraction(LboundDate, endDate);
					return SumOfVolSqr_;
				}
			}
		}
	}


#pragma endregion

//Calculate the NPV of the credit default swap
#pragma region CreditDefaultSwapV2_

	//constructor
	CreditDefaultSwapV2::CreditDefaultSwapV2(Protection::Side side,
		Real notional,
		const Date& protectionStart)
		: side_(side), notional_(notional), protectionStart_(protectionStart){
			calendar = TARGET();
			oAT1P_ = AT1P();
			lossGivenDefault = 0.6;
			risk_free = 0;
	}
	//constructor#2
	CreditDefaultSwapV2::CreditDefaultSwapV2(Protection::Side side,
		Real notional,
		Rate spread,
		Schedule schedule,
		const Date& protectionStart)
		: side_(side), notional_(notional), spread_(spread), schedule_(schedule), protectionStart_(protectionStart){
			calendar = TARGET();
			oAT1P_ = AT1P();
			lossGivenDefault = 0.6;
			risk_free = 0;
	}

	//reset the CDS premium
	void CreditDefaultSwapV2::SetPremium(Rate premium) {
		spread_ = premium;
	}
	//reset the schedule
	void CreditDefaultSwapV2::SetSchedule(Schedule schedule) {
		schedule_ = schedule;
	}
	//set the survivual probabilty for calculating the CDS NPV
	void CreditDefaultSwapV2::SetAT1PModel(AT1P oAT1P) {
		oAT1P_ = oAT1P;
	}

	//update the risk free rate during the calculation of the CDS NPV
	void CreditDefaultSwapV2::SetDis_fac(Real risk_freerate)
	{
		risk_free = risk_freerate;
	}
	
	//for the discrete calculation of the discount factor
	Real CreditDefaultSwapV2::pow(Real x, int power)
	{
		Real y = x;
		for(int i = 1; i< power; i++)
			x = x*y;

		return x;
	}

	//calculation of the discount factor, the risk free rate can be assumed to be the US treasury rate
	//http://www.treasury.gov/resource-center/data-chart-center/interest-rates/Pages/TextView.aspx?data=yield
	//After the advise from tutor Felix, the risk free rate is set to be flat at 5%
	Real CreditDefaultSwapV2::discount_factor(Date discount_date)
	{
		//continuous compounding
		//return exp(-5.3597/100*(discount_date-protectionStart_)/365);
		return exp(-risk_free/100*(discount_date-protectionStart_)/365);
	}
	
	//Calculate the PV of default leg 
		Real CreditDefaultSwapV2::DefaultLegPV(){
		
		Real defaultLegPV_ = 0.0;
		vector<Date> oCdsDates = schedule_.dates();
		int i;

		for (i = 1; i < schedule_.size() ; ++i)
		{
			if (i == 1) {
				defaultLegPV_ += notional_ * lossGivenDefault * (1 - oAT1P_.GetSurvivalProbability(oCdsDates[i])) * discount_factor(oCdsDates[i]);
			}
			else {
				defaultLegPV_ += notional_ * lossGivenDefault * (oAT1P_.GetSurvivalProbability(oCdsDates[i - 1]) - oAT1P_.GetSurvivalProbability(oCdsDates[i]))* discount_factor(oCdsDates[i]);
			}
		}

		return defaultLegPV_;
	}

		//Calculate the PV of coupon leg 
		Real CreditDefaultSwapV2::CouponLegPV(){
			
		Real couponLegPV_ = 0.0;
		vector<Date> oCdsDates = schedule_.dates();
		int j, k;
		Date startDate, endDate;

		for (j = 1; j < schedule_.size(); ++j){
			if (j == 1) {

				// for expected payoff of CDS coupon leg at CDS premium payment date
				couponLegPV_ += notional_ * (oCdsDates[j] - protectionStart_) / 365 * spread_
					* oAT1P_.GetSurvivalProbability(oCdsDates[j]) * discount_factor(oCdsDates[j]);
				// for expected payoff of CDS coupon legs in case underlying default between CDS premium payment dates
				k = 0;
				while (calendar.adjust(protectionStart_ + (k+1) * Months, Following) < oCdsDates[j]) 
				{	
					startDate = k == 0 ? protectionStart_ : calendar.adjust(protectionStart_ + k * Months, Following);
					
					endDate = calendar.adjust(protectionStart_ + (k + 1) * Months, Following);
					
					couponLegPV_ += notional_ * (endDate - protectionStart_) / 365 * spread_
						* (oAT1P_.GetSurvivalProbability(startDate) - oAT1P_.GetSurvivalProbability(endDate)) * discount_factor(oCdsDates[j]);
					//couponLegPV_1 += notional_ * (oCdsDates[j] - protectionStart_) / 365 * spread_ * (1-oAT1P_.GetSurvivalProbability(oCdsDates[j])) * discount_factor(oCdsDates[j]);
					k++;
				}
			}
			else {
				// for expected payoff of CDS coupon leg at CDS premium payment date
				
				couponLegPV_ += notional_ * (oCdsDates[j] - oCdsDates[j - 1]) / 365 * spread_
					* oAT1P_.GetSurvivalProbability(oCdsDates[j]) * discount_factor(oCdsDates[j]);
				
				// for expected payoff of CDS coupon legs in case underlying default between CDS premium payment dates
				k = 0;
	
				while (calendar.adjust(oCdsDates[j - 1] + (k+1) * Months, Following) < oCdsDates[j]) 
				{	
					startDate = k == 0 ? oCdsDates[j - 1] : calendar.adjust(oCdsDates[j - 1] + k * Months, Following);
					
					endDate = calendar.adjust(oCdsDates[j - 1] + (k+1) * Months, Following);
					
					couponLegPV_ += notional_ * (endDate - oCdsDates[j - 1]) / 365 * spread_ * (oAT1P_.GetSurvivalProbability(startDate) - oAT1P_.GetSurvivalProbability(endDate)) * discount_factor(oCdsDates[j]);
					
					k++;
				} 
				//couponLegPV_1 += notional_ * (oCdsDates[j] - protectionStart_) / 365 * spread_ * (1-oAT1P_.GetSurvivalProbability(oCdsDates[j])) * discount_factor(oCdsDates[j]);
			}
		}
		return couponLegPV_;
	}
	
	//calculate the NPV by netting the PV of default leg and coupon leg
	Real CreditDefaultSwapV2::NPV() {
		if (side_ == Protection::Buyer) {
			return DefaultLegPV() - CouponLegPV();
		}
		else {
			return -DefaultLegPV() + CouponLegPV();
		}

	}

#pragma endregion


#pragma region AT1P_

	AT1P::AT1P() {}

	AT1P::AT1P(Real V0, Real H0, Real B0, vector<Date> survDates, Date valueDate) : V0_(V0), H0_(H0), B0_(B0), survDates_(survDates), valueDate_(valueDate){
		sigma_ = VolatilityTS();
	}

	
	Real AT1P::GetSurvivalProbability(Date endDate) {
		// calculate the survivual probability by different tenors with calibrated volatilities
		//Analytical formula please refer to Brigo's book
		Real oQT;
		Real sumOfVolSqr = sigma_.GetSumOfVolSquare(endDate);
		d1 = (log(V0_ / H0_) + (2 * B0_ - 1) / 2 * sumOfVolSqr) / sqrt(sumOfVolSqr);
		d2 = (log(H0_ / V0_) + (2 * B0_ - 1) / 2 * sumOfVolSqr) / sqrt(sumOfVolSqr);
		boost::math::normal_distribution<>stdND(0, 1);
		oQT = cdf(stdND, d1) - pow((H0_ / V0_), (2 * B0_ - 1))*cdf(stdND, d2);
		return oQT;
	}

	//update the volatility term structure for calculating the survivual probability
	void AT1P::UpdateVolatilityTS(VolatilityTS oVolatilityTS) {
		sigma_ = oVolatilityTS;
	}

	VolatilityTS AT1P::get_sigma()
	{
		return sigma_;
	}
#pragma endregion


#pragma region FunctionCdsPv_

	FunctionCdsPv::FunctionCdsPv(VolatilityTS oVolTS, CreditDefaultSwapV2 oCDS, AT1P oAT1P) :
		oVolTS_1(oVolTS), oCDS_1(oCDS), oAT1P_1(oAT1P) {
			iTenor_ = 0;
	}

	FunctionCdsPv::FunctionCdsPv(VolatilityTS oVolTS, CreditDefaultSwapV2 oCDS, AT1P oAT1P, vector<Period> survTenor, vector<Rate> cdsPremiums, Date valueDate) :
		oVolTS_1(oVolTS), oCDS_1(oCDS), oAT1P_1(oAT1P), survTenor_(survTenor), cdsPremiums_(cdsPremiums), valueDate_(valueDate){
			iTenor_ = 0;
	}

	//reset the calibration tenor
	void FunctionCdsPv::SetTenorToCalibrate(int iTenor) {
		iTenor_ = iTenor;
	}
	
	//update the volatility term structure for calculating the CDS NPV and optimization
	void FunctionCdsPv::UpdateVolatilityTS(VolatilityTS oVolTS) {
		oVolTS_1 = oVolTS;
	}

	
	// Override the cost function 
	Real FunctionCdsPv::value(const Array& x) const{
		Real res = 0;
		Calendar calendar = TARGET();
		Schedule cdsSchedule;
		//setup the dummy classes for vol ts, at1p and cds npv calculation 
		VolatilityTS oVolTS_2 = oVolTS_1;
		CreditDefaultSwapV2 oCDS_2 = oCDS_1;
		AT1P oAT1P_2 = oAT1P_1;

		oVolTS_2.UpdateStepwiseVolatility(x, iTenor_);
		oAT1P_2.UpdateVolatilityTS(oVolTS_2);
		oCDS_2.SetAT1PModel(oAT1P_2);

	
		oCDS_2.SetPremium(cdsPremiums_[iTenor_]);
		cdsSchedule = MakeSchedule().from(valueDate_)
			.to(calendar.adjust(valueDate_ + survTenor_[iTenor_], Following))
			.withFrequency(Quarterly)
			.withCalendar(calendar)
			.withTerminationDateConvention(Unadjusted)
			.withRule(DateGeneration::Forward);
		oCDS_2.SetSchedule(cdsSchedule);

		//set the cost function as the NPV function of the CDS for optimitzation 
		//in order to the find the root (calibrated vol) for the NPV function
		res = abs(oCDS_2.NPV()); 
		return res;
	}

	Disposable<Array> FunctionCdsPv::values(const Array& x) const{
		Array res(1);
		res[0] = value(x);
		return res;
	}

#pragma endregion

#pragma region At1pCalibration_

//get the vector for the tenor of survival probabilities
vector<Period> At1pCalibration::get_survTenor()
{
	//print the results for testing purpose
	/*
	std::vector<Period>::iterator it;
	cout << "The survival probabilities generated from AT1P CDS calibration: " << endl;
	for (it = survTenor.begin(); it != survTenor.end(); ++it)
			cout << "survival tenor: " << *(it) << endl;*/
	return survTenor;
}

//get the vector of the discount factors
vector<Real> At1pCalibration::get_discount_rate()
{
	//print the results for testing purpose
	/*
	std::vector<Real>::iterator it;
	
	cout << "The discount factors for the CDS NPV: " << endl;
	for (it = discount_rate.begin(); it != discount_rate.end(); ++it)
			cout << "Discount rate: " << *(it) << endl;*/

	return discount_rate;
}

//get the vector of CDS premiums
vector<Rate> At1pCalibration::get_cdsPremiums()
{
	/*
	//print the results for testing purpose
	std::vector<Rate>::iterator it;
	cout << "The survival probabilities generated from AT1P CDS calibration: " << endl;
	for (it = cdsPremiums.begin(); it != cdsPremiums.end(); ++it)
			cout << "cds Premiums: " << *(it) << endl;*/

	return cdsPremiums;
}

//get the vector of calibrated vol from AT1P
vector<Real> At1pCalibration::get_cal_vol()
{
	return calibrated_vol;
}

//reset the vector of the tenor of survival probabilities
//vector<Period> At1pCalibration::set_survTenor(vector<Period> new_survTenor)
void At1pCalibration::set_survTenor(vector<Period> new_survTenor)
{
	survTenor.clear();
	for (Size i = 0; i < new_survTenor.size(); i++) {
		survTenor.push_back(new_survTenor[i]);
	}
	//return survTenor;
}

//reset the vector of the discount factors
//vector<Real> At1pCalibration::set_discount_rate(vector<Real> new_discount_rate)
void At1pCalibration::set_discount_rate(vector<Real> new_discount_rate)
{
	discount_rate.clear();
	for (Size i = 0; i < new_discount_rate.size(); i++) {
		discount_rate.push_back(new_discount_rate[i]);
	}
	//return discount_rate;
}

//reset the CDS premium input
//vector<Rate> At1pCalibration::set_cdsPremiums(vector<Rate> new_cdsPremiums)
void At1pCalibration::set_cdsPremiums(vector<Rate> new_cdsPremiums)
{
	cdsPremiums.clear();
	for (Size i = 0; i < new_cdsPremiums.size(); i++) {
		cdsPremiums.push_back(new_cdsPremiums[i]);
	}
	//return cdsPremiums;
}

//reset the value date
//Date At1pCalibration::set_valueDate(Date new_date)
void At1pCalibration::set_valueDate(Date new_date)
{
	valueDate = new_date;
}

//The implementation of the AT1P model to calibrate the CDS spreads into the volatilities term structure
vector<Real> At1pCalibration::AT1P_CDS_calibration()
{
	//clear the previous memories before a new calibration
	survDate.clear();
	at1p_prob.clear();
	calibrated_vol.clear();

	//calculate the actual date of the tenors
	for (Size i = 0; i < survTenor.size(); i++) {
		survDate.push_back(calendar.adjust(valueDate + survTenor[i], Following));
	}
	
	//invoke the volatility term structures
	oVolTS = VolatilityTS(valueDate, survDate);

	//invoke the class for CDS NPV calculation
	CreditDefaultSwapV2 oCDS = CreditDefaultSwapV2(Protection::Buyer, 100, valueDate);
	
	//invoke the class of AT1P model with the following setting: REC=40%, B=0, H=0.4
	oAT1P = AT1P(1, 0.4, 0, survDate, valueDate);

	//setup the cost function for optimization
	FunctionCdsPv oFunctionCdsPv = FunctionCdsPv(oVolTS, oCDS, oAT1P, survTenor, cdsPremiums, valueDate);

	Real solved_vol;
	
	//setup the constriants for performing the simplex optimaiztion 
		PositiveConstraint constraint_evt; //NoConstraint constraint_evt;
		EndCriteria end_crit_evt(10000, 1000, 1e-12, 1e-12, 1e-12);
		Simplex solver_evt(0.00001);
		EndCriteria::Type solved_evt;
		Problem problem_evt(oFunctionCdsPv, constraint_evt, Array(1, 0.2));

		// Run the optimization
	for (int k = 0; k < survDate.size(); k++) {
		oFunctionCdsPv.SetTenorToCalibrate(k);
		oCDS.SetDis_fac(discount_rate[k]);
		
		solved_evt = solver_evt.minimize(problem_evt, end_crit_evt);

		solved_vol = problem_evt.currentValue().at(0);

		calibrated_vol.push_back(solved_vol);
/*		
		//print the results
		cout << "Tenor: " << survDate[k] << ", CDS Spreads " << cdsPremiums[k] << ", Calibrated Vol " << solved_vol << endl;
		cout << endl;
*/		
		oVolTS.UpdateStepwiseVolatility(Array(1, solved_vol), k);  //calibrated vol is updated
		oFunctionCdsPv.UpdateVolatilityTS(oVolTS);  
		
		oAT1P.UpdateVolatilityTS(oVolTS);
		//get the survivual probability from the calibrated volatility
		at1p_prob.push_back(oAT1P.GetSurvivalProbability(survDate[k])); 
		
	}
/*
	std::vector<Real>::iterator it;
	//print the results
	cout << "The survival probabilities generated from AT1P CDS calibration: " << endl;
	for (it = at1p_prob.begin(); it != at1p_prob.end(); ++it)
			cout << "survival prob: " << *(it) << endl;
*/			
	return at1p_prob;
}

// ## Uniform random variable generator using seed
Real At1pCalibration::UniformDistributionGenerator(double range_from, double range_to) 
{
	Real uniform_r_v = 0.0;
	unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	static MersenneTwisterUniformRng unifMt(seed);
	uniform_r_v = unifMt.next().value;

	return uniform_r_v;
}

//Generate the default time for the counterparty by AT1P model and bisection method
Date At1pCalibration::Get_default_time_AT1P()
{
	vector<Date> survDate_complete;
	vector<Real> survProb_complete;

	//insert the first element to the default probability as well as the corresponding date
	survDate_complete.push_back(valueDate);
	survProb_complete.push_back(1);

	//setup the new vector for the survival probabilities and its corresponding vectors
	for (Size i = 0; i < survTenor.size(); i++) {
		survDate_complete.push_back(calendar.adjust(valueDate + survTenor[i], Following));
		survProb_complete.push_back(at1p_prob[i]);
	}

	Date valueDate_default(31, Dec, 2100);
	Real surv_Prob = UniformDistributionGenerator(0.0, 1.0); //random number simulated by uniform distribution
//	cout << surv_Prob << endl;

	int counter=0;

	//Check if the survival probabilities generated by AT1P for all tenor are all to the simulated random numbers
	if(survProb_complete[survProb_complete.size()-1] < surv_Prob)
	{
		for (int i = 0; i < survProb_complete.size(); i++)
		{   //When the survival probability generated by the AT1P model is first becoming negative 
			//(i.e. the root of (AT1P's survival prob - random number) is within this tenor
			if (survProb_complete[i] <= surv_Prob)
			{	
				counter = i;
				break;
			}
		}
		//use bisection method to calculate the root (i.e. default time) so that the AT1P survival prob is equal to the random number 
		return calendar.adjust(bisection(survDate_complete[counter-1],survDate_complete[counter],surv_Prob), Following);
	}
	else
		return valueDate_default;
}

//bisection method
Date At1pCalibration::bisection(Date startrange,Date endrange,Real parameters_input)
{
	Real sur_A, sur_B;
	Date a, b, mid;
	Real epsilon = 0.000001;		//Set the epsilon for the error tolerance for the bisection method
	
	mid = startrange + ((endrange - startrange)/2); 
	a = startrange; 
	b = endrange;
	sur_A = oAT1P.GetSurvivalProbability(startrange)-parameters_input; //setup initial end points for the interval of the root
	sur_B = oAT1P.GetSurvivalProbability(b)-parameters_input;

	while (abs(mid - a) > epsilon && abs(sur_A) > epsilon)
	{ //Bisection method: Test if the points input are close to the root
		sur_A = oAT1P.GetSurvivalProbability(a)-parameters_input; 
		sur_B = oAT1P.GetSurvivalProbability(mid)-parameters_input;	

		//if the mid point input to the function is on the same side of a, replace a, else replace b 
		if ((sur_A*sur_B) > 0)
			a = mid;
	
		if ((sur_A*sur_B) < 0)
			b = mid;

		//Make a new mid point with the new set of a and b
		mid = a + ((b - a) / 2);
	}
		return a;
}
#pragma endregion

	


	
	

	


