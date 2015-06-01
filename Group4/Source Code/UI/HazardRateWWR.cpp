#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include "HazardRateWWR.h"

using namespace std;
using namespace QuantLib;


#pragma region HazardRateWWRCalibrate_
// ====================================================================================================
// Class: HazardRateWWRCalibrate
// Description: Model Calibration class for calibrating a(ti) in WWR hazard rate model as proposed 
//				in academic paper <CVA AND WRONG WAY RISK> by John Hull and Alan White
// ====================================================================================================

HazardRateWWRCalibrate::HazardRateWWRCalibrate(Date valueDate) {
	 no_of_paths = 0;
	  iTenor_ = 0;
	 valueDate_ = valueDate;
}

void HazardRateWWRCalibrate::SetSensitivityB(Real b) {
	b_ = b;
}

void HazardRateWWRCalibrate::SetCdsPayDate(vector<Date> cds_maturities) {
	cds_maturities_ = cds_maturities;
}

void HazardRateWWRCalibrate::SetCdsSpread(vector<Real> cds_spreads) {
	cds_spreads_ = cds_spreads;
}

void HazardRateWWRCalibrate::SetRecoveryRate(Real recoveryRate){
	recoveryRate_ = recoveryRate;
}

void HazardRateWWRCalibrate::SetSimulatedPrices(vector<map <Date, Real> > simulatedPrices) {
	simulatedPrices_ = simulatedPrices;
	no_of_paths = simulatedPrices.size();
}

void HazardRateWWRCalibrate::SetTenorToCalibrate(int iTenor) {
	iTenor_ = iTenor;
}

void HazardRateWWRCalibrate::UpdateCalibratedParameter(vector<Real> a_t) {
	a_t_ = a_t;
}

void HazardRateWWRCalibrate::SetWWRHazardModel(WWRHazardModel iWWRHazardModel){
	iWWRHazardModel_ = iWWRHazardModel;
}

// Override the cost function 
Real HazardRateWWRCalibrate::value(const Array& x) const{
	Real res = 0.0;

	Date startDate;
	Date endDate;
	Real instantHazardRate;
	Real sumHazardRate;
	Real sumSurvProbHazardRate = 0.0;
	Real avgSurvProbHazardRate = 0.0;
	Real survProbCdsTS = 0.0;
	Real a_t;		// dummy a(ti)
	Real w_t_;		// EAD at time t

	for (int j = 0; j < no_of_paths; j++) {

		const map <Date, Real> & pricePathJ = simulatedPrices_[j];
		sumHazardRate = 0.0;

		for (int i = 0; i <= iTenor_; i++) {
			startDate = (i == 0) ? valueDate_ : cds_maturities_[i - 1];
			endDate = cds_maturities_[i];
			w_t_ = pricePathJ.find(cds_maturities_[i])->second;
			a_t = (i == iTenor_) ? x[0] : a_t_[i];

			switch (iWWRHazardModel_)
			{
				case EXPONENTIALLY:
					instantHazardRate = exp(a_t + b_ * w_t_);
					break;
				case LINEARLY:
					instantHazardRate = log(1 + exp(a_t + b_ * w_t_));
					break;
			}
			sumHazardRate += -instantHazardRate * Actual365Fixed().yearFraction(startDate, endDate);
		}

		sumSurvProbHazardRate += exp(sumHazardRate);		// sum the survival prob of path i calculated by model
	}

	// Expected survival probability calculated by HazardRate model
	avgSurvProbHazardRate = sumSurvProbHazardRate / no_of_paths;
	// Actual survival probability implied by stepwise CDS term structure
	survProbCdsTS = exp(-cds_spreads_[iTenor_] * Actual365Fixed().yearFraction(valueDate_, cds_maturities_[iTenor_]) / (1 - recoveryRate_));
	// Calibrate by minimizing the diff between model value and actual value 
	res = abs(avgSurvProbHazardRate - survProbCdsTS);   //<===== ###########

	return res;
}

Disposable<Array> HazardRateWWRCalibrate::values(const Array& x) const{
	Array res(1);
	res[0] = value(x);
	return res;
}

#pragma endregion


#pragma region HazardRateCalibrate_
// ====================================================================================================
// Class: HazardRateWWR
// Description: Execute the calibration and implementation of WWR Hazard Rate model as proposed 
//				in academic paper <CVA AND WRONG WAY RISK> by John Hull and Alan White
// ====================================================================================================

// Constructor #1
HazardRateWWR::HazardRateWWR(Date valueDate){
	calendar = TARGET();
	valueDate_ = valueDate;
	noOfPaths = 1000;	// number of MCGenerated paths used for calibration of a(ti)
}

void HazardRateWWR::SetSensitivityB(Real b){
	b_ = b;
}

void HazardRateWWR::SetCdsPayDate(vector<Date> cds_maturities) {
	cds_maturities_ = cds_maturities;
}

// Need to be input before doing PerformCalculation
void HazardRateWWR::SetCdsTenor(vector<Period> cds_tenors) {
	cds_tenors_ = cds_tenors;
	cds_maturities_.clear();
	for (Size i = 0; i < cds_tenors_.size(); i++) {
		cds_maturities_.push_back(calendar.adjust(valueDate_ + cds_tenors_[i], Following));
	}
}


void HazardRateWWR::SetCdsSpread(vector<Real> cds_spreads) {
	cds_spreads_ = cds_spreads;
}

void HazardRateWWR::SetRecoveryRate(Real recoveryRate){
	recoveryRate_ = recoveryRate;
}

void HazardRateWWR::SimulateStockPrices(MCgenerator oMCGenerator) {
	simulatedPrices.clear();
	oMCGenerator.SetForwardType(forwardType_);		// multiply payoff by -1.0 for short forward

	for (int i = 0; i < noOfPaths; i++) {
		simulatedPrices.push_back(oMCGenerator.SimulatePayoffPath_V1(cds_maturities_));
	}
}

void HazardRateWWR::SetWWRHazardModel(WWRHazardModel iWWRHazardModel){
	iWWRHazardModel_ = iWWRHazardModel;
}

void HazardRateWWR::PerformCalibration() {

	HazardRateWWRCalibrate oHazardRateWWRCalibrate = HazardRateWWRCalibrate(valueDate_);
	oHazardRateWWRCalibrate.SetSensitivityB(b_);
	oHazardRateWWRCalibrate.SetCdsPayDate(cds_maturities_);
	oHazardRateWWRCalibrate.SetCdsSpread(cds_spreads_);
	oHazardRateWWRCalibrate.SetRecoveryRate(recoveryRate_);
	oHazardRateWWRCalibrate.SetWWRHazardModel(iWWRHazardModel_);
	oHazardRateWWRCalibrate.SetSimulatedPrices(simulatedPrices);

	Real solved_a_t = 0.0;
	a_t_.clear();

	//cout << "Calibration result: " ;

	for (int i = 0; i < cds_maturities_.size(); i++) {
		// calibrate tentor by tenor
		oHazardRateWWRCalibrate.SetTenorToCalibrate(i);

		NoConstraint constraint_evt;
		EndCriteria end_crit_evt(500, 100, 1e-5, 1e-5, 1e-5);
		Problem problem_evt(oHazardRateWWRCalibrate, constraint_evt, Array(1, -1.0));
		Simplex solver_evt(0.01);
		EndCriteria::Type solved_evt = solver_evt.minimize(problem_evt, end_crit_evt);
		solved_a_t = problem_evt.currentValue().at(0);
		//cout << "a(t=" << i << ") was " << solved_a_t << ". " << endl;
		a_t_.push_back(solved_a_t);		// save the calibrated value a(ti) 

		oHazardRateWWRCalibrate.UpdateCalibratedParameter(a_t_);

	}
	//cout << endl;
}

Real HazardRateWWR::UniformDistributionGenerator(double range_from, double range_to) {
	Real uniform_r_v = 0.0;
	unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	static MersenneTwisterUniformRng unifMt(seed);
	uniform_r_v = unifMt.next().value;
	return uniform_r_v;
}

Real HazardRateWWR::GetSurvivalProbability(int iTenor, map<Date, Real> *payOffPath) {

	Real instantHazardRate;
	Real timePeriod;
	Real SumOfHazardRates = 0.0;
	Date startDate;
	Date endDate;
	Real ILongShort;

	ILongShort = (forwardType_ == Position::Short) ? -1.0 : 1.0;		// multiply payoff by -1.0 for short forward

	for (int i = 0; i <= iTenor; i++) {

		startDate = (i == 0) ? valueDate_ : cds_maturities_[i - 1];
		endDate = cds_maturities_[i];

		switch (iWWRHazardModel_)
		{
			case EXPONENTIALLY:
				instantHazardRate = exp(a_t_[i] + b_ * ILongShort * payOffPath->find(endDate)->second);
				break;
			case LINEARLY:
				instantHazardRate = log(1 + exp(a_t_[i] + b_ * ILongShort * payOffPath->find(endDate)->second));
				break;
		}

		timePeriod = Actual365Fixed().yearFraction(startDate, endDate);
		SumOfHazardRates += instantHazardRate * timePeriod;
	}

	return exp(-SumOfHazardRates);
}

void HazardRateWWR::SetForwardType(Position::Type forwardType) {
	forwardType_ = forwardType;
}

#pragma endregion