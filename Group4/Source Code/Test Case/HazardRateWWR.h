#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <boost/timer.hpp>
#include "MCgenerator.h"

using namespace std;
using namespace QuantLib;

// enum for Hazard Rate function
enum WWRHazardModel{ EXPONENTIALLY, LINEARLY };

// ====================================================================================================
// Class: HazardRateWWRCalibrate
// Description: Model Calibration class for calibrating a(ti) in WWR hazard rate model as proposed 
//				in academic paper <CVA AND WRONG WAY RISK> by John Hull and Alan White
// ====================================================================================================
class HazardRateWWRCalibrate : public CostFunction{

private:
	
	Date valueDate_;
	vector<Date> cds_maturities_;
	vector<Real> cds_spreads_;
	Real recoveryRate_;
	int no_of_paths;
	int iTenor_;
	WWRHazardModel iWWRHazardModel_;
	vector<map <Date, Real> > simulatedPrices_;

	vector<Real> h_t_;
	vector<Real> a_t_;		// term structure of a(ti) in WWR hazard rate model
	Real b_;				// fixed input b in in WWR hazard rate model

public:
	//Constructor #1
	HazardRateWWRCalibrate(Date valueDate);

	// Need to be input before doing calibration
	void SetSensitivityB(Real b);

	// Need to be input before doing calibration
	void SetCdsPayDate(vector<Date> cds_maturities);

	// Need to be input before doing calibration
	void SetCdsSpread(vector<Real> cds_spreads);

	// Need to be input before doing calibration
	void SetRecoveryRate(Real recoveryRate);

	// Need to be input before doing calibration
	void SetSimulatedPrices(vector<map <Date, Real> > simulatedPrices);

	// Need to be input before doing calibration for each tenor
	void SetTenorToCalibrate(int iTenor);

	// Need to be input before doing calibration for each tenor
	void UpdateCalibratedParameter(vector<Real> a_t);

	// Need to be set in order to calibrate either exponential or linear hazard rate model
	void SetWWRHazardModel(WWRHazardModel iWWRHazardModel);

	// Override the cost function
	Real value(const Array& x) const;

	Disposable<Array> values(const Array& x) const;
};

// ====================================================================================================
// Class: HazardRateWWR
// Description: Execute the calibration and implementation of WWR Hazard Rate model as proposed 
//				in academic paper <CVA AND WRONG WAY RISK> by John Hull and Alan White
// ====================================================================================================
class HazardRateWWR{

private:
	Calendar calendar;
	Date valueDate_;
	vector<Real> cds_spreads_;
	vector<Period> cds_tenors_;
	Real recoveryRate_;
	int noOfPaths;
	WWRHazardModel iWWRHazardModel_;
	vector<map <Date, Real> > simulatedPrices;
	Position::Type forwardType_;

	Real UniformDistributionGenerator(double range_from, double range_to);

public:

	vector<Date> cds_maturities_;
	vector<Real> a_t_;
	Real b_;

	HazardRateWWR(Date valueDate);

	void SetSensitivityB(Real b);

	// Need to be input before doing calibration/implementation
	void SetCdsPayDate(vector<Date> cds_maturities);

	// Need to be input before doing calibration/implementation
	void SetCdsTenor(vector<Period> cds_tenors);

	// Need to be input before doing calibration/implementation
	void SetCdsSpread(vector<Real> cds_spreads);

	// Need to be input before doing calibration/implementation
	void SetRecoveryRate(Real recoveryRate);

	// Need to be input before doing calibration/implementation
	void SimulateStockPrices(MCgenerator oMCGenerator);

	// Need to be set in order to calibrate either exponential or linear hazard rate model
	void SetWWRHazardModel(WWRHazardModel iWWRHazardModel);

	// will call class HazardRateWWRCalibrate to perform calibration of a(ti) in the model
	void PerformCalibration();

	// for implementing the calibrated hazard rate model
	// calculate survival probability up to tenor i
	// input parameter includes payOffPath as hazard rate model is path dependent
	Real GetSurvivalProbability(int iTenor, map<Date, Real> *payOffPath);

	void SetForwardType(Position::Type forwardType);

};