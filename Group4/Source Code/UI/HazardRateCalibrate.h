#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <boost/timer.hpp>


using namespace std;
using namespace QuantLib;

class HazardRateCalibrate {

/*  class "HazardRateCalibrate" is used for calibrating stepwise hazard rate term structure 
	given cds spread inputted by users.  */

private:
	Date valueDate_;
	boost::timer timer;
	Calendar calendar;
	boost::shared_ptr<Quote> flatRate;
	Handle<YieldTermStructure> tsCurve;
	vector<Real> cds_spreads_;
	vector<Period> cds_tenors_;
	vector<Date> cds_maturities_;
	
public:
	//Constructor #1
	HazardRateCalibrate(Date valueDate);

	//Constructor #2
	HazardRateCalibrate(Date valueDate, Real rate,Real recoveryRate,vector<Real> cds_spreads,vector<Period> cds_tenors)
	{
		calendar = TARGET();

		flatRate = boost::shared_ptr<Quote>(new SimpleQuote(rate));
		valueDate_ = valueDate;
		tsCurve = Handle<YieldTermStructure>(boost::shared_ptr<FlatForward>(
											new FlatForward(valueDate_, Handle<Quote>(flatRate),
											Actual365Fixed())));
		recovery_rate = recoveryRate;
		cds_spreads_ = cds_spreads;
		cds_tenors_ = cds_tenors;
		cds_maturities_.clear();
		for (Size i = 0; i < cds_tenors_.size(); i++) 
		{
		cds_maturities_.push_back(calendar.adjust(valueDate_ + cds_tenors_[i], Following));
		}
	}
	//These functions are built for implementing the function in class "BVA_calculator"
	Date get_valueday();
	vector<Date> get_cds_mat();
	vector<Date> get_hazard_dates();
	Calendar get_calendar();

	Real recovery_rate;
	vector<Date> oHazardDates;
	vector<Real> oHazardRates;

	void SetYieldRate(Real rate);
	void SetRecoveryRate(Real recoveryRate);
	void SetCdsSpread(vector<Real> cds_spreads);
	void SetCdsTenor(vector<Period> cds_tenors);
	
	void PerformCalibration();
};


