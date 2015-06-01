#include "HazardRateCalibrate.h"


//Return the value date
Date HazardRateCalibrate::get_valueday()
{
	return valueDate_;
}

//Return the vector of the CDS spreads input
vector<Date> HazardRateCalibrate::get_cds_mat()
{
	return cds_maturities_;
}

vector<Date> HazardRateCalibrate::get_hazard_dates()
{
	return oHazardDates;
}

//Return the calendar type set in this class
Calendar HazardRateCalibrate::get_calendar()
{
	return calendar;
}


// Constructor #1
HazardRateCalibrate::HazardRateCalibrate(Date valueDate) : valueDate_(valueDate) {
	calendar = TARGET();
}

// Need to be input before doing PerformCalculation
void HazardRateCalibrate::SetYieldRate(Real rate) {
	flatRate = boost::shared_ptr<Quote>(new SimpleQuote(rate));
	tsCurve = Handle<YieldTermStructure>(boost::shared_ptr<FlatForward>(
		new FlatForward(valueDate_, Handle<Quote>(flatRate),
		Actual365Fixed())));
}

// Need to be input before doing PerformCalculation
void HazardRateCalibrate::SetRecoveryRate(Real recoveryRate) {
	recovery_rate = recoveryRate;
}

// Need to be input before doing PerformCalculation
void HazardRateCalibrate::SetCdsSpread(vector<Real> cds_spreads) {
	cds_spreads_ = cds_spreads;
}

// Need to be input before doing PerformCalculation
void HazardRateCalibrate::SetCdsTenor(vector<Period> cds_tenors) {
	cds_tenors_ = cds_tenors;
	cds_maturities_.clear();
	for (Size i = 0; i < cds_tenors_.size(); i++) {
		cds_maturities_.push_back(calendar.adjust(valueDate_ + cds_tenors_[i], Following));
	}
}

// Calibrating hazard rates using CDS spread input
void HazardRateCalibrate::PerformCalibration() {

	std::vector<boost::shared_ptr<DefaultProbabilityHelper> > instruments;
	for (Size i = 0; i < cds_maturities_.size(); i++) 
	{
		instruments.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
			new SpreadCdsHelper(
			Handle<Quote>(boost::shared_ptr<Quote>(
			new SimpleQuote(cds_spreads_[i]))),
			cds_tenors_[i],
			0,
			calendar,
			Quarterly,
			Following,
			DateGeneration::Zero,
			Actual365Fixed(),
			recovery_rate,
			tsCurve)));
	}

	// Bootstrap hazard rates
	boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >
		hazardRateStructure(
		new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(
		valueDate_,
		instruments,
		Actual365Fixed()));

	vector<pair<Date, Real> > hr_curve_data = hazardRateStructure->nodes();

	oHazardDates.clear();
	oHazardRates.clear();

	for (Size i = 0; i < hr_curve_data.size(); i++) 
	{
		oHazardDates.push_back(hr_curve_data[i].first);
		oHazardRates.push_back(hr_curve_data[i].second);
	}

}
