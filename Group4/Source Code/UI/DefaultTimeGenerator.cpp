#include "DefaultTimeGenerator.h"

using namespace std;
using namespace QuantLib;


	// ## Uniform random variable generator using seed
	Real DefaultTimeGenerator::UniformDistributionGenerator(double range_from, double range_to) {
		Real uniform_r_v = 0.0;
		//unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		//static std::mt19937 generator(seed);
		//std::uniform_real_distribution<double> distribution(range_from, range_to);
		//uniform_r_v = distribution(generator);
		unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		static MersenneTwisterUniformRng unifMt(seed);
		uniform_r_v = unifMt.next().value;
		return uniform_r_v;
	}

	
	// ## Constructor #1
	DefaultTimeGenerator::DefaultTimeGenerator(Date valueDate):valueDate_(valueDate){
	calendar = TARGET();

	}

	//// ## Constructor #2: to be develop
	//void SetSurvivalProbability(vector<Date> oDates, vector<Real> oSurvProb) {
	//	oDates_ = oDates;
	//	oSurvProb_ = oSurvProb;
	//}

	// ## users need to input a pair of vectors : dates/hazard rates
	void DefaultTimeGenerator::SetHazardRates(vector<Date> oDates, vector<Real> oHazardRates) {
		oDates_ = oDates;
		oHazardRates_ = oHazardRates;
	}
	
	/* //Original method for the default time generator
	// ## function to generate default time of single individual party
	Date DefaultTimeGenerator::GetDefaultTime() {
		Real survProb = UniformDistributionGenerator(0.0, 1.0);
		Real accumulatedArea = 0.0;
		Real area_ = 0.0;
		Date startDate;
		Date endDate;
		int daysToAdd;
		for (int i = 1; i < oDates_.size(); ++i) {    // loop until the end of hazard rate term structure
			if (i != oDates_.size() - 1) {            
				startDate = (i == 1) ? valueDate_ : oDates_[i - 1];
				endDate = oDates_[i];
				area_ = oHazardRates_[i] * Actual365Fixed().yearFraction(startDate, endDate);    // equivalent to integration of hazard rate function up to time t
				if (exp(-accumulatedArea - area_) > survProb) {    // calculate the exponential of hazard rate integration
					accumulatedArea += area_;
				}
				else {
					
					daysToAdd = (int)((-log(survProb) - accumulatedArea) / oHazardRates_[i]) * 365;    // accept possible loss of data: conversion from double to int
					return calendar.adjust(oDates_[i] + daysToAdd * Days, Following);
				}
			}
			else {
				if (oDates_[i].year() + ((-log(survProb) - accumulatedArea) / oHazardRates_[i]) >= 2100) {    // quantlib limit the max year. We further limit it to be less than 2100 
					return Date(31, QuantLib::Dec, 2100);
				}
				else {
					daysToAdd = (int)((-log(survProb) - accumulatedArea) / oHazardRates_[i]) * 365;
					return calendar.adjust(oDates_[i] + daysToAdd * Days, Following);
				}
				
			}
		}
	}
	*/
	
	//New default time generator
	Date DefaultTimeGenerator::GetDefaultTime() {
		Real survProb = UniformDistributionGenerator(0.0, 1.0);
		Real accumulatedArea = 0.0;
		Real area_ = 0.0;
		Date endDate(31, Dec, 2100);
		DayCounter dc = Actual365Fixed();
		Real dummy_prob;

		vector<Date> new_date;
		
		new_date.push_back(valueDate_);

		//setup a new tenor vector which contains the value date as the beginning of the period 
		for (int i = 0; i < oDates_.size(); i++)
			new_date.push_back(oDates_[i]);
		
	//	cout << "sur :"<< survProb << endl;

		//accumulated the default probability by integrating the hazard rate until it reaches the random number
		for (int i = 0; i < new_date.size()-1; i++)
		{
			area_ = accumulatedArea;
			accumulatedArea = accumulatedArea + (oHazardRates_[i]*dc.yearFraction(new_date[i],new_date[i+1]));
			dummy_prob = exp(-accumulatedArea);
			
			//When the probability integrated by the hazard rate is become less than the random number (i.e. reaches the time that two numbers are eaual)
			if(dummy_prob<=survProb)
			{
			//	cout << "dummy :"<< dummy_prob << endl;
				endDate = new_date[i] + ((-log(survProb) - area_)/oHazardRates_[i]*365);
				endDate = calendar.adjust(endDate, Following);
				if(endDate.year() >= 2100)   // if the default date calculated is greater than 2100, set it back to the year 2100
					endDate = Date(31, Dec, 2100);
				break;
			}
		}

		return endDate;
	}