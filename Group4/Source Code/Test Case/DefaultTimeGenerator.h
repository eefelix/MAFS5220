#pragma warning(disable:4819)
#include <ql/quantlib.hpp>
#include <fstream>
#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>

using namespace std;
using namespace QuantLib;

class DefaultTimeGenerator {

/*  class "DefaultTimeGenerator" is used for generating default time of single individuals 
	given hazard rate term structure inputted by users.  */


private:
	//vector<Date> oDates_;
	vector<Real> oSurvProb_;
	vector<Real> oHazardRates_;
	Date valueDate_;
	Calendar calendar;

	// ## Uniform random variable generator using seed
	Real UniformDistributionGenerator(double range_from, double range_to);

public:
	vector<Date> oDates_;
	// ## Constructor #1
	DefaultTimeGenerator(Date valueDate);

	//## Constructor #2
	//The constructor has been changed for convenience on implementing the class "BVA_calculator" 
	DefaultTimeGenerator(Date valueDate,vector<Date> oDates, vector<Real> oHazardRates)
	{
		valueDate_ = valueDate;
		calendar = TARGET();
		oDates_ = oDates;
		oHazardRates_ = oHazardRates;
	}

	//// ## to be develop...
	//void SetSurvivalProbability(vector<Date> oDates, vector<Real> oSurvProb);


	// ## users need to input a pair of vectors : dates/hazard rates
	void SetHazardRates(vector<Date> oDates, vector<Real> oHazardRates);

	// ## function to generate default time of single individual party
	Date GetDefaultTime();


};
