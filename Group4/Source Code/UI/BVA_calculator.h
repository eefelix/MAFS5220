#include "MCgenerator.h"
#include "DefaultTimeGenerator.h"
#include "HazardRateCalibrate.h"
#include "AT1P.h"
#include "HazardRateWWR.h"
//#include "Collateral.h"
#include <string>

using namespace std;
using namespace QuantLib;

class BVA_calculator
{
public:
	//constructor of Hazard rate calibration without collateral
	BVA_calculator(MCgenerator *MCgen, HazardRateCalibrate *own, HazardRateCalibrate *counterparty, int date_path, int stock_path):
		gen(MCgen),myself(own),cpty(counterparty), no_of_Date_path(date_path), no_of_stock_path(stock_path)
	{}

	//constructor of AT1P without collateral
	BVA_calculator(MCgenerator *MCgen, At1pCalibration *own, At1pCalibration *counterparty, Real own_recover, Real counter_recover,
		int date_path, int stock_path):
		gen(MCgen),myself_AT1P(own),cpty_AT1P(counterparty), no_of_Date_path(date_path), no_of_stock_path(stock_path),own_recovery_rate(own_recover),counterparty_recovery_rate(counter_recover)
	{}

	//constructor of Hazard rate calibration with collateral
	BVA_calculator(MCgenerator *MCgen, HazardRateCalibrate *own, HazardRateCalibrate *counterparty, int date_path, int stock_path, Real coll_LGD_inv, Real coll_LGD_cpty):
		gen(MCgen),myself(own),cpty(counterparty), no_of_Date_path(date_path), no_of_stock_path(stock_path),
		own_collateral_LGD(coll_LGD_inv),cpty_collateral_LGD(coll_LGD_cpty)
	{}

	//constructor of AT1P withcollateral
	BVA_calculator(MCgenerator *MCgen, At1pCalibration *own, At1pCalibration *counterparty, Real own_recover, Real counter_recover,
		Real own_coll_LGD, Real counter_coll_LGD,int date_path, int stock_path):
		gen(MCgen),myself_AT1P(own),cpty_AT1P(counterparty), no_of_Date_path(date_path), no_of_stock_path(stock_path),
			own_recovery_rate(own_recover),	counterparty_recovery_rate(counter_recover),own_collateral_LGD(own_coll_LGD),cpty_collateral_LGD(counter_coll_LGD)
	{}

	//member functions
	//calculate the BVA (i.e. CVA+DVA) by using the Hazard rate calibration method
	vector<map<string,Real>> BVA_cal(Date maturity, char collateral_indicator, Real inv_col_LGD, Real cpty_col_LGD);
	

		//calculate the BVA (i.e. CVA+DVA) by using the AT1P CDS calibration method
	vector<map<string,Real>> BVA_cal_AT1P(Date maturity, char collateral_indicator, Real inv_col_LGD, Real cpty_col_LGD);
	
	Real maximun(Real left, Real right); //return the maximum if the two parameters

	Real minimun(Real left, Real right);//return the minimum if the two parameters

	//Calculate the collateral inclusive BVA 
	Real Coll_incl_BVA_cal(pair<Real,Real> forward_pv,char CVA_DVA_ind, Real own_LGD, Real cpt_LGD);

	//reset the collateral LGD ratio posted by investor
	void set_inv_coll_LGD(Real new_inv_coll_LGD);

	//reset the collateral LGD ratio posted by investor
	void set_cpty_coll_LGD(Real new_cpty_coll_LGD);

	Real own_collateral_LGD;  //Loss given default for the collateral account
	Real cpty_collateral_LGD; //same as above
private:
	Real CVA;
	Real DVA;
	Date investorDefaultDate;
	Date counterptyDefaultDate;
	vector<Date> maturities;
	MCgenerator *gen;				//pointers for the class of the BVA calculation input
	HazardRateCalibrate *myself;
	HazardRateCalibrate *cpty;
	int no_of_Date_path;
	int no_of_stock_path;
	At1pCalibration *myself_AT1P;
	At1pCalibration *cpty_AT1P;
	Real own_recovery_rate;
	Real counterparty_recovery_rate;
	//Real own_collateral_LGD;  //Loss given default for the collateral account
	//Real cpty_collateral_LGD; //same as above

};


// =====================================================================================
// Class: BVACalculatorWWR
// Description: used for calculate UDVA & UCVA
// =====================================================================================
class BVACalculatorWWR
{
public:
	// Constructor
	BVACalculatorWWR(MCgenerator *oMCgenerator, HazardRateWWR *investorHazardRate,
		HazardRateWWR *counterpartyHazardRate, int date_path, int stock_path) :
		oMCgenerator_(oMCgenerator), investorHazardRate_(investorHazardRate),
		counterpartyHazardRate_(counterpartyHazardRate), no_of_Date_path(date_path), no_of_stock_path(stock_path){
		calendar = TARGET();
	};

	// calculate the UCVA & UDVA
	map<string, Real> BVA_cal(Date maturity);

private:
	Real CVA;
	Real DVA;
	Date investorDefaultDate;
	Date counterptyDefaultDate;
	vector<Date> maturities;
	MCgenerator *oMCgenerator_;				//pointers for the class of the BVA calculation input
	HazardRateWWR *investorHazardRate_;
	HazardRateWWR *counterpartyHazardRate_;
	int no_of_Date_path;
	int no_of_stock_path;
	Calendar calendar;
};