#include "BVA_calculator.h"

//return the max of the two input parameters
Real BVA_calculator::maximun(Real left, Real right)
{
	if(left > right)
		return left;
	else
		return right;
}

//return the min of the two input parameters
Real BVA_calculator::minimun(Real left, Real right)
{
	if(left < right)
		return left;
	else
		return right;
}

//to calculate the BVA for the forward with the data of the counterparty and our own
vector<map<string,Real>> BVA_calculator::BVA_cal(Date maturity, char collateral_indicator, Real inv_col_LGD, Real cpty_col_LGD)
{
	Real CVA = 0.0;		
	Real DVA = 0.0;
	set_inv_coll_LGD(inv_col_LGD); //reset the investor's collateral LGD level
	set_cpty_coll_LGD(cpty_col_LGD); //reset the counterparty's collateral LGD level

	//initial the class for generating the default time of counterparty and our own based on the data in the class "HazardRateCalibrate"
	DefaultTimeGenerator oInvestorDefaultTime(myself->get_valueday(), myself->get_hazard_dates(), myself->oHazardRates);
	DefaultTimeGenerator oCounterptyDefaultTime(cpty->get_valueday(), cpty->get_hazard_dates(), cpty->oHazardRates);

	//get the Loss given default data for both counterparties
	Real Investor_LGD = 1-myself->recovery_rate;
	Real counterparty_LGD = 1-cpty->recovery_rate;

	// start the for loop to calculate expected EAD
	for (int iPath = 0; iPath < no_of_Date_path; ++iPath) 
	{
		//calculate a sample of default times for both counterparties
		investorDefaultDate = oInvestorDefaultTime.GetDefaultTime();
		counterptyDefaultDate = oCounterptyDefaultTime.GetDefaultTime();

		//First to default check
		// check if there's any default event before maturity 
		if (investorDefaultDate < maturity || counterptyDefaultDate < maturity)
		{	
			//check if both counterparties default before the term of the forward
			if (investorDefaultDate < maturity && counterptyDefaultDate < maturity)
			{
				if (investorDefaultDate < counterptyDefaultDate) 
				{   //if the collateral option is selected ("Y/y"),the montle carlo simulation will include the calculation of the collateral amount
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						DVA = DVA - Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, investorDefaultDate, 0),'d',Investor_LGD,counterparty_LGD);
					else
				//if our own default happens first, calculate the DVA (without collateral)
						DVA = DVA + (Investor_LGD*maximun(-(gen->MCgeneration(no_of_stock_path, investorDefaultDate, 0)),0));
				}
				else
				{//if the collateral option is selected ("Y/y"),the montle carlo simulation will include the calculation of the collateral amount
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						CVA = CVA + (Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, counterptyDefaultDate, 0),'c',Investor_LGD,counterparty_LGD));
					else
				//if counterparty default first, calculate the CVA (without collateral)
						CVA = CVA + (counterparty_LGD*maximun(gen->MCgeneration(no_of_stock_path, counterptyDefaultDate, 0),0));
				}
			}
			else
			{//check if there is only one of the counterparty default before the maturity of the forward
				if (investorDefaultDate < maturity)
				{
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						DVA = DVA - Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, investorDefaultDate, 0),'d',Investor_LGD,counterparty_LGD);
					else
						DVA = DVA + (Investor_LGD*maximun(-(gen->MCgeneration(no_of_stock_path, investorDefaultDate, 0)),0));
				}
				else
				{	
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						CVA = CVA + (Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, counterptyDefaultDate, 0),'c',Investor_LGD,counterparty_LGD));
					else
						CVA = CVA + (counterparty_LGD*maximun(gen->MCgeneration(no_of_stock_path, counterptyDefaultDate, 0),0));
				}
			}
		}
		// if no default event in path i, then skip to next iterator
		else 
		{
			continue;
		}
	}

	// Expected values for the sampled CVA and DVA
	DVA = DVA/ no_of_Date_path;
	CVA = CVA/ no_of_Date_path;

	//save the CVA and DVA into the output vector
	map<string,Real> dummy_var;
	vector<map<string,Real>> BVA;
	dummy_var.insert(pair<string,Real>("CVA",CVA));
	BVA.push_back(dummy_var);
	map<string,Real> dummy_var1;
	dummy_var1.insert(pair<string,Real>("DVA",DVA));
	BVA.push_back(dummy_var1);
	
	return BVA;
	
	
}


//to calculate the BVA for the forward with the data of the counterparty and our own
vector<map<string,Real>> BVA_calculator::BVA_cal_AT1P(Date maturity, char collateral_indicator, Real inv_col_LGD, Real cpty_col_LGD)
{
	Real CVA = 0.0;		
	Real DVA = 0.0;
	set_inv_coll_LGD(inv_col_LGD); //reset the investor's collateral LGD level
	set_cpty_coll_LGD(cpty_col_LGD); //reset the counterparty's collateral LGD level

	//generating the default probabilities of counterparty and our own based on the data in the class "At1pCalibration"
	vector<Real> own_at1p_prob = myself_AT1P->AT1P_CDS_calibration();
	vector<Real> counterparty_at1p_prob = cpty_AT1P->AT1P_CDS_calibration();

	//get the Loss given default data for both counterparties
	Real Investor_LGD = 1-own_recovery_rate;
	Real counterparty_LGD = 1-counterparty_recovery_rate;

	// start the for loop to calculate expected EAD
	for (int iPath = 0; iPath < no_of_Date_path; ++iPath) 
	{
		//calculate a sample of default times for both counterparties
		investorDefaultDate = myself_AT1P->Get_default_time_AT1P();
		counterptyDefaultDate = cpty_AT1P->Get_default_time_AT1P();

		//First to default check
		// check if there's any default event before maturity 
		if (investorDefaultDate < maturity || counterptyDefaultDate < maturity)
		{	
			//check if both counterparties default before the term of the forward
			if (investorDefaultDate < maturity && counterptyDefaultDate < maturity)
			{
				if (investorDefaultDate < counterptyDefaultDate) 
				{//if the collateral option is selected ("Y/y"),the montle carlo simulation will include the calculation of the collateral amount
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						DVA = DVA - Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, investorDefaultDate, 0),'d',Investor_LGD,counterparty_LGD);
					else
					//if our own default happens first, calculate the DVA
						DVA = DVA + (Investor_LGD*maximun(-(gen->MCgeneration(no_of_stock_path, investorDefaultDate, 0)),0));
				}
				else
				{//if the collateral option is selected ("Y/y"),the montle carlo simulation will include the calculation of the collateral amount
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						CVA = CVA + (Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, counterptyDefaultDate, 0),'c',Investor_LGD,counterparty_LGD));
					else
					//if counterparty default first, calculate the CVA 
						CVA = CVA + (counterparty_LGD*maximun(gen->MCgeneration(no_of_stock_path, counterptyDefaultDate, 0),0));
				}
			}
			else
			{//check if there is only one of the counterparty default before the maturity of the forward
				if (investorDefaultDate < maturity)
				{
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						DVA = DVA - Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, investorDefaultDate, 0),'d',Investor_LGD,counterparty_LGD);
					else
						DVA = DVA + (Investor_LGD*maximun(-(gen->MCgeneration(no_of_stock_path, investorDefaultDate, 0)),0));
				}
				else
				{	
					if(collateral_indicator == 'Y' || collateral_indicator == 'y')
						CVA = CVA + (Coll_incl_BVA_cal(gen->MCgeneration_collateral(no_of_stock_path, counterptyDefaultDate, 0),'c',Investor_LGD,counterparty_LGD));
					else
						CVA = CVA + (counterparty_LGD*maximun(gen->MCgeneration(no_of_stock_path, counterptyDefaultDate, 0),0));
				}
			}
		}
		// if no default event in path i, then skip to next iterator
		else 
		{
			continue;
		}
	}

	// Expected values for the sampled CVA and DVA
	DVA = DVA / no_of_Date_path;
	CVA = CVA / no_of_Date_path;

	//save the CVA and DVA into the output vector
	map<string,Real> dummy_var;
	vector<map<string,Real>> BVA;
	dummy_var.insert(pair<string,Real>("CVA",CVA));
	BVA.push_back(dummy_var);
	map<string,Real> dummy_var1;
	dummy_var1.insert(pair<string,Real>("DVA",DVA));
	BVA.push_back(dummy_var1);
	
	return BVA;
	
}

//Calculation of the collateral inclusive BVA at the termainal of each simulated path of montle carlo simulation 
Real BVA_calculator::Coll_incl_BVA_cal(pair<Real,Real> forward_pv,char CVA_DVA_ind, Real own_LGD, Real cpt_LGD)
{
	Real CVADVA = 0;
	Real forw = forward_pv.first;
	Real col_v = forward_pv.second;

	//For the equations below, please refer to chapter 13 of Brigo's book on page 314 with the assumption that the exposure is measure in mid-market 
	if(CVA_DVA_ind == 'c') //Collateral inclusive CVA
		CVADVA = cpt_LGD*(maximun((maximun(forw,0)-maximun(col_v,0)),0))+cpty_collateral_LGD*(maximun(minimun(forw,0)-minimun(col_v,0),0));

	if(CVA_DVA_ind == 'd')//Collateral inclusive DVA
		CVADVA = own_LGD*(minimun((minimun(forw,0)-minimun(col_v,0)),0))+own_collateral_LGD*(minimun(maximun(forw,0)-maximun(col_v,0),0));


	return CVADVA;
}

#pragma region BVACalculatorWWR_

// =====================================================================================
// Class: BVACalculatorWWR
// Description: used for calculate UDVA & UCVA
// =====================================================================================

// to calculate the UCVA & UDVA with WWR of the counterparty and investor 
map<string, Real> BVACalculatorWWR::BVA_cal(Date maturity)
{
	Real CVA = 0.0;
	Real DVA = 0.0;
	Real iPath_CVA;
	Real iPath_DVA;
	Real investorDefaultProb;
	Real counterpartyDefaultProb;
	Real v_t_previous;
	Real v_t_current;

	//get the Loss given default data for both counterparties
	Real Investor_LGD = 0.60;
	Real counterparty_LGD = 0.60;
	vector<Date> oDates;
	map<Date, Real> payoffPath;

	// constructing vector<> oDates: assuming product maturity is among the list of CDS tenors 
	for (Size i = 0; i < investorHazardRate_->cds_maturities_.size(); i++) {
		if (investorHazardRate_->cds_maturities_[i] <= maturity) {
			oDates.push_back(investorHazardRate_->cds_maturities_[i]);
		}
		else {
			break;
		}
	}

	// start the for loop to calculate expected UCVA & UDVA
	for (int iPath = 0; iPath < no_of_Date_path; ++iPath)
	{
		payoffPath = oMCgenerator_->SimulatePayoffPath_V1(oDates);
		iPath_DVA = 0.0;
		iPath_CVA = 0.0;

		for (int jDate = 0; jDate < oDates.size(); jDate++) {
			switch (jDate)
			{
			case 0:
				investorDefaultProb = 1 - investorHazardRate_->GetSurvivalProbability(jDate, &payoffPath);
				counterpartyDefaultProb = 1 - counterpartyHazardRate_->GetSurvivalProbability(jDate, &payoffPath);
				v_t_previous = 0;
				break;
			default:
				investorDefaultProb = investorHazardRate_->GetSurvivalProbability(jDate - 1, &payoffPath) - investorHazardRate_->GetSurvivalProbability(jDate, &payoffPath);
				counterpartyDefaultProb = counterpartyHazardRate_->GetSurvivalProbability(jDate - 1, &payoffPath) - counterpartyHazardRate_->GetSurvivalProbability(jDate, &payoffPath);
				v_t_previous = payoffPath.find(oDates[jDate - 1])->second;
				break;
			}
			v_t_current = payoffPath.find(oDates[jDate])->second;

			iPath_DVA += Investor_LGD * investorDefaultProb * std::max(-0.5*(v_t_previous + v_t_current), 0.0);
			iPath_CVA += counterparty_LGD * counterpartyDefaultProb * std::max(0.5*(v_t_previous + v_t_current), 0.0);		// CVA formula according to paper by John Hull with exposure = 0.5 x [v(t-1) + v(t)]
		}

		DVA += iPath_DVA;
		CVA += iPath_CVA;
	}

	// Expected values for the sampled UCVA and UDVA
	DVA = DVA / no_of_Date_path;
	CVA = CVA / no_of_Date_path;

	//save the UCVA and UDVA into the output vector
	map<string, Real> valueAdjustment;
	valueAdjustment["UDVA"] = DVA;
	valueAdjustment["UCVA"] = CVA;
	return valueAdjustment;

}

#pragma endregion

//reset the collateral LGD ratio posted by investor
void BVA_calculator::set_inv_coll_LGD(Real new_inv_coll_LGD)
{
	own_collateral_LGD = new_inv_coll_LGD;
}

//reset the collateral LGD ratio posted by investor
void BVA_calculator::set_cpty_coll_LGD(Real new_cpty_coll_LGD)
{
	cpty_collateral_LGD = new_cpty_coll_LGD;
}