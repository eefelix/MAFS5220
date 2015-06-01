#include "pch.h"
#include "IRSUCVAWithWWR.hpp"

using namespace QuantLib;

//initialize fixedLeg and floatingLeg vectors, which store cashflow exchange date
//the first element of each vector is the startDate, which has no cashflow exchange
//the cashflow exchange starts from the second element
void CVAwithWWR::initialization(){
	for (int i = 0; i <= fixedFrequency*tenor.length(); i++)
		fixedLeg.push_back(startDate + i*Years / fixedFrequency);
	for (int i = 0; i <= floatingFrequency*tenor.length(); i++)
		floatingLeg.push_back(startDate + i*Years / floatingFrequency);
	endDate = startDate + tenor;
}

// A and B are used to calculate P, which is the discount factor
Real CVAwithWWR::A(Date t1, Date t2){
	Real sigma12 = sigma1*sigma1;
	Real h = std::sqrt(k*k+2.0*sigma12);
	Real numerator = 2.0*h*std::exp(0.5*(k + h)*dc.yearFraction(t1,t2));
	Real denominator = 2.0*h + (k + h)*(std::exp(dc.yearFraction(t1,t2)*h) - 1.0);
	Real value = std::log(numerator / denominator)*
		2.0*k*theta / sigma12;
	return std::exp(value);
}

Real CVAwithWWR::B(Date t1, Date t2){
	Real h = std::sqrt(k*k + 2.0*sigma1*sigma1);
	Real temp = std::exp(dc.yearFraction(t1,t2)*h) - 1.0;
	Real numerator = 2.0*temp;
	Real denominator = 2.0*h + (k + h)*temp;
	Real value = numerator / denominator;
	return value;
}

Real CVAwithWWR::P(Date t1, Date t2, Real rt){
	return A(t1, t2)*exp(-B(t1,t2)*rt);
}

//calculate the npv of time tao
Real CVAwithWWR::npv(Date tao){
	Real npv = 0.0;
	//firstFixedIndex is used to store the index of the first cashflow exchange date after default of FixedLeg
	//firstFloatingIndex is used to store the index of the first cashflow exchange date after default of FloatingLeg
	int firstFixedIndex=0;
	int firstFloatingIndex = 0;
	//fixednpv stores the new present value of remaining cash flow of fixed leg at default time tao
	Real fixednpv = 0.0;
	//floatingnpv stores the new present value of remaining cash flow of floating leg at default time tao
	Real floatingnpv = 0.0;

	for (Size i = 1; i < fixedLeg.size(); i++)
	{
		if (fixedLeg[i] >= tao){
			firstFixedIndex = i;
			break;
		}
	}

	for (Size i = firstFixedIndex; i < fixedLeg.size(); i++){
		fixednpv += fixedRate*P(tao, fixedLeg[i], rtao)*dc.yearFraction(fixedLeg[i - 1], fixedLeg[i]);
	}

	for (Size i = 1; i < floatingLeg.size(); i++)
	{
		if (floatingLeg[i] >= tao){
			firstFloatingIndex = i;
			break;
		}
	}

	//the libor part of floatingLeg 
	floatingnpv += P(tao, floatingLeg[firstFloatingIndex], rtao) - P(tao, floatingLeg[floatingLeg.size()-1], rtao);
	//the spread part of floatingLeg
	for (Size i = firstFloatingIndex; i < floatingLeg.size(); i++){
		floatingnpv += P(tao, floatingLeg[i], rtao)*spread*dc.yearFraction(floatingLeg[i - 1], floatingLeg[i]);
	}

	if (type == Payer)
		return floatingnpv - fixednpv;
	else return fixednpv - floatingnpv;
}

//calculate the ucva of single path 
Real CVAwithWWR::ucvaPath(Date referenceDate, Date endDate){
	BigInteger seed = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMt(seed);
	//expRandomNum is a random number conforming exponential distribution, which is used to calculate the default time
	Real expRandomNum = -log(1 - unifMt.next().value);

	seed = SeedGenerator::instance().get();
	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGauss(unifMt);

	BigInteger seed2 = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMt2(seed2);
	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGauss2(unifMt2);

	Real dt = dc.yearFraction(referenceDate, referenceDate + 1*Days);
	int n = endDate - referenceDate;
	Real interestIntegral = 0.0;
	Real hazardIntegral = 0.0;

	Real rt = r0;
	Real temp1 = r0;

	Real temp2 = lamda0;
	Real lamdat = lamda0;
	Real discount = 0.0;

	bool defaulting = false;
	Date defaultTime;
	Calendar calendar = TARGET();

	for (int i = 0; i < n; i++){
		if (!calendar.isBusinessDay(referenceDate + i * Days))
			continue;
		if (hazardIntegral<expRandomNum){
			Real rm1 = bmGauss.next().value;
			Real rm2 = bmGauss2.next().value;
			//interestIntegral is used to calculate the discount factor from reference date to default time
			interestIntegral += rt*dt;
			rt += k*(theta - temp1)*dt + sigma1*std::sqrt(temp1)*rm1*std::sqrt(dt);
			temp1 = rt;

			hazardIntegral += lamdat*dt;
			lamdat += a*(b - temp2)*dt + sigma2*std::sqrt(temp2)*(rho*rm1 + sqrt(1 - rho*rho)*rm2)*std::sqrt(dt);
			temp2 = lamdat;
		}
		//when default happens, it turns to the following else branch
		else{
			discount = std::exp(-interestIntegral);
			//rtao is used to calculate P
			rtao = rt;
			defaulting = true;
			defaultTime = referenceDate + i * Days;
			break;
		}
	}

	if (!defaulting) { return 0.0; }
	// when default, it turns to the following branch
	else if (defaultTime > endDate) { return 0.0; }
	else
	{
		if (npv(defaultTime) > 0)
		{
			return  (1 - recoveryRate) *nominal* npv(defaultTime) * discount;
		}
		else { return 0.0; }
	}
}

//simulate paths whose number is pathNum
Real CVAwithWWR::cvaWithWWRCalculation(int pathNum){
	initialization();
	Real ucva = 0.0;
	for (int i = 0; i < pathNum; i++) {
		ucva += ucvaPath(referenceDate, endDate);
	}
	return ucva / pathNum;
}

