#include "pch.h"
#include "IRSBVA.hpp"

using namespace QuantLib;

void IRSBVA::initialization(){
	for (int i = 0; i <= fixedFrequency*tenor.length(); i++)
		fixedLeg.push_back(startDate + i*Years / fixedFrequency);
	for (int i = 0; i <= floatingFrequency*tenor.length(); i++)
		floatingLeg.push_back(startDate + i*Years / floatingFrequency);
	endDate = startDate + tenor;
}

Real IRSBVA::A(Date t1, Date t2){
	Real sigma12 = sigmar*sigmar;
	Real h = std::sqrt(ar*ar + 2.0*sigma12);
	Real numerator = 2.0*h*std::exp(0.5*(ar + h)*dc.yearFraction(t1, t2));
	Real denominator = 2.0*h + (ar + h)*(std::exp(dc.yearFraction(t1, t2)*h) - 1.0);
	Real value = std::log(numerator / denominator)*
		2.0*ar*br / sigma12;
	return std::exp(value);

}

Real IRSBVA::B(Date t1, Date t2){
	Real h = std::sqrt(ar*ar + 2.0*sigmar*sigmar);
	Real temp = std::exp(dc.yearFraction(t1, t2)*h) - 1.0;
	Real numerator = 2.0*temp;
	Real denominator = 2.0*h + (ar + h)*temp;
	Real value = numerator / denominator;
	return value;
}

Real IRSBVA::P(Date t1, Date t2, Real rt){
	return A(t1, t2)*exp(-B(t1, t2)*rt);
}

Real IRSBVA::npv(Date tao){
	Real npv = 0.0;
	int firstFixedIndex = 0;
	int firstFloatingIndex = 0;
	Real fixednpv = 0.0;
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

	floatingnpv += P(tao, floatingLeg[firstFloatingIndex], rtao) - P(tao, floatingLeg[floatingLeg.size() - 1], rtao);

	for (Size i = firstFloatingIndex; i < floatingLeg.size(); i++){
		floatingnpv += P(tao, floatingLeg[i], rtao)*spread*dc.yearFraction(floatingLeg[i - 1], floatingLeg[i]);
	}

	if (type == Payer)
		return floatingnpv - fixednpv;
	else return fixednpv - floatingnpv;
}

Real IRSBVA::bvaPath(Date referenceDate, Date endDate){
	BigInteger seedb = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMtb(seedb);
	Real expRandomNumb = -log(1 - unifMtb.next().value);

	BigInteger seedc = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMtc(seedc);
	Real expRandomNumc = -log(1 - unifMtc.next().value);

	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGaussb(unifMtb);
	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGaussc(unifMtc);

	BigInteger seedr = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMtr(seedr);
	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGaussr(unifMtr);

	Real dt = dc.yearFraction(referenceDate, referenceDate + 1 * Days);
	int n = endDate - referenceDate;
	Real interestIntegral = 0.0;
	Real hazardbIntegral = 0.0;
	Real hazardcIntegral = 0.0;

	Real rt = r0;
	Real tempr = r0;
	Real discount = 0.0;

	Real tempb = lamda0b;
	Real lamdatb = lamda0b;

	Real tempc = lamda0c;
	Real lamdatc = lamda0c;

	bool defaultingb = false;
	bool defaultingc = false;
	Date defaultTime;
	Calendar calendar = TARGET();

	//the following parameter is used for cholesky decomposition
	Real a21 = rhobr;
	Real a22 = sqrt(1-a21*a21);
	Real a31 = rhocr;
	Real a32 = (rhobc - a31*a21) / sqrt(1-a21*a21);
	Real a33 = sqrt(1-rhocr*rhocr-a32*a32);

	for (int i = 0; i < n; i++){
		if (!calendar.isBusinessDay(referenceDate + i * Days))
			continue;
		if (hazardbIntegral<expRandomNumb && hazardcIntegral<expRandomNumc){
			//the following three parameters are independent normal random number
			Real rmr = bmGaussr.next().value;
			Real rmb = bmGaussb.next().value;
			Real rmc = bmGaussc.next().value;
			interestIntegral += rt*dt;
			rt += ar*(br - tempr)*dt + sigmar*std::sqrt(tempr)*rmr*std::sqrt(dt);
			tempr = rt;

			hazardbIntegral += lamdatb*dt;
			lamdatb += ab*(bb - tempb)*dt + sigmab*std::sqrt(tempb)*(a21*rmr + a22*rmb)*std::sqrt(dt);
			tempb = lamdatb;

			hazardcIntegral += lamdatc*dt;
			lamdatc += ac*(bc - tempc)*dt + sigmac*std::sqrt(tempc)*(a31*rmr + a32*rmb+a33*rmc)*std::sqrt(dt);
			tempc = lamdatc;
		}
		//if B defaults
		else if (hazardbIntegral>=expRandomNumb)
		{
			discount = std::exp(-interestIntegral);
			defaultingb = true;
			rtao = rt;
			defaultTime = referenceDate + i * Days;
			break;
		}
		//if Counterparty C defaults
		else{
			discount = std::exp(-interestIntegral);
			defaultingc = true;
			rtao = rt;
			defaultTime = referenceDate + i * Days;
			break;
		}
	}
	
	// if no default happens
	if ((!defaultingb) && (!defaultingc)) { return 0.0; }
	//if B defaults
	else if (defaultTime > endDate){ return 0.0; }
	else if (defaultingb){
		Real bva = 0.0;
		//To calculat DVA, we only include npv that is negative
		if (npv(defaultTime)<0)
			bva = (1 - recoveryRateb) *(-npv(defaultTime)) * discount;
		return bva;
	}
	//if C defaults
	else{
		Real bva = 0.0;
		//To calculat CVA, we only include npv that is positive
		if (npv(defaultTime)>0) 
			bva = -(1 - recoveryRatec) * npv(defaultTime) * discount;
		return bva;
	}
}

Real IRSBVA::bvaCalculation(int pathNum){
	initialization();
	Real bva=0.0;
	for (int i = 0; i < pathNum; i++) {
		 bva += bvaPath(referenceDate, endDate);
	}
	return bva / pathNum;
}