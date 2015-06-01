#include "pch.h"
#include"CIRprocess.hpp"

using namespace QuantLib;

CIRprocess::CIRprocess(Rate r01, Real theta1, Real k1, Real sigma1){
	r0 = r01;
	theta = theta1;
	k = k1;
	sigma = sigma1;
}

Real CIRprocess::P(Date startDate, Date endDate, Rate rt) const{
	Real value = A(startDate, endDate)*std::exp(-B(startDate, endDate)*rt);
	return value;
}

Real CIRprocess::D(Date startDate, Date endDate) const{
	BigInteger seed = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMt(seed);
	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGauss(unifMt);

	Real dt = dc.yearFraction(startDate, startDate + 1 * Days);
	int n = endDate - startDate;
	Real interestIntegral = 0.0;
	Real rt = r0;
	Real temp1 = r0;

	//int n = endDate - startDate;

	for (int i = 0; i < n; i++)
	{
		interestIntegral += rt*dt;
		rt += k*(theta - temp1)*dt + sigma*std::sqrt(temp1)*bmGauss.next().value*std::sqrt(dt);
		temp1 = rt;
	}
	Real  discountFactor = std::exp(-interestIntegral);
	return discountFactor;
}

Real CIRprocess::A(Date startDate, Date endDate) const {
	Real sigma2 = sigma*sigma;
	Real h = std::sqrt(k*k + 2.0*sigma2);
	Real numerator = 2.0*h*std::exp(0.5*(k + h)*(endDate - startDate));
	Real denominator = 2.0*h + (k + h)*(std::exp((endDate - startDate)*h) - 1.0);
	Real value = std::log(numerator / denominator)*
		2.0*k*theta / sigma2;
	return std::exp(value);
}

Real CIRprocess::B(Date startDate, Date endDate) const {
	Real h = std::sqrt(k*k + 2.0*sigma*sigma);
	Real temp = std::exp((endDate - startDate)*h) - 1.0;
	Real numerator = 2.0*temp;
	Real denominator = 2.0*h + (k + h)*temp;
	Real value = numerator / denominator;
	return value;
}