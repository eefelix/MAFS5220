#include <pch.h>
#include <boost/math/distributions/normal.hpp>
#include <ql/exercise.hpp>
#include "analytickoueuropeanengine.hpp"
//! \file analytickoueuropeanengine.cpp

using namespace std;
using namespace QuantLib;

AnalyticKouEuropeanEngine::AnalyticKouEuropeanEngine(const boost::shared_ptr<KouProcess>&process, const Real &tolerance /* = 0.00001 */)
	: process_(process), tolerance_(tolerance)
{
	registerWith(process_);
}

void AnalyticKouEuropeanEngine::calculate() const
{

	QL_REQUIRE(arguments_.exercise->type() == Exercise::European,
		"not an European option");

	boost::shared_ptr<StrikedTypePayoff> payoff =
		boost::dynamic_pointer_cast<StrikedTypePayoff>(arguments_.payoff);

	QL_REQUIRE(payoff, "non-striked payoff given");

	/*
	*calculate to be implemented
	*/
	Option::Type type = payoff->optionType();
	Real eta1 = process_->posJumpMean();
	QL_REQUIRE(eta1 > 1, "invalid eta1 input!!");
	Real eta2 = process_->negJumpMean();
	QL_REQUIRE(eta2 > 0, "invalid eta2 input!!");
	Real lambda = process_->jumpIntensity();
	Real p = process_->posProbability();

	//here assume the three term structure have the same daycounter
	DayCounter daycounter = process_->riskFreeRate()->dayCounter();
	Real sigma = process_->blackVolatility()->blackVol(arguments_.exercise->lastDate(), payoff->strike());
	Real variance = sigma*sigma;
	// the paper ignore the case when dividends are payed out
	Real effectiveR =
		process_->riskFreeRate()->zeroRate(arguments_.exercise->lastDate(), daycounter, Continuous)
		- process_->dividendYield()->zeroRate(arguments_.exercise->lastDate(), daycounter, Continuous);
	Real S = process_->stateVariable()->value();
	QL_REQUIRE(S > 0.0, "negative or null underlying given");
	Real K = payoff->strike();
	Time T = daycounter.yearFraction(process_->riskFreeRate()->referenceDate(),
		arguments_.exercise->lastDate());

	/*
	*auxillary variables
	*/
	Real zeta = p*eta1 / (eta1 - 1) + (1 - p)*eta2 / (eta2 + 1) - 1;
	Real adjustEta1 = eta1 - 1;
	Real adjustEta2 = eta2 + 1;
	Real adjustLambda = process_->jumpIntensity()*(zeta + 1);
	Real effectiveP = p / (1 + zeta)*eta1 / (eta1 - 1);

	//price of call option
	Real price =
		S*process_->dividendYield()->discount(arguments_.exercise->lastDate())
		*KouHelper::Gamma(effectiveR + 0.5*sigma*sigma - lambda*zeta, sigma, adjustLambda, effectiveP, adjustEta1, adjustEta2, log(K / S), T, tolerance_) -
		K*process_->riskFreeRate()->discount(arguments_.exercise->lastDate())
		*KouHelper::Gamma(effectiveR - 0.5*sigma*sigma - lambda*zeta, sigma, lambda, p, eta1, eta2, log(K / S), T, tolerance_);

	if (type == Option::Type::Call)
	{
		results_.value = price;
	}

	else results_.value = price - S*process_->dividendYield()->discount(arguments_.exercise->lastDate())
		+ K*process_->riskFreeRate()->discount(arguments_.exercise->lastDate());

	//more results if needed
}

const int KouHelper::combination(int n, int k)
{
	QL_REQUIRE(k >= 0 && k <= n && n >= 0,
		"k should between 0 and n, n should not be smaller than zero.");

	return int(Factorial::get(n) / (Factorial::get(n - k)*Factorial::get(k)));
}

const double KouHelper::Hh(int n, double x)
{
	QL_REQUIRE(n >= -1,
		"n should be a integer not less than minus one.");

	boost::math::normal Normal;
	if (n == -1)
		return sqrt(2 * M_PI)*boost::math::pdf(Normal, x);
	else if (n == 0)
		return sqrt(2 * M_PI)*boost::math::cdf(Normal, -x);
	else
		return (Hh(n - 2, x) - x*Hh(n - 1, x)) / n;
}

const double KouHelper::P(int n, int k, double eta1, double eta2, double p)
{
	QL_REQUIRE(k >= 1 && k <= n && n >= 1,
		"k should between 1 and n-1.");
	double s = 0.;

	double a = eta1 / (eta1 + eta2);
	double b = 1 - a;
	if (k == n)
	{
		return pow(p, n);
	}
	else
	{
		for (int i = k; i <= n - 1; i++)
		{
			s += combination(n - k - 1, i - k)*combination(n, i)*pow(a, i - k)
				*pow(b, n - i)*pow(p, i)*pow(1 - p, n - i);
		}
		return s;
	}
}

const double KouHelper::Q(int n, int k, double eta1, double eta2, double p)
{
	QL_REQUIRE(k >= 1 && k <= n && n >= 1,
		"k should between 1 and n-1.");
	double s = 0.0;

	double a = eta1 / (eta1 + eta2);
	double b = 1 - a;
	if (k == n) {
		return pow(1 - p, n);
	}
	else
	{
		for (int i = k; i <= n - 1; i++)
		{
			s += combination(n - k - 1, i - k)*combination(n, i)*pow(a, n - i)*pow(b, i - k)
				*pow(p, n - i)*pow(1 - p, i);
		}
		return s;
	}
}

const double KouHelper::I(int n, double c, double alpha, double beta, double delta)
{
	double S = 0.;
	boost::math::normal Normal;
	if (beta > 0 && alpha != 0 && n >= -1)
	{
		for (int i = 0; i <= n; i++)
		{
			S = S + pow(beta / alpha, n - i)*Hh(i, beta*c - delta);
		}
		S = -S*exp(alpha*c) / alpha;
		S = S + pow(beta / alpha, n + 1)*sqrt(2. * M_PI)
			/ beta*exp(alpha*delta / beta + alpha*alpha / 2. / beta / beta)
			*boost::math::cdf(Normal, -beta*c + delta + alpha / beta);
	}
	else if (beta < 0 && alpha < 0 && n >= -1)
	{
		for (int i = 0; i <= n; i++)
		{
			S = S + pow(beta / alpha, n - i)*Hh(i, beta*c - delta);
		}
		S = -S*exp(alpha*c) / alpha;
		S = S - pow(beta / alpha, n + 1)*sqrt(2. * M_PI)
			/ beta*exp(alpha*delta / beta + alpha*alpha / 2. / beta / beta)
			*boost::math::cdf(Normal, beta*c - delta - alpha / beta);
	}
	return S;
}

const double KouHelper::Pi(int n, double lambda, double T)
{
	QL_REQUIRE(n >= 0 && T > 0 && lambda >= 0,
		"parameters invalid!");
	return exp(-lambda*T)*pow(lambda*T, n) / Factorial::get(n);
}

const double KouHelper::Gamma(
	double mu, double sigma, double lambda,
	double p, double eta1, double eta2, double a,
	double T, double tolerance)
{
	double S = 0;
	double s, s1, s2;
	double t = sigma*sqrt(T);
	boost::math::normal Normal;
	int n = 1;
	double a1 = exp(0.5*eta1*eta1*t*t) / (t*sqrt(2 * M_PI));
	double a2 = exp(0.5*eta2*eta2*t*t) / (t*sqrt(2 * M_PI));

	S = S + Pi(0, lambda, T)*boost::math::cdf(Normal, -(a - mu*T) / t);
	do
	{
		s = 0;
		s1 = 0;
		s2 = 0;
		for (int k = 1; k <= n; k++)
		{
			s1 = s1 + P(n, k, eta1, eta2, p)*pow(t*eta1, k)*I(k - 1, a - mu*T, -eta1, -1 / t, -eta1*t);
		}
		s1 = s1*Pi(n, lambda, T)*a1;
		for (int k = 1; k <= n; k++)
		{
			s2 = s2 + Q(n, k, eta1, eta2, p)*pow(t*eta2, k)*I(k - 1, a - mu*T, eta2, 1 / t, -eta2*t);
		}
		s2 = s2*Pi(n, lambda, T)*a2;
		s = s1 + s2;
		S = S + s;
		n++;
	} while (abs(s) > S * tolerance);
	return S;
}