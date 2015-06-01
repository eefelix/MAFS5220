
#include "pch.h"
#include "CIRcalibration.hpp"

void CIRProblemFunction::Initialization(std::vector<Real> quoteVector)
{
	termNum = quoteVector.size();
	quote.clear();
	quote.insert(quote.end(), quoteVector.begin(), quoteVector.end());

	Real sum = 0;
	for (int i = 0; i != quoteVector.size(); i++)
		sum += quoteVector[i];
	gamma = sum / quoteVector.size();

	for (int i = 0; i != quoteVector.size(); i++)
		rt.push_back(quoteVector[i] - gamma);
}

Real CIRProblemFunction::Calculation(Real rt_, Real rt_1, Real phi1, Real gamma1) const{
	return (rt_ - phi1*rt_1)*(rt_ - phi1*rt_1) / (rt_ + gamma1);
}

Real CIRProblemFunction::value(const Array& x) const
{
	QL_REQUIRE(x.size() == 1, "phi is 1-dim");
	Real phi = x[0];
	Real sum = 0.0;
	Real tmp;
	for (int i = 1; i < termNum; i++){
		tmp = Calculation(rt[i], rt[i - 1], phi, gamma);
		sum += tmp;
	}
	return sum;
}

Disposable<Array> CIRProblemFunction::values(const Array& x) const
{
	QL_REQUIRE(x.size() == 1, "phi is 1-dim");
	Array res(1);
	res[0] = value(x);
	return res;
}


void CIRcalibration::Initialization(std::vector<Real> quoteVector)
{
	optimizeFunc.Initialization(quoteVector);
	termNum = quoteVector.size();
	quote.clear();
	quote.insert(quote.end(), quoteVector.begin(), quoteVector.end());

	this->r0 = quoteVector.back();
	this->gamma = optimizeFunc.getGamma();
}

void CIRcalibration::StartCalibrate()
{
	Size maxIterations = 1000;
	Size minStatIterations = 100;
	Real rootEpsilon = 1e-12;
	Real functionEpsilon = 1e-12;
	Real gradientNormEpsilon = 1e-12;
	EndCriteria myEndCrit(maxIterations, minStatIterations, rootEpsilon,
		functionEpsilon, gradientNormEpsilon);

	Array pi(1);
	pi[0] = 0.01;

	NoConstraint constraint;
	Problem CIRcalibrationProblem(optimizeFunc, constraint, pi);
	Simplex solver(0.1);
	EndCriteria::Type solvedCrit = solver.minimize(CIRcalibrationProblem, myEndCrit);

	Array resultedVol = CIRcalibrationProblem.currentValue();
	phi = resultedVol[0];
	RSS = CIRcalibrationProblem.functionValue();
}

void CIRcalibration::CalculateParameters()
{
	DiscreteSigma = sqrt(RSS / (termNum - 2));
	k = log(1 / phi);
	ContinousSigma = sqrt(DiscreteSigma*DiscreteSigma * 2 * k / (1 - exp(-2 * k)));
}