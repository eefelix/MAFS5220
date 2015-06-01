#include <pch.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include "cirdefaultmodel.hpp"

//! \file cirdefaultmodel.cpp
//!
using namespace QuantLib;

CirDefaultModel::CirDefaultModel(Rate r0, Real theta,
	Real k, Real sigma)
	: CoxIngersollRoss(r0,theta,k,sigma){
}

const Probability CirDefaultModel::defaultProbability(Time t) const 
{
	return 1 - A(0, t)*std::exp(-B(0, t)*x0());
}
const Time CirDefaultModel::defaultTime(const Path& path) const 
{
	static boost::mt19937 uRng;
	static boost::exponential_distribution<> expDist(1);
	static boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> expRng(uRng, expDist);

	Real objCumIntensity = expRng();
	Real realizedCumIntensity = 0.0;
	TimeGrid time = path.timeGrid();
	Time previousTime = 0.;

	for (auto itr = time.begin(); itr != time.end(); ++itr){
		realizedCumIntensity += path[itr - time.begin()]*(*itr-previousTime);
		previousTime = *itr;

		if (realizedCumIntensity >= objCumIntensity)
			return previousTime;
	}

	return Null<Time>();
}

boost::shared_ptr<StochasticProcess1D> CirDefaultModel::process() const
{
	process_ = boost::make_shared<CIRprocess>(theta(), k(), sigma(), x0() );
	return process_;
}

void CirDefaultModel::calibrate(
	const std::vector<boost::shared_ptr<CalibrationHelper> >& helper,
	OptimizationMethod& method,
	const EndCriteria& endCriteria,
	const Constraint& constraint/* = Constraint()*/,
	const std::vector<Real>& weights/* = std::vector<Real>()*/,
	const std::vector<bool>& fixParameters /*= std::vector<bool>()*/)
{
	CalibratedModel::calibrate(helper, method, endCriteria, constraint, weights, fixParameters);
}