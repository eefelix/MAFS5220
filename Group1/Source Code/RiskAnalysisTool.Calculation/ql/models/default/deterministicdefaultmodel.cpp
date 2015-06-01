#include <pch.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <ql/models/defaultmodel.hpp>
#include <ql/methods/montecarlo/path.hpp>

#include "deterministicdefaultmodel.hpp"

//! \file deterministicdefaultmodel.cpp
//!
using namespace QuantLib;

const Time DeterministicDefaultModel::defaultTime(const Path& path) const{
	static boost::mt19937 uRng;
	static boost::exponential_distribution<> expDist(1);
	static boost::variate_generator<boost::mt19937&, boost::exponential_distribution<>> expRng(uRng, expDist);

	Real objCumIntensity = expRng();
	Real realizedCumIntensity = 0.0;
	TimeGrid time = path.timeGrid();
	Time previousTime = 0.;

	for (auto itr = time.begin(); itr != time.end(); ++itr){
		realizedCumIntensity += hazardRateStructure_->defaultDensity(*itr,true) * (*itr - previousTime);
		previousTime = *itr;

		if (realizedCumIntensity >= objCumIntensity)
			return previousTime;
	}

	return Null<Time>();
}