#pragma once

#include <Calculation/Calculation.h>
#include <boost/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/models/model.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/methods/montecarlo/multipath.hpp>

namespace QuantLib {
	class _RISKANALYSISTOOL_CALCULATION_API ExposureModel  {
	public:
		virtual Real NPV(const MultiPath& path, Time defaultTime) const = 0;
		virtual std::vector<boost::shared_ptr<StochasticProcess1D>> processes() const = 0;
	};
}
