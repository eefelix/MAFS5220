#pragma once 

#include <ql/types.hpp>
#include <ql/math/copulas/gaussiancopula.hpp>
#include <ql/math/statistics/statistics.hpp>
//! \file exactsimulation.hpp

namespace QuantLib {

    template <template <class> class SIM, class RNG, class S = Statistics, typename C = GaussianCopula>
    class ExactSimulation {
        inline void calculate(Time t) const {
            QL_FAIL("Not implemented!");
        }


    };
}