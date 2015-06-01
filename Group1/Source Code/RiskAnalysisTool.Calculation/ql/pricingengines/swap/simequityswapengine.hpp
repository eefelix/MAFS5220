#pragma once

#include <Calculation/Calculation.h>
#include <ql/instruments/swap/equityswap.hpp>
#include <ql/pricingengines/swap/discountingswapengine.hpp>

namespace QuantLib {

    class _RISKANALYSISTOOL_CALCULATION_API SimEquitySwapEngine : public EquitySwap::engine {
    public:
        SimEquitySwapEngine(size_t n_ = 1000); /* default simulate time of the engine is ten thousand */
        ~SimEquitySwapEngine() {
        };
        void calculate() const;
    private:
        size_t n;
		mutable boost::shared_ptr<QuantLib::BlackScholesMertonProcess> process;
    };


}
