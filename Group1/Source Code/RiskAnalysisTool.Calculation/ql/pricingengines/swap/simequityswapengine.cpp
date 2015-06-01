#include <pch.h>
#include "simequityswapengine.hpp"
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/math/randomnumbers/boxmullergaussianrng.hpp>
#include <ql/cashflows/dividend.hpp>

using namespace QuantLib;

SimEquitySwapEngine::SimEquitySwapEngine(size_t n_) : n(n_) {
	
}

void SimEquitySwapEngine::calculate() const {
	process = boost::make_shared<QuantLib::BlackScholesMertonProcess>(
		QuantLib::Handle<QuantLib::Quote>(boost::make_shared<QuantLib::SimpleQuote>(arguments_.spotPrice)),
		QuantLib::Handle<QuantLib::YieldTermStructure>(
		boost::make_shared<QuantLib::FlatForward>(arguments_.referenceDate, arguments_.dividend, arguments_.discountCurve->dayCounter())),
		QuantLib::Handle<QuantLib::YieldTermStructure>(arguments_.discountCurve),
		QuantLib::Handle<QuantLib::BlackVolTermStructure>(
			boost::make_shared<QuantLib::BlackConstantVol>(arguments_.referenceDate, arguments_.discountCurve->calendar(), arguments_.sigma, arguments_.discountCurve->dayCounter())));

    double dividend = arguments_.dividend;
    double leg1 = 0, leg2 = 0;
    double s = 0;
    double xt = arguments_.spotPrice;
    double dt = 0.004;
    double t;
    double t0;
    double delta;

	
    long seed1 = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng rnd1(seed1);
    BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> normrnd1(rnd1);
    double dw = 0;

    std::vector<double> floatingleg(arguments_.count, 0);
    // simulate ten thousand times
    for (size_t k = 1; k <= n; ++k) {
		xt = arguments_.spotPrice;
        /* for every path simulate the underlying process
        * and calculate the dividend payment every section
        */
        for (size_t i = 0; i < arguments_.count; ++i) {
            s = 0;
            t = arguments_.yearfraction[i + 1];
            t0 = arguments_.yearfraction[i];

            do {
                dw = normrnd1.next().value;
                xt = process->evolve(t0, xt, dt, dw);
                delta = std::min(dt, t - t0);
                s += xt*delta*dividend;
                t0 = std::min(t0 + dt, t);
            } while (t0 < (t - 1e-6));
            if (i == 0) {
                floatingleg[i] += (arguments_.cumulativeDividend + s);
            } else if (i < arguments_.count - 1) {
                floatingleg[i] += s;
            } else floatingleg[i] += (s + xt);
        }

    }

    for (size_t i = 0; i < arguments_.count; ++i) {
        arguments_.legs[0].push_back(boost::make_shared<QuantLib::FixedDividend>(floatingleg[i] / n, arguments_.schedule->at(i)));
        leg1 += arguments_.legs[0].at(i)->amount()*arguments_.discountCurve->discount(arguments_.legs[0].at(i)->date());
        leg2 += arguments_.legs[1].at(i)->amount()*arguments_.discountCurve->discount(arguments_.legs[1].at(i)->date());
    }
    results_.legNPV.push_back(leg1);
    results_.legNPV.push_back(leg2);
    if (arguments_.type == EquitySwap::Type::Payer) {
        results_.value = arguments_.amount*(leg1 - leg2);
    } else results_.value = arguments_.amount*(leg2 - leg1);
}