/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Fang Xiong

This file is part of QuantLib, a free-software/open-source library
for financial quantitative analysts and developers - http://quantlib.org/

QuantLib is free software: you can redistribute it and/or modify it
under the terms of the QuantLib license.  You should have received a
copy of the license along with this program; if not, please email
<quantlib-dev@lists.sf.net>. The license is also available online at
<http://quantlib.org/license.shtml>.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#pragma once
#include <Calculation/Calculation.h>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/processes/eulerdiscretization.hpp>
#include <ql/math/distributions/poissondistribution.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <ql/math/randomnumbers/seedgenerator.hpp>
//! \file jcirprocess.hpp
//! \brief JCIR process class
//!
namespace QuantLib {
    //!	Jump diffusion CIR process class
    /*!	This class describes a jump diffusion CIR process governed by
		\f[
		dx = a (b - x_t) dt + \sigma \sqrt{x_t} dW_t + dJ_t.
		\f]
		It is an extension of the classical CIR model. The jumps of the JCIR are introduced with the help
		of a pure-jump Levy process \f$J_t\f$.
    */
    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    class /*_RISKANALYSISTOOL_CALCULATION_API*/ GeneralizedJcirProcess : public StochasticProcess1D
    {
    public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of CIR process using the four parameters: long-term mean, 
		//! revert speed, volatility, initial value and the additional jump information
        GeneralizedJcirProcess(
            Real mean, Real speed, Real jumpIntensity, Real jumpMean, Volatility sigma, Real x0 = 0.0,
            const __URng_Poisson_Type &URng_Poisson = __URng_Poisson_Type(SeedGenerator::instance().get()),
            const __URng_Exp_Type &URng_Exp = __URng_Poisson_Type(SeedGenerator::instance().get()),
            const boost::shared_ptr<discretization> &d = boost::make_shared<EulerDiscretization>()
            );
		//@}

		//! \name Public interface
		//{@
		//! Return the drift term of the process
        const Real drift(Time t, Real x) const;
		//! Return the diffusion of the process
        const Real diffusion(Time t, Real x) const;
		//! Return simulation value of the process at time t0+dt
        const Real evolve(Time t0, Real x0, Time dt, Real dw) const;

        const Real x0() const;
        const Real mean() const;
        const Real speed() const;
        const Real jumpIntensity() const;
        const Real jumpMean() const;
        const Volatility volatility() const;
		//@}
    private:
        Real x0_, mean_, speed_, jumpIntensity_, jumpMean_;
        Volatility volatility_;
        __URng_Poisson_Type	URng_Poisson_;
        __URng_Exp_Type	URng_Exp_;
    };

    typedef GeneralizedJcirProcess<MersenneTwisterUniformRng, MersenneTwisterUniformRng> JcirProcess;

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::GeneralizedJcirProcess(
        Real mean, Real speed, Real jumpIntensity, Real jumpMean, Volatility sigma, 
        Real x0,
        const __URng_Poisson_Type &URng_Poisson,
        const __URng_Exp_Type &URng_Exp,
        const boost::shared_ptr<discretization>& d
        ) 
        : StochasticProcess1D(d), mean_(mean), speed_(speed), jumpIntensity_(jumpIntensity)
        , jumpMean_(jumpMean), volatility_(sigma), x0_(x0) {
        // Empty
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::x0() const {
        return x0_;
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::mean() const {
        return mean_;
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::speed() const {
        return speed_;
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::jumpMean() const {
        return jumpMean_;
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::jumpIntensity() const {
        return jumpIntensity_;
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Volatility GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::volatility() const {
        return volatility_;
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::drift(Time, Real x) const {
        return speed_*(mean_ - x);
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::diffusion(Time, Real x) const {
        return volatility_*std::sqrt(x);
    }

    template<typename __URng_Poisson_Type, typename __URng_Exp_Type>
    inline const Real GeneralizedJcirProcess<__URng_Poisson_Type, __URng_Exp_Type>::evolve(Time t0, Real x0, Time dt, Real dw) const {
        Real xContinuous, totalJump = 0.;
        xContinuous = apply(expectation(t0, x0, dt), stdDeviation(t0, x0, dt)*dw);

        InverseCumulativePoisson invP(jumpIntensity_*dt);
        int Pt = (int)invP(URng_Poisson_.next().value);

        for (int i = 0; i < Pt; ++i) {
            totalJump += -mean*std::log(1 - URng_Exp_.next().value);
        }

        return apply(xContinuous, totalJump);
    }
}