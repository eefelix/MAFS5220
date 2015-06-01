/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Xiang, GAO

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
#include <ql/types.hpp>
#include <ql/methods/montecarlo/mctraits.hpp>
//! \file multivaluemctraits.hpp
//! \brief Traits class for multi-state Monte Carlo Simulation. 
//! 
namespace QuantLib {
	//! \brief Traits class for multi-state Monte Carlo Simulation. 
	//! 
    template <class RNG>
    struct MultiValueMultiVariate {
        typedef RNG rng_traits;
        typedef MultiPath path_type;
        typedef PathPricer<path_type, Array> path_pricer_type;
        typedef typename RNG::rsg_type rsg_type;
        typedef MultiPathGenerator<rsg_type> path_generator_type;
        enum {
            allowsErrorEstimate = RNG::allowsErrorEstimate
        };
    };
}