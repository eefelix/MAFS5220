/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
Copyright (C) 2015 Tian, XIE

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
#include <memory>
#include <vector>
#include <ql/math/optimization/all.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
//! \file kouprocesscalibrator.hpp
//! \brief Kou process calibrator class
//!
namespace QuantLib {
	//!	A Kou process calibrator class which is used to do calibration using sets of vanilla options
	/*!	This class is used to calibrate kou process using sets of vanilla options with different tenors
		on the underlying stock, it is implemented by inheriting the CostFunction and using the
		Optimizing framework in QuantLib.
	*/
    class _RISKANALYSISTOOL_CALCULATION_API KouProcessCalibrator : public QuantLib::CostFunction {
    public:
		//! \name Constructors & Destructors
		//{@
		//! Creates an instance of Kou process calibrator using sets of vanilla options information
		//! and some market observable information of underlying stock
        KouProcessCalibrator(
            const QuantLib::Date& referenceDate, const QuantLib::Calendar& calendar, const QuantLib::DayCounter& dayCounter,
            double riskFreeRate, double spotPrice, double dividend,
            std::shared_ptr<std::vector<double>> strike,
            std::shared_ptr<std::vector<QuantLib::Date>> maturityDate,
            std::shared_ptr<std::vector<double>> optionPrice,
            std::shared_ptr<QuantLib::EndCriteria> endcriteria
            );
		//@}

		//! \name Public interface
		//{@
		//! Return the black volatility term structure
        QuantLib::BlackConstantVol getVol();
        double getPosJumpMean();
        double getNegJumpMean();
        double getPosProb();
        double getJumpIntensity();

		//!	Calibrate the Kou process
        void calibrate();
		//! Sum the all difference term given from values method
        QuantLib::Real value(const QuantLib::Array& x) const;
		//! Caculate the difference between the model value and market value
        QuantLib::Disposable<QuantLib::Array> values(const QuantLib::Array& x) const;
		//@}
    private:
        std::shared_ptr<std::vector<double>> optionPrice_;
        std::shared_ptr<std::vector<double>> strike_;
        std::shared_ptr<std::vector<QuantLib::Date>> maturityDate_;
        mutable QuantLib::Date referenceDate_;
        QuantLib::Calendar calendar_;
        QuantLib::DayCounter dayCounter_;

        double riskFreeRate_;
        double spotPrice_;
        double dividend_;

        std::shared_ptr<QuantLib::EndCriteria> endcriteria_;

        double volImplied_;
        double posJumpMeanImplied_;
        double negJumpMeanImplied_;
        double posProbability_;
        double jumpIntensity_;

    };
}
