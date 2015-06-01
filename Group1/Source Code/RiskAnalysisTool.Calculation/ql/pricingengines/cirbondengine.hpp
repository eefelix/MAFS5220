#pragma once

#include <Calculation/Calculation.h>
#include <ql/instruments/bonds/zerocouponbond.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>

//! \file cirbondengine.hpp
//! \brief Bond pricing engine based on CIR model
//!

namespace QuantLib {
	//! \brief A CIRBondEngine is a bond pricing engine based on CIR model.
	/*!	This class is to calculate bond prices using CIR model, which are used by the calibration helpers
		of CoxIngersollRoss model.
		Internally, this class calculate the bond price using
		\f[ e^{b(t)-a(t)r_t} \f]

		\see CoxIngersollRoss
		\see ZerocouponbondHelper
	*/
	//! 
	class _RISKANALYSISTOOL_CALCULATION_API CIRBondEngine
		:public GenericModelEngine<CoxIngersollRoss, ZeroCouponBond::arguments, ZeroCouponBond::results>
	{
	public:
		//! \name Constructors & Destructors
		//@{
		//! Creates an instance of CIRBondEngine using a CIR model instance and a yield term structure.
		CIRBondEngine(const boost::shared_ptr<CoxIngersollRoss>& model,
			const Handle<YieldTermStructure>& termStructure);
		//@}
		//! Calculates the model value of bond price using CIR model
		void calculate() const;
	private:
		Handle<YieldTermStructure> termStructure_;
	};
}