#pragma once

#include <Calculation/Calculation.h>
#include <ql/pricingengine.hpp>
//! \file creditvarengine.hpp
//! \brief base class for credit VaR engine. 
//!
//! This class is base class for credit VaR calculation, which is used to set some basic information.
//! Each derived class should implement its own method to calculate credit VaR.
//! 
//!	\date 2015-04-15
namespace QuantLib {
	//! \brief base class for credit VaR engine. 
	//!
	//! This class is base class for credit VaR calculation, which is used to set some basic information.
	//! Each derived class should implement its own method to calculate credit VaR.
    template <typename INST>
    class CreditVaREngine : public INST::engine {
    public:
		//! \name Constructors & Destructors
		//{@
        CreditVaREngine(
            const boost::shared_ptr<const INST> instrument, const boost::shared_ptr<const Counterparty> issuer,
            const Time &endTime, const Handle<YieldTermStructure> &discountCurve
            ) : instrument_(instrument), issuer_(issuer), discountCurve_(discountCurve), endTime_(endTime) {
        }

        virtual ~CreditVaREngine() {
        }
		//@}
    public:
        Time endTime() const {
            return endTime_;
        }
    protected:
        const Time endTime_;
        const boost::shared_ptr<const INST> instrument_;
        const boost::shared_ptr<const Counterparty> issuer_;
        const Handle<YieldTermStructure> discountCurve_;
    };
}