#pragma once

#include <Calculation/Calculation.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <ql/types.hpp>
#include <ql/instrument.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/models/shortratetermstructure.hpp>

//! \file creditvarengine.hpp
//! \brief base class for credit var engine.

namespace QuantLib {
	template<typename INST>
	class CreditVaREngine : public INST::engine {
	public:
		CreditVaREngine(
			Time endTime,
			const std::vector<boost::shared_ptr<const Counterparty>> &counterparties,
			const std::vector<boost::shared_ptr<const MCExposureModel>> &creditvarModels,
			const Handle<YieldTermStructure> &domesticCurve, const boost::shared_ptr<OneFactorAffineModel>& domesticCIRModel
			) : endTime_(endTime), counterparties_(counterparties), CreditVarModels_(creditvarModels), domesticCurve_(domesticCurve), domesticCIRModel_(domesticCIRModel)
		{
			domesticShortRateModel_ = boost::make_shared<ShortRateTermStructure>(domesticCurve_->calendar(), domesticCurve_->dayCounter(), domesticCIRModel_);
		}

		virtual ~CreditVaREngine() {
		}

	public:
		std::vector<boost::shared_ptr<const Counterparty>> counterparties() const {
			return counterparties_;
		}

		/*std::vector<boost::shared_ptr<const Instrument>> portfolio() {
		return portfolio_;
		}*/

		const Time endTime() const { return endTime_; }
		const Handle<YieldTermStructure> & domesticCurve() const {
			return domesticCurve_;
		}

		const boost::shared_ptr<ShortRateTermStructure>& domesticShortRateModel() const{
			return domesticShortRateModel_;
		}

		const boost::shared_ptr<OneFactorAffineModel>& domesticCIRModel() const{
			return domesticCIRModel_;
		}

		const std::vector<boost::shared_ptr<const MCExposureModel>>& models() const{
			return CreditVarModels_;
		}

	protected:
		//std::vector<boost::shared_ptr<const Instrument>> portfolio_;
		std::vector<boost::shared_ptr<const Counterparty>> counterparties_;
		Handle<YieldTermStructure> domesticCurve_;
		boost::shared_ptr<OneFactorAffineModel> domesticCIRModel_;
		boost::shared_ptr<ShortRateTermStructure> domesticShortRateModel_;
		std::vector<boost::shared_ptr<const MCExposureModel>> CreditVarModels_;
		Time endTime_;
	};
}