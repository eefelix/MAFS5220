#include <pch.h>
#include "mcvanillaswapucvaengine.hpp"
#include "ucvapathpricer.hpp"
#include <ql/pricingengines/blackformula.hpp>

using namespace std;

//! \file mcvanillaswapucvaengine.cpp
//!
namespace QuantLib {
	//! \brief UCVA Path Pricer for vanilla swap 
	//!
	class VanillaSwapUCVAPathPricer : public UCVAPathPricer < MCVanillaSwapUCVAEngine >{
	public:
		//! \name Constructors & Destructors
		//@{
		VanillaSwapUCVAPathPricer(const boost::shared_ptr<const MCVanillaSwapUCVAEngine> &engine)
			: UCVAPathPricer(engine){

		}
		//@}
	protected:
		//! \brief Calculate the NPV of the remaining cash flow at default time \f$ \tau \f$.
		//!
		//! The calculation method used here is based on CIR model.
		virtual Real defaultNPV(const MultiPath& path, Time defaultTime) const {
			//find short rate at default time
			TimeGrid timeGrid = path[0].timeGrid();

			auto defaultPos = timeGrid.index(defaultTime);
			Real shortRate = path[0][defaultPos];

			//calculate discount factor from reference time to the default time
			Time previousTime = timeGrid[0];
			Real discountFactor = 1.0;
			for (auto i = 0; i <= defaultPos; ++i){
				discountFactor *= std::exp(-path[0].value(i) * (timeGrid[i] - previousTime));
				previousTime = timeGrid[i];
			}

			boost::shared_ptr<const MCVanillaSwapUCVAEngine> engine = engine_.lock();
			VanillaSwap::arguments* arguments = dynamic_cast<VanillaSwap::arguments*>(engine->getArguments());
			Handle<YieldTermStructure> discountCurve = engine->discountCurve();

			Date referenceDate = discountCurve->referenceDate();
			std::vector<Date> fixedDates = arguments->fixedPayDates;
			std::vector<Date> floatingDates = arguments->floatingPayDates;
			// transfer payment Dates to Times
			std::vector<Real> fixedTimes;
			std::vector<Real> floatingTimes;
			for (int i = 0; i != fixedDates.size(); ++i){
				fixedTimes.push_back(engine->instrument()->fixedDayCount().yearFraction(referenceDate, fixedDates[i]));
			}
			for (int i = 0; i != floatingDates.size(); ++i){
				floatingTimes.push_back(engine->instrument()->floatingDayCount().yearFraction(referenceDate, floatingDates[i]));
			}
			// calculate the corresponding discount factors of each leg to default time
			std::vector<Real> fixedDiscountFactor;
			std::vector<Real> floatingDiscountFactor;
			int fixedIndex = 0;
			int floatingIndex = 0;
			for (int i = 0; i != fixedTimes.size(); ++i){
				if (fixedTimes[i] >= defaultTime){
					fixedDiscountFactor.push_back(engine->Model()->discountBond(0.0, fixedTimes[i] - defaultTime, shortRate));
				}
				else{
					++fixedIndex;
				}
			}
			for (int i = 0; i != floatingTimes.size(); ++i){
				if (floatingTimes[i] >= defaultTime){
					floatingDiscountFactor.push_back(engine->Model()->discountBond(0.0, floatingTimes[i] - defaultTime, shortRate));
				}
				else{
					++floatingIndex;
				}
			}

			//calculate the NPV(tau)
			Real NPVtau = 0.0;
			NPVtau += arguments->nominal * (1. - floatingDiscountFactor.back());
			for (int m = 0; m != floatingDiscountFactor.size(); ++m){
				NPVtau += arguments->nominal * arguments->floatingSpreads[m + floatingIndex] * floatingDiscountFactor[m];
			}
			for (int m = 0; m != fixedDiscountFactor.size(); ++m){
				NPVtau -= arguments->fixedCoupons[m + fixedIndex] * fixedDiscountFactor[m];
			}

			if (arguments->type == QuantLib::VanillaSwap::Payer){
				NPVtau = NPVtau;
			}
			else{
				NPVtau = -NPVtau;
			}
			return NPVtau * discountFactor;
		}
	};

	MCVanillaSwapUCVAEngine::MCVanillaSwapUCVAEngine(
		const Calendar &calender, const DayCounter &dayCounter, Date referenceDate,
		const Handle<YieldTermStructure> &riskFreeTermStructure,
		const boost::shared_ptr<const VanillaSwap> &swap,
		const boost::shared_ptr<const Counterparty> &issuer,
		const boost::shared_ptr<const CoxIngersollRoss> &cirModel,
		const Matrix &corr,
		Size timeStepsPerYear/* = 360*/, bool antitheticVariate/* = true*/,
		Size requiredSamples/* = 50000*/, Real requiredTolerance/* = 0.0001*/, Size maxSamples/* = QL_MAX_INTEGER*/,
		BigNatural seed/* = SeedGenerator::instance().get()*/)
		: MCUCVAEngine(swap, issuer, corr,
		dayCounter.yearFraction(swap->startDate(), swap->maturityDate(), referenceDate, referenceDate),
		riskFreeTermStructure, Null<Size>(), 360, true, requiredSamples, requiredTolerance, maxSamples, seed),
		cirModel_(cirModel)
	{
	}

	MCVanillaSwapUCVAEngine::~MCVanillaSwapUCVAEngine() {
	}

	boost::shared_ptr<MCVanillaSwapUCVAEngine::path_pricer_type> MCVanillaSwapUCVAEngine::pathPricer() const {
		return boost::make_shared<VanillaSwapUCVAPathPricer>(shared_from_this());
	}
}