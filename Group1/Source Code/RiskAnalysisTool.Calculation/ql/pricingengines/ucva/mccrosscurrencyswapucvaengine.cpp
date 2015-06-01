#include <pch.h>
#include "mccrosscurrencyswapucvaengine.hpp"
#include "ucvapathpricer.hpp"
#include <ql/math/randomnumbers/boxmullergaussianrng.hpp>

using namespace std;

//! \file mccrosscurrencyswapucvaengine.cpp
//!
namespace QuantLib {
	//! \brief UCVA Path Pricer for fixed-fixed cross currency swap 
	//!
	class CrossCurrencySwapUCVAPathPricer : public UCVAPathPricer < MCCrossCurrencySwapUCVAEngine > {
	public:
		//! \name Constructors & Destructors
		//@{
		CrossCurrencySwapUCVAPathPricer(const boost::shared_ptr<const MCCrossCurrencySwapUCVAEngine> &engine)
			: UCVAPathPricer(engine) {
		}

		~CrossCurrencySwapUCVAPathPricer() {
		}
		//@}
	protected:
		//! \brief Calculate the NPV of the remaining cash flow at default time \f$ \tau \f$.
		virtual Real defaultNPV(const MultiPath& path, Time defaultTime) const {


			// find fx rate at default time
			TimeGrid timeGrid = path[0].timeGrid();
			std::vector<Time>::const_iterator defaultPos = std::find(timeGrid.begin(), timeGrid.end(), defaultTime);
			double domesticRate = path[0][defaultPos - timeGrid.begin()];
			double foreignRate = path[1][defaultPos - timeGrid.begin()];

			boost::shared_ptr<const MCCrossCurrencySwapUCVAEngine> engine = engine_.lock();
			CrossCurrencySwap::arguments* arguments = dynamic_cast<CrossCurrencySwap::arguments*>(engine->getArguments());
			Handle<YieldTermStructure> discountCurve = engine->discountCurve();

			double fxRate = arguments->fxRate;
			double fxVol = arguments->fxVol;
			fxRate *= std::exp(-0.5*fxVol*fxVol*(defaultTime));
			Time previousTime = timeGrid[0];

			for (auto itr = timeGrid.begin(); itr <= defaultPos; ++itr)
			{
				fxRate *= std::exp(path[0][itr - timeGrid.begin()] - path[1][itr - timeGrid.begin()] * (*itr - previousTime));
				previousTime = *itr;
			}

			static long seed = QuantLib::SeedGenerator::instance().get();
			static MersenneTwisterUniformRng rnd(seed);
			static BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> normrnd(rnd);
			fxRate *= std::exp(fxVol*normrnd.next().value*std::sqrt(defaultTime));

			// calculate default date
			Date today = discountCurve->referenceDate();
			Date defaultDate = today
				+ int(defaultTime / arguments->payDayCount.yearFraction(today, arguments->payDates.back())
				*arguments->payDayCount.dayCount(arguments->payDates.front(), arguments->payDates.back()));

			// calculate remaining pay and receive dates
			std::function<bool(const Date&)> dateCompare = [&defaultDate](const Date& d){ return d > defaultDate; };
			std::vector<Date>::iterator
				payDateItr = std::find_if(arguments->payDates.begin(), arguments->payDates.end(), dateCompare),
				receiveDateItr = std::find_if(arguments->receiveDates.begin(), arguments->receiveDates.end(), dateCompare);

			Real payNPV = 0., receiveNPV = 0.;
			if (arguments->type == CrossCurrencySwap::Type::payDomestic){
				while (payDateItr != arguments->payDates.end()){
					payNPV +=
						arguments->payCoupons[payDateItr - arguments->payDates.begin()] * discountCurve->discount(*payDateItr, true);
					++payDateItr;
				}

				while (receiveDateItr != arguments->receiveDates.end()) {
					receiveNPV +=
						arguments->receiveCoupons[receiveDateItr - arguments->receiveDates.begin()]
						* fxRate*std::exp((domesticRate - foreignRate)*arguments->receiveDayCount.yearFraction(defaultDate, *receiveDateItr))
						*discountCurve->discount(*receiveDateItr, true);
					++receiveDateItr;
				}
			}
			else{
				while (payDateItr != arguments->payDates.end()){
					payNPV +=
						arguments->payCoupons[payDateItr - arguments->payDates.begin()]
						* fxRate*std::exp((domesticRate - foreignRate)*arguments->payDayCount.yearFraction(defaultDate, *payDateItr))
						*discountCurve->discount(*payDateItr, true);
					++payDateItr;
				}

				while (receiveDateItr != arguments->receiveDates.end()) {
					receiveNPV +=
						arguments->receiveCoupons[receiveDateItr - arguments->receiveDates.begin()]
						* discountCurve->discount(*receiveDateItr, true);
					++receiveDateItr;
				}

			}

			return receiveNPV - payNPV;
		}
	};

	MCCrossCurrencySwapUCVAEngine::MCCrossCurrencySwapUCVAEngine(
		const Calendar &calender, const DayCounter &dayCounter, Date referenceDate,
		const Handle<YieldTermStructure> &riskFreeTermStructure,
		const boost::shared_ptr<const CrossCurrencySwap> &swap,
		const boost::shared_ptr<const Counterparty> &issuer,
		Rate domesticRate, Real domesticRateSpeed, Real domesticRateMean, Volatility domesticRateVol,
		Rate foreignRate, Real foreignRateSpeed, Real foreignRateMean, Volatility foreignRateVol, const Matrix &corr,
		Size timeStepsPerYear /*= 360*/, bool antitheticVariate /*= true*/,
		Size requiredSamples /*= 50000*/, Real requiredTolerance /*= 0.0001*/, Size maxSamples/* = QL_MAX_INTEGER*/,
		BigNatural seed /*= SeedGenerator::instance().get()*/)
		: MCUCVAEngine(swap, issuer, corr,
		dayCounter.yearFraction(swap->startDate(), swap->maturityDate(), referenceDate, referenceDate),
		riskFreeTermStructure, Null<Size>(), 360, true, requiredSamples, requiredTolerance, maxSamples, seed)
		, domesticRate_(domesticRate), domesticRateMean_(domesticRateMean), domesticRateSpeed_(domesticRateSpeed), domesticRateVol_(domesticRateVol)
		, foreignRate_(foreignRate), foreignRateMean_(foreignRateMean), foreignRateSpeed_(foreignRateSpeed), foreignRateVol_(foreignRateVol)
	{
	}

	MCCrossCurrencySwapUCVAEngine::~MCCrossCurrencySwapUCVAEngine() {
	}

	boost::shared_ptr<MCCrossCurrencySwapUCVAEngine::path_pricer_type> MCCrossCurrencySwapUCVAEngine::pathPricer() const {
		return boost::make_shared<CrossCurrencySwapUCVAPathPricer>(shared_from_this());
	}

	//std::vector<boost::shared_ptr<StochasticProcess1D>> MCCrossCurrencySwapBVAEngine::instrumentProcess() const {
	//    QL_FAIL("Not Implemented!");
	//}
}