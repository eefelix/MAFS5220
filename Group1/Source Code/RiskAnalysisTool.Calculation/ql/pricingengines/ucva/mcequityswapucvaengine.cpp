#include <pch.h>
#include "mcequityswapucvaengine.hpp"
#include <ql/pricingengines/swap/analyticequityswapengine.hpp>
#include <ql/pricingengines/ucva/ucvapathpricer.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>

//! \file mcequityswapucvaengine.cpp
//!
using namespace std;

namespace QuantLib {

	//! \brief UCVA Path Pricer for equity swap
	class EquitySwapUCVAPathPricer : public UCVAPathPricer<MCEquitySwapUCVAEngine> {
	public:

		//! \name Constructors & Destructors
		//@{
		EquitySwapUCVAPathPricer(const boost::shared_ptr<const MCEquitySwapUCVAEngine> &engine)
			: UCVAPathPricer(engine) {
		}

		~EquitySwapUCVAPathPricer() {
		}
		//@}

	protected:

		//! \brief Calculate the NPV of the remaining cash flow at default time \f$ \tau \f$.
		virtual Real defaultNPV(const MultiPath& path, Time defaultTime) const {
			TimeGrid timeGrid = path[0].timeGrid();
			std::vector<Time>::const_iterator defaultPos = std::find(timeGrid.begin(), timeGrid.end(), defaultTime);
			double spotprice = path[0][defaultPos - timeGrid.begin()];



			boost::shared_ptr<const MCEquitySwapUCVAEngine> engine = engine_.lock();
			EquitySwap::arguments* arguments = dynamic_cast<EquitySwap::arguments*>(engine->getArguments());
			Handle<YieldTermStructure> discountCurve = engine->discountCurve();

			double interestrate = arguments->discountCurve->forwardRate(
				defaultTime, defaultTime + 0.000001, QuantLib::Continuous);

			double cumulativedividend = 0;
			double cumulativecoupon = 0;
			size_t defaultgrid = std::min<size_t>(timeGrid.closestIndex(defaultTime), timeGrid.size() - 2);
			for (std::vector<double>::reverse_iterator it = arguments->yearfraction.rbegin();
				it != arguments->yearfraction.rend(); ++it){
				if (*it < defaultTime){
					cumulativecoupon = (defaultTime - *it)
						*arguments->fixedRate*arguments->startPrice;
					for (size_t j = timeGrid.closestIndex(*it); j <= defaultgrid; ++j){
						cumulativedividend += path[0][j] * arguments->dividend*timeGrid.dt(j);
					}
					if (it == (arguments->yearfraction.rend() - 1)){
						cumulativecoupon += arguments->cumulativeCoupoon;
						cumulativedividend += arguments->cumulativeDividend;
					}
					break;
				}
			}

			// calculate default date
			Date today = discountCurve->referenceDate();
			Date defaultDate = today
				+ int(defaultTime / arguments->discountCurve->timeFromReference(arguments->schedule->endDate())
				*arguments->discountCurve->dayCounter().dayCount
				(arguments->startDate, arguments->schedule->endDate()));

			EquitySwap ers(arguments->type, arguments->startDate, defaultDate,
				arguments->startPrice, spotprice, 1,
				cumulativecoupon, cumulativedividend,
				arguments->fixedRate, arguments->dividend, arguments->maturity, arguments->sigma,
				boost::shared_ptr<QuantLib::YieldTermStructure>(
				new QuantLib::FlatForward(defaultDate, interestrate, arguments->discountCurve->dayCounter())));
			/* the risk free rate should be adapted to the time level tau*/
			ers.setPricingEngine(boost::shared_ptr<AnalyticESEngine>(new AnalyticESEngine()));
			// cva is the positive part of NPV
			return ers.NPV();
		}
	};

	MCEquitySwapUCVAEngine::MCEquitySwapUCVAEngine(
		Handle<YieldTermStructure> riskFreeTermStructure,
		const boost::shared_ptr<const EquitySwap> &swap,
		const boost::shared_ptr<const Counterparty> &issuer,
		const Matrix &corr,
		Size timeStepsPerYear, bool antitheticVariate, Size requiredSamples, Real requiredTolerance, Size maxSamples, BigNatural seed)
		: MCUCVAEngine(swap, issuer, corr, swap->getMaturity(), riskFreeTermStructure, Null<Size>(),
		/* TODO: use actual # of days according to Day Counter*/ 360,
		true, requiredSamples, requiredTolerance, maxSamples, seed) {
	}

	MCEquitySwapUCVAEngine::~MCEquitySwapUCVAEngine() {
	}

	boost::shared_ptr<MCEquitySwapUCVAEngine::path_pricer_type> MCEquitySwapUCVAEngine::pathPricer() const {
		return boost::make_shared<EquitySwapUCVAPathPricer>(shared_from_this());
	}

	std::vector<boost::shared_ptr<StochasticProcess1D>> MCEquitySwapUCVAEngine::instrumentProcess() const {

		process = boost::make_shared<BlackScholesMertonProcess>(
			Handle<Quote>(boost::make_shared<SimpleQuote>(instrument_->getSpotPrice())),
			Handle<YieldTermStructure>(boost::make_shared<FlatForward>(
				discountCurve_->referenceDate(), instrument_->getDividendYield(), discountCurve_->dayCounter())),
			discountCurve_,
			Handle<BlackVolTermStructure>(boost::make_shared<BlackConstantVol>(
				discountCurve_->referenceDate(), instrument_->getCalendar(), instrument_->getVolatility(), discountCurve_->dayCounter()))
			);

		return std::vector<boost::shared_ptr<StochasticProcess1D>>({process});
	}

}