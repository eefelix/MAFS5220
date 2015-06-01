#include <vector>
#include <boost\smart_ptr\shared_ptr.hpp>
#include <ql\instrument.hpp>
#include <ql\pricingengine.hpp>
#include <ql\credit\counterparty.hpp>
#include <ql\pricingengines\credit\mcexposuremodel.hpp>
#include <ql\models\shortrate\onefactormodels\coxingersollross.hpp>
#include <ql\stochasticprocess.hpp>
#include <map>

struct Product {
	std::vector<boost::shared_ptr<QuantLib::Instrument>> underlyings_;
	std::vector<boost::shared_ptr<const QuantLib::Counterparty>> counterparties_;
	std::vector<boost::shared_ptr<const QuantLib::MCExposureModel>> models_;
	boost::shared_ptr<QuantLib::PricingEngine> engine_;
	boost::shared_ptr<QuantLib::Instrument>	portfolioManager_;
	boost::shared_ptr<QuantLib::YieldTermStructure> domesticCurve_;
	boost::shared_ptr<QuantLib::CoxIngersollRoss> domesticCIR_;
	static std::map<std::string, boost::shared_ptr<QuantLib::StochasticProcess1D>> symToProcess;
};

