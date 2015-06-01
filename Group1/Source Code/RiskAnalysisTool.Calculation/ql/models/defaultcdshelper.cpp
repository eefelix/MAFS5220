#include <pch.h>

#include <ql/models/defaultcdshelper.hpp>
#include <ql/models/defaultmodel.hpp>

using namespace QuantLib;

DefaultCdsHelper::DefaultCdsHelper(
	Time cdsTime,
	const Handle<Quote>& impliedDefaultProb,
	const Handle<YieldTermStructure>& termStructure,
	const boost::shared_ptr<CalibratedModel> & model
	) :CalibrationHelper(impliedDefaultProb, termStructure), cdsTime_(cdsTime),model_(model)
{
	registerWith(model_);
}

DefaultCdsHelper::~DefaultCdsHelper()
{
}

Real DefaultCdsHelper::modelValue() const 
{
	boost::shared_ptr<DefaultModel>	defaultModel = boost::dynamic_pointer_cast<DefaultModel>(model_);
	if (defaultModel){
		return defaultModel->defaultProbability(cdsTime_);
	}
	else
	{
		QL_FAIL("Cannot use default cds helper to calibrate non-default model.");
	}
}

Real DefaultCdsHelper::blackPrice(Volatility impliedDefaultProb) const 
{
	return impliedDefaultProb;
}

void DefaultCdsHelper::addTimesTo(std::list<Time>& times) const 
{
	QL_FAIL("Not implemented in DefaultCdsHelper.");
}