#include <pch.h>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <ql/termstructures/volatility/equityfx/blackvariancecurve.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/termstructures/yield/flatforward.hpp>

//! \file at1pdefaultmodel.cpp

using namespace QuantLib;

QuantLib::Integer AT1Pmodel::betaT(QuantLib::Time T) const{
	int s = 0;
	QL_REQUIRE(T >= 0, "T is greater than zero!!!");
	for (int i = Cdsmaturities->size() - 2; i >= 0; i--) {
		if (Cdsmaturities->at(i) < T) {
			s = i;
			break;
		}
	}
	return s;
}

QuantLib::Real AT1Pmodel::SigmaT(QuantLib::Time T, QuantLib::Array sigma) const{
	return sigma[betaT(T)];
}

QuantLib::Real AT1Pmodel::SigmaInt(QuantLib::Time T, QuantLib::Array sigma) const{
	QuantLib::Real S = 0;
	int n = betaT(T);
	QL_REQUIRE(sigma.size() == (Cdsmaturities->size() - 1) && T >= 0, "Invalid Data Input!!!");
	if (T == 0);
	else{
		for (int i = 0; i < n; i++) {
			S = S + sigma[i] * sigma[i] * (Cdsmaturities->at(i + 1) - Cdsmaturities->at(i));
		}
		S = S + (T - Cdsmaturities->at(n))*sigma[n] * sigma[n];
	}
	return S;
}

QuantLib::Real AT1Pmodel::QtauGreaterThanT(QuantLib::Time T, QuantLib::Array sigma) const{
	QuantLib::CumulativeNormalDistribution normal;
	QuantLib::Real d1, d2;
	QuantLib::Real sigmaInt = SigmaInt(T, sigma);
	d1 = normal((-log(H_V) + (B - 0.5)*sigmaInt) / sqrt(sigmaInt));
	d2 = normal((log(H_V) + (B - 0.5)*sigmaInt) / sqrt(sigmaInt));
	return d1 - pow(H_V, 2 * B - 1)*d2;
}

AT1Pmodel::AT1Pmodel(
	boost::shared_ptr<const std::vector<QuantLib::Period>> periods,
	boost::shared_ptr<QuantLib::YieldTermStructure> discountcurve,
	const QuantLib::Real dividend_,
	const QuantLib::Real h_v,
	const QuantLib::Real b,
	const QuantLib::Real x0_,
	boost::shared_ptr<QuantLib::Array> volatility_) :
	CalibratedModel(1),
	Volatility(arguments_[0]), H_V(h_v), B(b), x0(x0_), dividend(dividend_),
	Periods(periods), DiscountCurve(discountcurve)
{
	
	Volatility = PiecewiseConstantParameter(
		std::vector<double>(Periods->size()-1, 0.2),
		NonhomogeneousBoundaryConstraint(Array(Periods->size(), 0), Array(Periods->size(), 1)));
	setParams(Array(Periods->size(), 0.2));

	Cdsmaturities = boost::shared_ptr<std::vector<QuantLib::Time>>(
		new std::vector<QuantLib::Time>(Periods->size() + 1, 0));
	CDSmaturityDates = boost::shared_ptr<std::vector<QuantLib::Date>>(
		new std::vector<QuantLib::Date>(Periods->size(), Date()));

	for (size_t i = 0; i < Periods->size(); ++i){
		CDSmaturityDates->at(i) = DiscountCurve->referenceDate() + Periods->at(i);
		Cdsmaturities->at(i + 1) = DiscountCurve->timeFromReference(CDSmaturityDates->at(i));
	}
}

boost::shared_ptr<QuantLib::BlackVolTermStructure> AT1Pmodel::getVolTermSturcture() const
{

	std::vector<double> variance;
	for (size_t i = 1; i <= Volatility.size(); ++i){
		variance.push_back(sqrt(SigmaInt(Cdsmaturities->at(i), Volatility.params()) / Cdsmaturities->at(i)));
	}
	VolTermStructure = boost::shared_ptr<QuantLib::BlackVarianceCurve>(
		new QuantLib::BlackVarianceCurve(
		DiscountCurve->referenceDate(), *CDSmaturityDates, variance, DiscountCurve->dayCounter(), true));
	return VolTermStructure;
}


QuantLib::Probability AT1Pmodel::defaultBarrier(QuantLib::Time t) const
{
	if (process_ == nullptr) throw("Calibrate First!!!");
	else
		return process_->x0()*process_->dividendYield()->discount(t) / process_->riskFreeRate()->discount(t)*
		exp(-B*SigmaInt(t, Volatility.params()))*H_V;
}

const QuantLib::Time AT1Pmodel::defaultTime(const Path &path) const {
	TimeGrid timeGrid = path.timeGrid();
	int i = 0;
	while (i < timeGrid.size()){
		if (path[i] <= defaultBarrier(timeGrid[i]))
			return timeGrid[i];
		++i;
	}
	return Null<Time>();
}

boost::shared_ptr<StochasticProcess1D> AT1Pmodel::process() const
{
	getVolTermSturcture();
	process_ = boost::shared_ptr<BlackScholesMertonProcess>(new BlackScholesMertonProcess(
		Handle<Quote>(boost::make_shared<SimpleQuote>(x0)),
		Handle<YieldTermStructure>(boost::make_shared<FlatForward>(DiscountCurve->referenceDate(), dividend, DiscountCurve->dayCounter())),
		Handle<YieldTermStructure>(DiscountCurve),
		Handle<BlackVolTermStructure>(VolTermStructure)
		));

	return process_;
}

const QuantLib::Probability AT1Pmodel::defaultProbability(const QuantLib::Time T) const
{
	return (1 - QtauGreaterThanT(T, Volatility.params()));
}

void AT1Pmodel::calibrate(
	const std::vector<boost::shared_ptr<CalibrationHelper> >& helpers,
	OptimizationMethod& method,
	const EndCriteria& endCriteria,
	const Constraint& constraint,
	const std::vector<Real>& weights /*= std::vector<Real>()*/,
	const std::vector<bool>& fixParameters /*= std::vector<bool>()*/)
{
	CalibratedModel::calibrate(helpers, method, endCriteria, 
		NonhomogeneousBoundaryConstraint(Array(helpers.size(), 0), Array(helpers.size(), 1)),
		weights, fixParameters);
}