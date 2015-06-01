#include <pch.h>
#include "simcrosscurrencyswapucvaengine.hpp"
//
//using namespace QuantLib;
//
//int CCSEngine::defaultNumbers = 0;
//
//CCSEngine::CCSEngine(
//	double exchangeRate,
//	const Date& npvDate,
//	double rd,
//	double rf,
//	boost::shared_ptr<QuantLib::YieldTermStructure>	discountCurve,
//	long nPath /*= 50000*/,
//	long nSteps /*= 1000*/)
//	:exchangeRate_(exchangeRate), npvDate_(npvDate), domesticRate_(rd), foreignRate_(rf),
//	discountCurve_(discountCurve), nPath_(nPath), nSteps_(nSteps)
//{
//	registerWith(discountCurve_);
//}
//
//CCSEngine::~CCSEngine()
//{
//}
//
//void CCSEngine::calculate() const
//{
//	QL_REQUIRE(discountCurve_->referenceDate() == npvDate_,
//		"discount reference date should be the npv date");
//	results_.value = 0.0;
//	arguments_.validate();
//
//	double payNPV = 0.0;
//	double receiveNPV = 0.0;
//
//	pathPricer(receiveNPV, payNPV, exchangeRate_, npvDate_);
//
//	results_.value = receiveNPV - payNPV;
//	results_.legNPV.push_back(payNPV);
//	results_.legNPV.push_back(receiveNPV);
//}
//
//void CCSEngine::pathPricer(
//	double & receiveNPV,
//	double & payNPV,
//	const double & fxRate,
//	const Date& startDate) const
//{	
//	std::vector<Date>::iterator
//		payDateItr = arguments_.payDates.begin(),
//		receiveDateItr = arguments_.receiveDates.begin();
//	std::vector<Real>::iterator
//		payCouponItr = arguments_.payCoupons.begin(),
//		receiveCouponItr = arguments_.receiveCoupons.begin();
//
//	while (*payDateItr < startDate)
//	{
//		++payDateItr;
//		++payCouponItr;
//	}
//	while (*receiveDateItr < startDate)
//	{
//		++receiveDateItr;
//		++receiveCouponItr;
//	}
//
//
//	if (arguments_.type == CrossCurrencySwap::Type::payDomestic){
//		while (payDateItr != arguments_.payDates.end()){
//			payNPV += (*payCouponItr)*discountCurve_->discount(*payDateItr, true);
//			++payCouponItr;
//			++payDateItr;
//		}
//
//		while (receiveDateItr != arguments_.receiveDates.end()) {
//			receiveNPV +=
//				(*receiveCouponItr)
//				*fxRate*std::exp((domesticRate_ - foreignRate_)*arguments_.receiveDayCount.yearFraction(startDate, *receiveDateItr))
//				*discountCurve_->discount(*receiveDateItr, true);
//			++receiveCouponItr;
//			++receiveDateItr;
//		}
//	}
//	else{
//		while (payDateItr != arguments_.payDates.end()){
//			payNPV +=
//				(*payCouponItr)*fxRate*std::exp((domesticRate_ - foreignRate_)
//				*arguments_.receiveDayCount.yearFraction(startDate, *receiveDateItr))
//				*discountCurve_->discount(*payDateItr, true);
//			++payCouponItr;
//			++payDateItr;
//		}
//
//		while (receiveDateItr != arguments_.receiveDates.end()) {
//			receiveNPV +=
//				(*receiveCouponItr)
//				*discountCurve_->discount(*receiveDateItr, true);
//			++receiveCouponItr;
//			++receiveDateItr;
//		}
//
//	}
//
//};