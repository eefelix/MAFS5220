#include <pch.h>
#include <ql/pricingengines/credit/mcexposuremodel.hpp>
#include <ql/instruments/credit/tva.hpp>
#include <ql/instruments/credit/bva.hpp>
#include <ql/models/defaultmodel.hpp>
#include <ql/models/default/at1pdefaultmodel.hpp>
#include <ql/credit/counterparty.hpp>
#include <ql/pricingengines/credit/mccrosscurrencyswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mcvanillaswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mcequityswapexposuremodel.hpp>
#include <ql/pricingengines/credit/mctvaengine.hpp>
#include <ql/pricingengines/credit/mcbvaengine.hpp>
#include <ql/models/shortratetermstructure.hpp>
#include <ql/instruments/credit/var.hpp>
#include <ql/pricingengines/credit/mccreditvarengine.hpp>
#include "Utilities/DateHelper.hpp"
#include "InstrumentFactory.hpp"
#include "ModelFactory.hpp"
#include "PortfolioManagerFactory.hpp"
#include "EngineFactory.hpp"
#include "CounterpartyFactory.hpp"
#include <ql/termstructures/yield/discountcurve.hpp>
#include <boost/make_shared.hpp>
#include <utility>

#include <Calculation/Utilities/cliext.hpp>
#include <ql/models/shortrate/onefactormodels/coxingersollross.hpp>
#include <ql/models/calibrationhelper.hpp>
#include <ql/pricingengines/cirbondengine.hpp>
#include <ql/models/shortrate/calibrationhelpers/zerocouponbondhelper.hpp>
#include <ql/math/optimization/simplex.hpp>
#include "Product.hpp"
#include <ql/instruments/vanillaswap.hpp>

// This is a managed class that will be exposed to WebAPI, we don't need .h/.hpp file because managed classes are exposed via metadata.
using namespace System;
using namespace System::Linq;
using namespace System::Collections::Generic;
using namespace RiskAnalysisTool::Instruments;
using namespace RiskAnalysisTool::Requests;
using namespace RiskAnalysisTool::Results;
using namespace std;

namespace RiskAnalysisTool {
	namespace Calculation {
		public ref class ComputeEngine sealed {
			// Make this class stateless, thread-safe and re-entrance, so NO non-static read/write member variables
			ref class Builder;

		public:
			ComputeEngine() {
			}

			~ComputeEngine() {
				this->!ComputeEngine();
			}

			!ComputeEngine() {
			}

			ComputationResult^ calculate(ComputationRequest^ request);
		};

		

		ref class ComputeEngine::Builder {
		public:
			Builder() {
				product_ = new Product();
			}

			~Builder() {
				this->!Builder();
			}

			!Builder() {
				if (!product_)
					delete product_;
			}

			Product* getProduct(ComputationRequest^ request) {
				createDiscountCurve(request);
				createInstrument(request);
				createCounterparties(request);
				createEngine(request);
				createPortfolioManager(request);
				Product::symToProcess.clear();
				//product_->portfolioManager_->setPricingEngine(product_->engine_);

				Product* temp = product_;
				product_ = 0;
				return temp;
			}
		protected:
			void createCounterparties(ComputationRequest^ request) {
				CounterpartyFactory^ counterpartyFactory = gcnew CounterpartyFactory(&product_->domesticCurve_);
				boost::shared_ptr<std::vector<boost::shared_ptr<QuantLib::Counterparty>>> cPtr(counterpartyFactory->CreateCounterparties(request));

				product_->counterparties_.push_back(cPtr->at(0));
				product_->counterparties_.push_back(cPtr->at(1));
				Product::symToProcess.emplace(std::make_pair("ISSUER", cPtr->at(0)->createDefaultProcess()));
				Product::symToProcess.emplace(std::make_pair("INVESTOR", cPtr->at(1)->createDefaultProcess()));
			}

			void createEngine(ComputationRequest^ request) {
				EngineFactory^ eFactory = gcnew EngineFactory(&product_->domesticCurve_, &product_->counterparties_, &product_->models_, &product_->domesticCIR_);
				product_->engine_ = eFactory->CreateEngine(request);
			}

			void createPortfolioManager(ComputationRequest^ request) {
				PortfolioManagerFactory^ pmFactory = gcnew PortfolioManagerFactory;
				product_->portfolioManager_.reset(pmFactory->CreatePortfolioManager(request));
			}
			void createInstrument(ComputationRequest^ request) {
				InstrumentFactory^ instFactory = gcnew InstrumentFactory(request, &product_->domesticCurve_);
				for each(Instrument^ instrument in request->Portfolio) {
					boost::shared_ptr<QuantLib::Instrument> instPtr(instFactory->CreateQLInstrument(instrument));
					product_->underlyings_.push_back(instPtr);

					ModelFactory^ modelFactory = gcnew ModelFactory(request, &product_->domesticCurve_, &instPtr);
					boost::shared_ptr<QuantLib::MCExposureModel> modelPtr(modelFactory->CreateQLModel(instrument));
					product_->models_.push_back(modelPtr);
				}
			}

			void createDiscountCurve(ComputationRequest^ request) {
				QuantLib::Date referenceDate = Utilities::DateHelper::ToQLDate(request->EvaluationData);
				String^ domesticCurrency = request->DomesticCurrency;

				vector<QuantLib::Date> domesticBondDate;
				vector<double> domesticDiscountor;

				for each (auto bond in Enumerable::OfType<Bond ^>(request->MarketData)) {
					String^ symbol = bond->Symbol;
					if (symbol->Contains(domesticCurrency)){
						domesticBondDate.push_back(Utilities::DateHelper::ToQLDate(bond->MaturityDate));
						domesticDiscountor.push_back(bond->Price / bond->Principle);
					}
				}

				if (std::find(domesticBondDate.begin(), domesticBondDate.end(), referenceDate) == domesticBondDate.end()){
					domesticBondDate.push_back(referenceDate);
					domesticDiscountor.push_back(1.);
				}

				product_->domesticCurve_
					= boost::make_shared<QuantLib::DiscountCurve>(domesticBondDate, domesticDiscountor, QuantLib::Actual360(), QuantLib::TARGET());
				product_->domesticCIR_ = boost::make_shared<QuantLib::CoxIngersollRoss>();
				
				std::vector<boost::shared_ptr<QuantLib::CalibrationHelper>> zcbonds;
				boost::shared_ptr<QuantLib::PricingEngine> engine(
					new QuantLib::CIRBondEngine(product_->domesticCIR_, QuantLib::Handle<QuantLib::YieldTermStructure>(product_->domesticCurve_)));

				std::vector<QuantLib::Date> maxDate;
				for each(RiskAnalysisTool::Instruments::Instrument^ instrument in request->Portfolio) {
					maxDate.push_back(Utilities::DateHelper::ToQLDate(instrument->MaturityDate));
				}

				QuantLib::Date mDate = *max_element(maxDate.begin(), maxDate.end());

				for (int i = 0; i < domesticBondDate.size(); ++i){
					if (domesticBondDate[i] != referenceDate && domesticBondDate[i]<=mDate){
						zcbonds.push_back(boost::shared_ptr<QuantLib::ZerocouponbondHelper>(new QuantLib::ZerocouponbondHelper(
							0,
							product_->domesticCurve_->calendar(),
							1,
							domesticBondDate[i],
							QuantLib::Following,
							100.0,
							referenceDate,
							QuantLib::Handle<QuantLib::Quote>(boost::make_shared<QuantLib::SimpleQuote>(0.1)),
							QuantLib::Handle<QuantLib::YieldTermStructure>(product_->domesticCurve_))));
						zcbonds.back()->setPricingEngine(engine);
					}
				}
				QuantLib::Simplex solver(0.001);
				const QuantLib::Size maxIteration = 1000;
				const QuantLib::Size minStatIteration = 50;
				const QuantLib::Real rootEpsilon = 1e-8;
				const QuantLib::Real FunctionEpsilon = 1e-8;
				const QuantLib::Real gradientNormEpsilon = 1e-8;
				const QuantLib::EndCriteria endcriteria = QuantLib::EndCriteria(maxIteration, minStatIteration, rootEpsilon, FunctionEpsilon, gradientNormEpsilon);

				product_->domesticCIR_->calibrate(zcbonds, solver, endcriteria, *(product_->domesticCIR_->constraint()));

				Product::symToProcess.emplace(std::make_pair(
					Utilities::StringHelper::ToNativeString(request->DomesticCurrency),
					boost::make_shared<QuantLib::CIRprocess>(
					product_->domesticCIR_->params()[0]
					, product_->domesticCIR_->params()[1]
					, product_->domesticCIR_->params()[2]
					, product_->domesticCIR_->params()[3])));
			}

		protected:
			Product* product_;
		};

		ComputationResult^ ComputeEngine::calculate(ComputationRequest^ request) {
			QuantLib::Settings::instance().evaluationDate() = Utilities::DateHelper::ToQLDate(request->EvaluationData);

			Type^ requestType = request->GetType();
			Builder^ builder = gcnew Builder();
			boost::shared_ptr<Product> product;
			if (requestType == TvaRequest::typeid) {
				TvaResult^ result = gcnew TvaResult();
				product.reset(builder->getProduct(request));
				assert(product != 0);
				boost::shared_ptr<QuantLib::TVA> tvaManager = boost::dynamic_pointer_cast<QuantLib::TVA>(product->portfolioManager_);
				assert(tvaManager != 0);
				tvaManager->setPricingEngine(product->engine_);
				result->Fva = tvaManager->FVA();

				return result;
			}
			else if (requestType == BvaRequest::typeid) {
				BvaResult^ result = gcnew BvaResult();
				product.reset(builder->getProduct(request));
				boost::shared_ptr<QuantLib::BVA> bvaManager = boost::dynamic_pointer_cast<QuantLib::BVA>(product->portfolioManager_);
				bvaManager->setPricingEngine(product->engine_);
				result->Cva = bvaManager->CVA();
				result->Dva = bvaManager->DVA();

				return result;
			}
			else if (requestType == UCvaRequest::typeid) {
				UCvaResult^ result = gcnew UCvaResult();
				product.reset(builder->getProduct(request));
				boost::shared_ptr<QuantLib::BVA> ucvaManager = boost::dynamic_pointer_cast<QuantLib::BVA>(product->portfolioManager_);
				ucvaManager->setPricingEngine(product->engine_);
				result->UCva = ucvaManager->CVA();

				return result;
			}
			else if (requestType == CreditVaRRequest::typeid) {
				CreditVaRResult^ result = gcnew CreditVaRResult();
				product.reset(builder->getProduct(request));
				boost::shared_ptr<QuantLib::VAR> varManager = boost::dynamic_pointer_cast<QuantLib::VAR>(product->portfolioManager_);
				varManager->setPricingEngine(product->engine_);
				result->CreditVaR = varManager->CreditVAR();

				return result;
			}
			else {
				throw gcnew NotSupportedException(String::Format("Unsupported Request Type."));
			}
		}
	}
}