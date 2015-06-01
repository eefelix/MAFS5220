#pragma once
#include <Calculation/Calculation.h>
#include <memory>
#include <vector>
namespace RiskAnalysisTool {
    namespace Calculation {
        namespace Sample {
            class _RISKANALYSISTOOL_CALCULATION_API SampleClass {
            public:
                SampleClass();
                virtual ~SampleClass();
            public:
                double Sum(const std::shared_ptr<const std::vector<double>> &data);
            };
        }
    }
}
