#include <pch.h>

#include "SampleClass.hpp"

using namespace std;
using namespace RiskAnalysisTool::Calculation::Sample;

SampleClass::SampleClass() {
}

SampleClass::~SampleClass() {
}

double SampleClass::Sum(const shared_ptr<const vector<double>> &data) {
    double result = 0;
    for (auto it = data->begin(); it != data->end(); ++it) {
        result += *it;
    }
    return result;
}