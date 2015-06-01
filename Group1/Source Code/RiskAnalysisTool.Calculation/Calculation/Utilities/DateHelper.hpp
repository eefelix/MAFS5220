#pragma once
#include <ql/types.hpp>
#include <ql/math/matrix.hpp>
#include <ql/time/date.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/timeunit.hpp>

namespace RiskAnalysisTool {
	namespace Calculation {
		namespace Utilities {
#ifdef _M_CEE
			public ref class DateHelper abstract sealed {
			public:
				static QuantLib::Date ToQLDate(System::DateTime date) {
					return QuantLib::Date(
						(QuantLib::Day)date.Day,
						(QuantLib::Month)date.Month,
						(QuantLib::Year)date.Year);
				}

				static QuantLib::Period ToQLPeriod(RiskAnalysisTool::Time::Period^ period) {
					QuantLib::TimeUnit unit;
					switch (period->Unit) {
					case RiskAnalysisTool::Time::PeriodUnit::Day:
						unit = QuantLib::Days;
						break;
					case RiskAnalysisTool::Time::PeriodUnit::Month:
						unit = QuantLib::Months;
						break;
					case RiskAnalysisTool::Time::PeriodUnit::Year:
						unit = QuantLib::Years;
						break;
					default:
						throw gcnew System::ArgumentException("Invalid period unit.", "period");
					}

					return QuantLib::Period(period->Value, unit);
				}
			};

			public ref class CorrelationHelper abstract sealed{
			public:
				static ::QuantLib::Matrix ToQLMatrix(cli::array<double, 2>^ m){
					::QuantLib::Matrix correlationMatrix(m->GetLength(0), m->GetLength(1));
					for (int i = 0; i < correlationMatrix.rows(); ++i)
						for (int j = 0; j < correlationMatrix.columns(); ++j)
							correlationMatrix[i][j] = m[i, j];

					return correlationMatrix;
				}
			};

			public ref class EndTimeHelper abstract sealed{
			public:
				static double GetEndTime(const QuantLib::DayCounter& dc, RiskAnalysisTool::Requests::ComputationRequest ^request) {
					std::vector<double> endTime;
					for each(RiskAnalysisTool::Instruments::Instrument^ instrument in request->Portfolio) {
						endTime.push_back(dc.yearFraction(
							Utilities::DateHelper::ToQLDate(instrument->EvaluationData), Utilities::DateHelper::ToQLDate(instrument->MaturityDate)));
					}

					return *max_element(endTime.begin(), endTime.end());
				}
			};

			public ref class StringHelper abstract sealed{
			public:
				static std::string ToNativeString(System::String^ str){
					array<System::Byte>^ bytes = System::Text::Encoding::UTF8->GetBytes(str);
					pin_ptr<System::Byte> ptr = &bytes[0];
					return std::string((char *)ptr);
				}
			};
#endif
		}
	}
}