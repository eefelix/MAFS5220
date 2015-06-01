#pragma once 


namespace RiskAnalysisTool {
	namespace Calculation {
		namespace Utilities {
#ifdef _M_CEE
			generic <typename T> 
				where T : ref class
			public ref class DynamicVisitor {
			public:
				void Visit(T target);
			};
#endif
		}
	}
}

