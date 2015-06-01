#include <pch.h>
#include "DynamicVisitor.hpp"

using namespace System;
using namespace System::Linq;
using namespace System::Collections::Generic;
using namespace System::Reflection;
using namespace RiskAnalysisTool::Calculation::Utilities;



generic <typename T>
void DynamicVisitor<T>::Visit(T target) {
	Type ^visitorType = this->GetType();

	visitorType->InvokeMember(
		"OnVisit",
		BindingFlags::InvokeMethod | BindingFlags::NonPublic | BindingFlags::Public | BindingFlags::Instance,
		nullptr,
		this,
		gcnew array < Object ^ > {target});
}