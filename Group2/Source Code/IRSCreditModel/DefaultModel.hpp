#ifndef defaultmodel_hpp
#define defaultmodel_hpp

#include <ql\time\date.hpp>
#include <ql\types.hpp>

using namespace QuantLib;


/**
* Class Name: DefaultModel.
* Usage:
*	1. Abstract class for default model
*	2. Derived class for AT1P model and Intensity model
*/
class DefaultModel {
public :
	DefaultModel(){}
	~DefaultModel(){}
	virtual Probability SurvivalProbability(Date date) = 0;
};

#endif