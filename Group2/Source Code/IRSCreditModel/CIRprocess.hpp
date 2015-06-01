#ifndef CIRprocess_h
#define CIRprocess_h

#include <ql\types.hpp>
#include <ql\time\date.hpp>
using namespace QuantLib;

/**
* Class Name: CIRprocess.
* Usage:
*	1. Model for CIR process
*/
class CIRprocess {
public:
	CIRprocess(Rate r01 = 0.05,
		Real theta1 = 0.1,
		Real k1 = 0.1,
		Real sigma1 = 0.1);
	~CIRprocess(){};

	Real P(Date startDate, Date endDate, Rate rt) const;
	Real D(Date startDate, Date endDate) const;

protected:

	Real A(Date startDate, Date endDate) const;
	Real B(Date startDate, Date endDate) const;

private:
	Real theta;
	Real k;
	Real sigma;
	Rate r0;
	Date startDate;
	Date endDate;
	const DayCounter dc = Actual360();
};


#endif