#include <math.h>
#include <ql/quantlib.hpp>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include"MatrixOperations.h"


using namespace QuantLib;
using namespace std;
using namespace boost::numeric::ublas;

/**
* Class Name: LeastSquareMC.
* Usage:
*	1. Implementation of American Monte Carlo, useful in FVA calculation
*/
class LeastSquareMC : public CMatrixOperations{
public:
	enum Type { Put = -1, Call = 1 };
	LeastSquareMC(){}
	LeastSquareMC(int M, int N, double initPrice, double T, double strike,
		double rate, double vol, Type type);

	/**
		Utility functions to determine the early exercise boundary
		Details refer to the explanations of cpp file
	*/

	double computeLSM(int time, matrix<double>& Paths, matrix<double>& CashFlow, matrix<double>& Exercise);
	double computeValue(matrix<double>& CashFlow);
	double max(double x, double y);
	/**
	* Function Name: calcLeastSquareMC
	* Usage: calculate American call option price
	*/
	double calcLeastSquareMC(); 

private:
	Type type_;
	long M_;  //M_ represents simulation path number
	long N_;
	double rate_;
	double vol_;
	double initPrice_;
	double T_;
	double dt_;
	double strike_;
	std::vector<double> prices[10000];
	//each element of array prices is a vector
	//number in bracket should at least be equal to M
};