#ifndef CIRcalibration_hpp
#define CIRcalibration_hpp

#include <vector>
#include <math.h>

using namespace QuantLib;

/**
* Class Name: CIRProblemFunction.
* Usage:
*	1. Derived from CostFunction class to do the optimization.
*/
class CIRProblemFunction : public CostFunction{
public:
	CIRProblemFunction()
	{
		dc = Actual360();
	}
	~CIRProblemFunction(){};

	void Initialization(std::vector<Real> quoteVector);
	Real value(const Array& x) const;
	Disposable<Array> values(const Array& x) const;

	/**
	* Function Name: calculation.
	* Return: <Real>
	* Usage:
	*	1. Get the residual terms using in the calibration.
	*/
	Real Calculation(Real rt_, Real rt_1, Real phi1, Real gamma1) const;
	Real getGamma() { return gamma; }

private:
	DayCounter dc;
	int termNum;
	Real phi;
	Real gamma;
	std::vector<Real> rt;
	std::vector<Real> quote;

};


/**
* Class Name: CIRcalibration.
* Usage:
*	1. Calculate the CIR parameters.
*/
class CIRcalibration{
public:
	CIRcalibration()
	{
		dc = Actual360();
	};
	~CIRcalibration(){};

	void Initialization(std::vector<Real> quoteVector);

	/**
	* Function Name: StartCalibration.
	* Return: Null
	* Usage:
	*	1. Use Optimizer to calibrate the optimal volatilities.
	*/
	void StartCalibrate();
	/**
	* Function Name: CalibrationParameters.
	* Return: Null
	* Usage:
	*	1. Calculate the parameters used in CIRprocess
	*/
	void CalculateParameters();
	/**
	* Parameter Name: optimizeFunc.
	* Type: <CIRProblemFunction>
	* Usage:
	*	1. Fitting CIR parameters
	*/
	CIRProblemFunction optimizeFunc;
	Real getFinalizedVolatilty() { return ContinousSigma; }
	Real getFinalizedDrift() { return k; }
	Real getFinalizedCenter() { return gamma; }

private:
	DayCounter dc;
	int termNum;
	std::vector<Real> quote;
	Real phi;
	Volatility DiscreteSigma;
	Volatility ContinousSigma;
	Real k;
	Real gamma;
	Real RSS;
	Real r0;
};

#endif