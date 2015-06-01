#pragma once

#include <vector>

// Structure to store upper and lower matrices
typedef struct{
	std::vector<std::vector<double> > L;
	std::vector<std::vector<double> > U;
}LUstruct;

/**
* Class Name: CMatrixOperations.
* Usage:
*	1. Utility tools used in American Monte Carlo
*/
class CMatrixOperations {
public:
	CMatrixOperations(void);
	~CMatrixOperations(void);

	// LU decomposition;
	LUstruct LU(std::vector<std::vector<double> > A);

	// Inverse of an upper triangular matrix
	std::vector<std::vector<double> >  MatUpTriangleInv(std::vector<std::vector<double> > U);

	// Inverse of a lower triangular matrix
	std::vector<std::vector<double> >  MatLowTriangleInv(std::vector<std::vector<double> > L);

	// Multiply two matrices together
	// First matrix is (n x k), Second is (k x m)
	// Resulting matrix is (n x m)
	std::vector<std::vector<double> > MMult(std::vector<std::vector<double> >A, std::vector<std::vector<double> >B, int n, int k, int m);

	// Inverse of a matrix through LU decomposition
	std::vector<std::vector<double> >  MInvLU(std::vector<std::vector<double> > A);

	// Transpose of a matrix
	std::vector<std::vector<double> > MTrans(std::vector<std::vector<double> >A);

	// Multiply a matrix by a vector
	std::vector<double> MVecMult(std::vector<std::vector<double> > X, std::vector<double> Y);

	// Regression beta coefficients
	std::vector<double> betavec(std::vector<std::vector<double> > X, std::vector<double> Y);

};

