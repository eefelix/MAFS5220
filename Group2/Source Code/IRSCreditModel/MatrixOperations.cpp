// Matrix inversion, multiplication, LU decomposition, and other operations
// by Fabrice Douglas Rouah, FRouah.com and Volopta.com

#include "pch.h"
#include "MatrixOperations.h"

CMatrixOperations::CMatrixOperations(void) { }
CMatrixOperations::~CMatrixOperations(void) { }

// LU decomposition;
LUstruct CMatrixOperations::LU(std::vector<std::vector<double> > A) {
	int N = A.size();
	std::vector<std::vector<double> >  B(N, std::vector<double>(N));
	for (int i = 0; i <= N - 1; i++)
		for (int j = 0; j <= N - 1; j++)
			B[i][j] = A[i][j];

	for (int k = 0; k <= N - 2; k++) {
		for (int i = k + 1; i <= N - 1; i++)
			B[i][k] = B[i][k] / B[k][k];
		for (int j = k + 1; j <= N - 1; j++)
			for (int i = k + 1; i <= N - 1; i++)
				B[i][j] = B[i][j] - B[i][k] * B[k][j];
	}
	std::vector<std::vector<double> >  L(N, std::vector<double>(N));
	std::vector<std::vector<double> >  U(N, std::vector<double>(N));
	for (int i = 0; i <= N - 1; i++)	{
		L[i][i] = 1.0;
		for (int j = 0; j <= N - 1; j++) {
			if (i>j)
				L[i][j] = B[i][j];
			else
				U[i][j] = B[i][j];
		}
	}
	LUstruct Mats;
	Mats.L = L;
	Mats.U = U;
	return Mats;
}

// Inverse of an upper triangular matrix
std::vector<std::vector<double> >  CMatrixOperations::MatUpTriangleInv(std::vector<std::vector<double> > U) {
	int N = U.size();
	std::vector<std::vector<double> >  V(N, std::vector<double>(N));
	for (int j = N - 1; j >= 0; j--) {
		V[j][j] = 1.0 / U[j][j];
		for (int i = j - 1; i >= 0; i--)
			for (int k = i + 1; k <= j; k++)
				V[i][j] -= 1.0 / U[i][i] * U[i][k] * V[k][j];
	}
	return V;
}

// Inverse of a lower triangular matrix
std::vector<std::vector<double> >  CMatrixOperations::MatLowTriangleInv(std::vector<std::vector<double> > L) {
	int N = L.size();
	std::vector<std::vector<double> >  V(N, std::vector<double>(N));
	for (int i = 0; i <= N - 1; i++)	{
		V[i][i] = 1.0 / L[i][i];
		for (int j = i - 1; j >= 0; j--)
			for (int k = i - 1; k >= j; k--)
				V[i][j] -= 1.0 / L[i][i] * L[i][k] * V[k][j];
	}
	return V;
}

// Multiply two matrices together
// First matrix is (n x k), Second is (k x m)
// Resulting matrix is (n x m)
std::vector<std::vector<double> > CMatrixOperations::MMult(std::vector<std::vector<double> >A, std::vector<std::vector<double> >B, int n, int k, int m) {
	std::vector<std::vector<double> >  C(n, std::vector<double>(m));
	for (int j = 0; j <= m - 1; j++)
		for (int i = 0; i <= n - 1; i++) {
			C[i][j] = 0;
			for (int r = 0; r <= k - 1; r++)
				C[i][j] += A[i][r] * B[r][j];
		}
	return C;
}

// Inverse of a matrix through LU decomposition
std::vector<std::vector<double> >  CMatrixOperations::MInvLU(std::vector<std::vector<double> > A) {
	LUstruct Mats;
	Mats = LU(A);
	std::vector<std::vector<double> >  L = Mats.L;
	std::vector<std::vector<double> >  U = Mats.U;
	std::vector<std::vector<double> >  Linv = MatLowTriangleInv(L);
	std::vector<std::vector<double> >  Uinv = MatUpTriangleInv(U);
	int n = Linv.size();
	int k = Linv[0].size();
	int m = Uinv[0].size();
	std::vector<std::vector<double> >  Ainv = MMult(Uinv, Linv, n, k, m);
	return Ainv;
}
// Transpose of a matrix
std::vector<std::vector<double> > CMatrixOperations::MTrans(std::vector<std::vector<double> >A) {
	int rows = A.size();
	int cols = A[0].size();
	std::vector<std::vector<double> >  C(cols, std::vector<double>(rows));
	for (int i = 0; i<cols; i++)
		for (int j = 0; j<rows; j++)
			C[i][j] = A[j][i];
	return C;
}

// Multiply a matrix by a vector
std::vector<double> CMatrixOperations::MVecMult(std::vector<std::vector<double> > X, std::vector<double> Y) {
	int rows = X.size();
	int cols = X[0].size();
	std::vector<double> XY(rows);
	for (int r = 0; r<rows; r++) {
		XY[r] = 0.0;
		for (int c = 0; c<cols; c++)
			XY[r] += X[r][c] * Y[c];
	}
	return XY;
}

// Regression beta coefficients
std::vector<double> CMatrixOperations::betavec(std::vector<std::vector<double> > X, std::vector<double> Y) {
	int rows = X.size();
	int cols = X[0].size();
	std::vector<std::vector<double> > Xt = MTrans(X);
	std::vector<std::vector<double> > XtX = MMult(Xt, X, cols, rows, cols);
	std::vector<std::vector<double> > XtXinv = MInvLU(XtX);
	std::vector<double> XtY = MVecMult(Xt, Y);
	return MVecMult(XtXinv, XtY);
}



