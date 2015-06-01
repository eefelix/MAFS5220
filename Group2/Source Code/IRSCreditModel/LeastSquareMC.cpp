#include "pch.h"
#include"LeastSquareMC.h"

using namespace std;

LeastSquareMC::LeastSquareMC(int M, int N, double initPrice, double T, double strike,
	double rate, double vol, Type type)
	: M_(M), N_(N), initPrice_(initPrice), T_(T), strike_(strike), rate_(rate), vol_(vol), type_(type)
{
	dt_ = T / N;
}

double LeastSquareMC::max(double x, double y){
	return x > y ? x : y;
}

//generate prices at every time step along each path
double LeastSquareMC::calcLeastSquareMC(){
	double lnS, S;
	int i, j;
	BigInteger seed = SeedGenerator::instance().get();
	MersenneTwisterUniformRng unifMt(seed);
	BoxMullerGaussianRng<MersenneTwisterUniformRng> bmGauss(unifMt);
	matrix<double> Paths(M_, N_);     //record stock price
	matrix<double> Exercise(M_, N_);  //record exercise payoff
	matrix<double> CashFlow(M_, N_);
	double mu_ = (rate_ - vol_*vol_ / 2)*dt_;
	double voldt = vol_*sqrt(dt_);

	for (i = 0; i < M_; i++)
	{
		lnS = log(initPrice_);
		for (j = 0; j < N_; j++)
		{
			Real normNumber = bmGauss.next().value;
			lnS = lnS + mu_ + voldt*normNumber;
			S = exp(lnS);
			prices[i].push_back(S);
		}
	}

	for (i = 0; i < M_; i++)
		for (j = 0; j < N_; j++)
		{
			Paths(i, j) = prices[i][j];
			Exercise(i, j) = max((prices[i][j] - strike_)*type_, 0);
		}

	//compute cash flows at maturity
	for (i = 0; i < M_; i++)
		CashFlow(i, N_ - 1) = Exercise(i, N_ - 1);

	return computeLSM(N_ - 2, Paths, CashFlow, Exercise);
}


double LeastSquareMC::computeLSM(int time, matrix<double>& Paths, matrix<double>& CashFlow, matrix<double>& Exercise)
{
	double val = 0.0;
	std::vector<int> num(M_); //the initial value of every element in num is 0
	int k = 0; // k records the number of path that is in the money
	for (int j = 0; j < M_; j++){
		prices[j].pop_back();
		val = prices[j].back();
		if ((val - strike_)*type_> 0)
		{
			num[j] = j + 1;
			k++;
		}
	}

	std::vector<double> Y(k); //dependent variable
	std::vector<double> B(k); //coefficient
	std::vector<double> C(k); //continuation
	std::vector<std::vector<double> >  X(k, std::vector<double>(3));

	int i = 0;
	for (int j = 0; j < M_; j++){
		if (i < k && Exercise(j, time) > 0){
			double cashAmount = 0.0;
			int temp = 0;
			for (int t = time + 1; t < N_; t++){
				if (CashFlow(j, t) != 0)
					cashAmount = CashFlow(j, t);
				temp = t;
			}

			if (temp == 0)
				Y[i] = 0.0;
			else
				Y[i] = cashAmount*exp(-rate_*dt_*(temp - time)); //problem

			X[i][0] = 1;
			X[i][1] = Paths(j, time);
			X[i][2] = Paths(j, time)*Paths(j, time);

			i++;
		}
	}

	B = betavec(X, Y);
	C = MVecMult(X, B);

	i = 1;
	for (int j = 0; j < M_; j++){
		if (num[j] != 0){
			if (Exercise(j, time) > C[i]){
				CashFlow(j, time) = Exercise(j, time);
				for (int t = time + 1; t < N_; t++){
					CashFlow(j, t) = 0.0; //bug
				}
			}
			else CashFlow(j, time) = 0.0;
		}
		else CashFlow(j, time) = 0.0;
	}
	if (time > 0)
		computeLSM(time - 1, Paths, CashFlow, Exercise);

	return computeValue(CashFlow);
}

double LeastSquareMC::computeValue(matrix<double>& CashFlow){
	int i, j;
	double discValue = 0.0;
	for (i = 0; i < M_; i++){
		for (j = 0; j < N_; j++){
			if (CashFlow(i, j) > 0)
				discValue = discValue + CashFlow(i, j)*exp(-rate_*dt_*(j + 1)); // bug
		}
	}
	return discValue / M_;
}