// HKUST MAFS 5220 Spring 2014-15
// Group 5 Project
// Project #8: Compute the unilateral CVA of a FX forward
// Arthor: CHEONG, Siu Sing (20279891); TSUI, Hin Kan (08109319); XU, Zhou Lei (20279889) 
// Description: Calculate the unilateral CVA of a FX forward contract

#pragma warning(disable: 4819)
#pragma warning (disable:4244)
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <ql/methods/montecarlo/path.hpp>
#include <ql/methods/montecarlo/pathgenerator.hpp>
#include <ql/stochasticprocess.hpp>
#include <ql/processes/geometricbrownianprocess.hpp>
#include <ql/math/randomnumbers/seedgenerator.hpp>
#include <ql/math/randomnumbers/randomsequencegenerator.hpp>
#include <ql/math/randomnumbers/boxmullergaussianrng.hpp>
#include <ql/math/randomnumbers/mt19937uniformrng.hpp>
#include <vector>
#include <ql/processes/ornsteinuhlenbeckprocess.hpp>
#include <fstream>
#include <numeric>
#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>
#include <time.h>

using namespace QuantLib;
using namespace std;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

	Integer sessionId() { return 0; }

}
#endif

//  Function to calculate the fair value of the FX forward
Real FX_MTM(Real not, Real p, Real s_0, Real df_for, Real strike, Real df_do)
{
	Real mtm;
	mtm = not * p * (s_0 * df_for - strike * df_do);
	return mtm;
}

int main(int, char*[]){

	boost::timer timer;

	clock_t tStart = clock();

	// Calendar
	Calendar calendar = TARGET();

	// Day count convention
	DayCounter dayCounter = Actual365Fixed();

	// ********** //
	// User input //
	// ********** //

	// Valuation date
	Date t0(31, Dec, 2014);
	t0 = calendar.adjust(t0);

	// Maturity date
	Date T(30, Jun, 2016);
	Real dt = dayCounter.yearFraction(t0, T);

	// Position of the FX forward - long (1) or short (-1)
	Real Pos = 1;

	// FX spot at valuation date
	Real S0 = 0.8198;

	// Contract rate of the FX forward
	Real K = 0.7;

	// Notional amount of the FX forward
	Real Notional = 1000000;

	// Drift and the sigma in the Geometric Brownian Motion of FX rates
	// They are derived by using 1 year historical log return data of the FX rate
	Real mu = 0.000324;
	Volatility sigma = 0.0055;

	Settings::instance().evaluationDate() = t0;

	// Recovery Assumption
	Real recovery_rate = 0.4;

	// Loss given default
	Real LGD = 1 - recovery_rate;

	// CDS Quotes (1Y,2Y,3Y,4Y,5Y,7Y,10Y)
	Real quoted_spreads[] = { 0.00675, 0.00845, 0.01005, 0.0120, 0.013757, 0.016152, 0.017866 };

	// ***************** //
	// End of user input //
	// ***************** //


	// ******************************* //
	// Bootstrapping of the zero rates //
	// ******************************* //

	//collect the BondPrice, Coupon and tenor data
	ifstream filein("BootStrapping.csv");
	int t;
	double d_BondPrice;
	double d_Coupon;
	double f_BondPrice;
	double f_Coupon;
	double tenor;

	//construct the vectors
	vector<Real> d_BondPrice_indexes;
	vector<Real> d_Coupon_indexes;
	vector<Real> f_BondPrice_indexes;
	vector<Real> f_Coupon_indexes;
	vector<Real> tenor_indexes;

	for (string line; getline(filein, line);)
	{
		stringstream lineStream(line);
		string cell;
		int cell_idx = 0;

		while (getline(lineStream, cell, ',')){

			//classify the cell into different fields
			switch (cell_idx){

			case 0:
				t = atoi(cell.c_str());
				if (t < 1 || t>32) t = -1; //Mark the invalid entry as -1;
				break;
			case 1:
				d_BondPrice = atof(cell.c_str());
				break;
			case 2:
				d_Coupon = atof(cell.c_str());
				break;
			case 3:
				f_BondPrice = atof(cell.c_str());
				break;
			case 4:
				f_Coupon = atof(cell.c_str());
				break;
			case 5:
				tenor = atof(cell.c_str());
				break;
			}
			cell_idx++;
		}

		//skip the invalid data
		if (t == -1) continue;
		d_BondPrice_indexes.push_back(d_BondPrice);
		d_Coupon_indexes.push_back(d_Coupon);
		f_BondPrice_indexes.push_back(f_BondPrice);
		f_Coupon_indexes.push_back(f_Coupon);
		tenor_indexes.push_back(tenor);
	}

	vector<Real> d_ZeroRate_indexes(32);
	vector<Real> f_ZeroRate_indexes(32);

	// bootstrapping the domestic zero rates
	d_ZeroRate_indexes[0] = log(d_BondPrice_indexes[0] / 100) / -tenor_indexes[0];
	d_ZeroRate_indexes[1] = log(d_BondPrice_indexes[1] / 100) / -tenor_indexes[1];
	d_ZeroRate_indexes[2] = log(d_BondPrice_indexes[2] / 100) / -tenor_indexes[2];

	double d_acc = exp(-d_ZeroRate_indexes[2] * tenor_indexes[2]);

	for (Size i = 3; i < d_BondPrice_indexes.size(); i++) {
		d_ZeroRate_indexes[i] = log((d_BondPrice_indexes[i] - 0.5*d_Coupon_indexes[i] * d_acc) / (d_Coupon_indexes[i] * 0.5 + 100)) / -tenor_indexes[i];
		if (i == d_BondPrice_indexes.size() - 1) break;
		d_acc = d_acc + exp(-d_ZeroRate_indexes[i] * tenor_indexes[i]);
	}

	// bootstrapping the foreign zero rates
	f_ZeroRate_indexes[0] = log(f_BondPrice_indexes[0] / 100) / -tenor_indexes[0];
	f_ZeroRate_indexes[1] = log(f_BondPrice_indexes[1] / 100) / -tenor_indexes[1];
	f_ZeroRate_indexes[2] = log(f_BondPrice_indexes[2] / 100) / -tenor_indexes[2];

	double f_acc = exp(-f_ZeroRate_indexes[2] * tenor_indexes[2]);

	for (Size i = 3; i < f_BondPrice_indexes.size(); i++) {
		f_ZeroRate_indexes[i] = log((f_BondPrice_indexes[i] - 0.5*f_Coupon_indexes[i] * f_acc) / (f_Coupon_indexes[i] * 0.5 + 100)) / -tenor_indexes[i];
		if (i == d_BondPrice_indexes.size() - 1) break;
		f_acc = f_acc + exp(-f_ZeroRate_indexes[i] * tenor_indexes[i]);
	}


	// ******************************* //
	// End of bootstrapping zero rates //
	// ******************************* //
	

	// **************************************** //
	// Fair value calculation of the FX forward //
	// **************************************** //

	// Linear Interpolation of the domestic zero rates
	LinearInterpolation IntZero(tenor_indexes.begin(), tenor_indexes.end(), d_ZeroRate_indexes.begin());
	Real r;

	// Assume constant rate when dt is out of range
	if (dt < tenor_indexes[0])
	{
		r = d_ZeroRate_indexes[0];
	}
	else if (dt > tenor_indexes[31])
	{
		r = d_ZeroRate_indexes[31];
	}
	else
	{
		r = IntZero(dt);
	}


	// Linear interpolation of the foreign zero rates
	LinearInterpolation IntZero_foreign(tenor_indexes.begin(), tenor_indexes.end(), f_ZeroRate_indexes.begin());
	Real r_foreign;

	// Assume constant rate when dt is out of range
	if (dt < tenor_indexes[0])
	{
		r_foreign = f_ZeroRate_indexes[0];
	}
	else if (dt > tenor_indexes[31])
	{
		r_foreign = f_ZeroRate_indexes[31];
	}
	else
	{
		r_foreign = IntZero_foreign(dt);
	}

	// Discount factor of the domestic currency
	Real df = std::exp(-1 * r*dt);

	// Discount factor of the foreign currency
	Real df_foreign = std::exp(-1 * r_foreign*dt);

	// The fair value of the FX forward at the valuation date
	Real FX_Value = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

	// Output the fair value of the FX forward at the valuation date
	std::cout << "Fair Value of the FX forward = " << FX_Value << std::endl;


	// *********************************************** //
	// End of fair value calculation of the FX forward //
	// *********************************************** //


	// **************************** //
	// Calibration of Vasicek model //
	// **************************** //

	double delta_t = 0.0833333333333333;  //1 month time increment between observations

	// Common statistics for maximum likelihood method
	double Sx = 0;
	double Sy = 0;
	double Sxx = 0;
	double Sxy = 0;
	double Syy = 0;

	// Load the interest rate data from a csv file
	ifstream filein2("data2.csv");
	int z;
	double Int_rate;
	double Int_rate_lag;

	// Construct the interest rate indexes and interest rate lag indexes
	vector<Real> int_rate_indexes;
	vector<Real> int_rate_lag_indexes;

	for (string line; getline(filein2, line);)
	{
		stringstream lineStream(line);
		string cell;
		int cell_idx = 0;

		while (getline(lineStream, cell, ',')){

			// Classify the cell into different fields
			switch (cell_idx){

			case 0:
				z = atoi(cell.c_str());
				if (z < 1 || z>253) z = -1; //Mark the invalid entry as -1;
				break;
			case 1:
				Int_rate = atof(cell.c_str());
				break;
			case 2:
				Int_rate_lag = atof(cell.c_str());
				break;
			}
			cell_idx++;
		}

		//skip the invalid data
		if (z == -1) continue;
		int_rate_indexes.push_back(Int_rate);
		int_rate_lag_indexes.push_back(Int_rate_lag);
	}

	// Calculate the common statistcs for the MLL method
	for (Size i = 0; i < int_rate_indexes.size() - 1; i++) {
		Sx = Sx + int_rate_indexes[i];
		Sxx = Sxx + (int_rate_indexes[i] * int_rate_indexes[i]);

	}

	for (Size i = 0; i < int_rate_lag_indexes.size() - 1; i++) {
		Sy = Sy + int_rate_lag_indexes[i];
		Syy = Syy + (int_rate_lag_indexes[i] * int_rate_lag_indexes[i]);

	}

	for (Size i = 0; i < int_rate_indexes.size() - 1; i++) {
		Sxy = Sxy + (int_rate_indexes[i] * int_rate_lag_indexes[i]);
	}

	//maximum log-likelihood function method
	double alpha;
	double lamda_mll;
	double mu_mll;
	double sigmah2_mll;
	double sigma_mll;
	int n2;
	n2 = int_rate_indexes.size() - 1;
	mu_mll = (Sy*Sxx - Sx*Sxy) / (n2*(Sxx - Sxy) - (Sx*Sx - Sx*Sy));
	lamda_mll = -log((Sxy - mu_mll*Sx - mu_mll*Sy + n2*mu_mll*mu_mll) / (Sxx - 2 * mu_mll*Sx + n2*mu_mll*mu_mll)) / delta_t;
	alpha = exp(-lamda_mll*delta_t);
	sigmah2_mll = (Syy - 2 * alpha*Sxy + alpha*alpha * Sxx - 2 * mu_mll*(1 - alpha)*(Sy - alpha*Sx) + n2*mu_mll*mu_mll * (1 - alpha)*(1 - alpha)) / n2;
	sigma_mll = sqrt(sigmah2_mll * 2 * lamda_mll / (1 - alpha*alpha));
	

	// *********************************** //
	// End of calibration of Vasicek model //
	// *********************************** //


	// ******************//
	// Hazard Rate Model //
	// ******************//
	// Vector for PD estimation

	vector<Period> tenors;
	tenors.push_back(1 * Years);
	tenors.push_back(2 * Years);
	tenors.push_back(3 * Years);
	tenors.push_back(4 * Years);
	tenors.push_back(5 * Years);
	tenors.push_back(7 * Years);
	tenors.push_back(10 * Years);

	vector<Date> dates;
	vector<Rate> rates;
	dates.push_back(t0); rates.push_back(d_ZeroRate_indexes[0]); // O/N rate
	dates.push_back(t0 + 6 * Months); rates.push_back(d_ZeroRate_indexes[2]); // 6M Zero Rate
	dates.push_back(t0 + 1 * Years); rates.push_back(d_ZeroRate_indexes[3]); // 1Y Zero Rate
	dates.push_back(t0 + 2 * Years); rates.push_back(d_ZeroRate_indexes[5]); // 2Y Zero Rate
	dates.push_back(t0 + 3 * Years); rates.push_back(d_ZeroRate_indexes[7]); // 3Y Zero Rate
	dates.push_back(t0 + 4 * Years); rates.push_back(d_ZeroRate_indexes[9]); // 4Y Zero Rate
	dates.push_back(t0 + 5 * Years); rates.push_back(d_ZeroRate_indexes[11]); // 5Y Zero Rate
	dates.push_back(t0 + 7 * Years); rates.push_back(d_ZeroRate_indexes[15]); // 7Y Zero Rate
	dates.push_back(t0 + 10 * Years); rates.push_back(d_ZeroRate_indexes[21]); // 10Y Zero Rate
	dates.push_back(t0 + 15 * Years); rates.push_back(d_ZeroRate_indexes[31]); // 15Y Zero Rate 

	Handle<YieldTermStructure> irTermStructure(boost::shared_ptr<InterpolatedZeroCurve<Linear> >(new InterpolatedZeroCurve<Linear>(dates, rates, dayCounter)));

	// DefaultProbabilityHelper
	std::vector<boost::shared_ptr<DefaultProbabilityHelper> > instruments;
	
	for (int k = 0; k < 7; k++) {
		instruments.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
			new SpreadCdsHelper(
			Handle<Quote>(boost::shared_ptr<Quote>(
			new SimpleQuote(quoted_spreads[k]))),
			tenors[k],
			0,
			calendar,
			Quarterly,
			Following,
			DateGeneration::TwentiethIMM,
			dayCounter,
			recovery_rate,
			irTermStructure, 1, 1)));
	}

	// Bootstrap hazard rates
	boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >	hazardRateStructure(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(t0, instruments, dayCounter));


	// *********************** //
	// End of Hazard Rate Model//
	// *********************** //

	
	// ************* //
	// PD Estimation //
	// ************* //

	Size timeSteps = dayCounter.dayCount(t0, T);

	std::vector < Real > PD(timeSteps + 1);

	Date todaysDate = t0;

	// PD at t (i.e. PD(t-1 < t < t+1))
	for (Size i = 1; i <= timeSteps; i++)
	{
		todaysDate = t0 + i * Days;
		PD[i] = hazardRateStructure->hazardRate(todaysDate) / 365;
	}

	// ******************** //
	// End of PD Estimation //
	// ******************** //


	// ***************** //
	// Simulation of CVA //
	// ***************** //

	// No of simulation
	int n_sim = 1000;

	// Construct the vectors and variable used for EAD simulation
	Date defaultDate;
	Real r0;
	Real dt_1d;
	Real dt_sim;
	Real dt_df;
	Real r_df;
	Real df_df;
	Real sumNPV;
	Real ir_simt0;
	Real ir_simt1;
	Real fx_simt0;
	Real fx_simt1;
	Real sim_cf;
	Real sim_ir;
	Real sim_df;
	Real sim_EAD;
	std::vector < Real > sim_NPV(timeSteps);
	std::vector < Real > CVA_sim(n_sim);	

	// Simulate CVA for each path
	for (int j = 0; j < n_sim; j++)
	{
		sumNPV = 0;

		// 1 day
		dt_1d = dayCounter.yearFraction(t0, t0 + 1 * Days);

		typedef BoxMullerGaussianRng<MersenneTwisterUniformRng> MersenneBoxMuller;

		// Geometric Brownian Motion process
		const boost::shared_ptr<StochasticProcess>& gbm = boost::shared_ptr<StochasticProcess>(new GeometricBrownianMotionProcess(S0, mu, sigma));

		// Vasicek interest rate process
		r0 = d_ZeroRate_indexes[0];
		const boost::shared_ptr<StochasticProcess>& ir = boost::shared_ptr<StochasticProcess>(new OrnsteinUhlenbeckProcess(lamda_mll, sigma_mll, r0, mu_mll));

		// Simulate 1 day ahead fx rate
		fx_simt0 = 0;

		// Simulate 1 day ahead short rate
		ir_simt0 = 0;

		for (int z = 0; z < n_sim; z++)
		{
			BigInteger seed_fx = SeedGenerator::instance().get();
			MersenneTwisterUniformRng mersenneRng_fx(seed_fx);
			MersenneBoxMuller boxMullerRng_fx(mersenneRng_fx);
			RandomSequenceGenerator<MersenneBoxMuller> gsg_fx(1, boxMullerRng_fx);
			PathGenerator<RandomSequenceGenerator<MersenneBoxMuller>> fxPathGenerator_t0(gbm, dt_1d, 1, gsg_fx, false);
			const Path& fxpath_t0 = fxPathGenerator_t0.next().value;

			fx_simt0 = fx_simt0 + fxpath_t0[1];

			BigInteger seed_ir = SeedGenerator::instance().get();
			MersenneTwisterUniformRng mersenneRng_ir(seed_ir);
			MersenneBoxMuller boxMullerRng_ir(mersenneRng_ir);
			RandomSequenceGenerator<MersenneBoxMuller> gsg_ir(1, boxMullerRng_ir);
			PathGenerator<RandomSequenceGenerator<MersenneBoxMuller>> irPathGenerator_t0(ir, dt_1d, 1, gsg_ir, false);
			const Path& irpath_t0 = irPathGenerator_t0.next().value;

			ir_simt0 = ir_simt0 + irpath_t0[1];
			
		}

		fx_simt0 = fx_simt0 / n_sim;
		ir_simt0 = ir_simt0 / n_sim;

		// Simulate the EAD on a daily basis
		for (Size i = 1; i < timeSteps; i++)
		{
			
			defaultDate = t0 + i * Days;
			defaultDate = calendar.adjust(defaultDate);
			Settings::instance().evaluationDate() = defaultDate;

			// Default date to maturity
			dt_sim = dayCounter.yearFraction(defaultDate, T);

			Size timeSteps_sim = dayCounter.dayCount(defaultDate, T);

			// Update the Geometric Brownian Motion process
			const boost::shared_ptr<StochasticProcess>& gbm_sim = boost::shared_ptr<StochasticProcess>(new GeometricBrownianMotionProcess(fx_simt0, mu, sigma));

			// Update the Vasicek Interest Rate Model
			const boost::shared_ptr<StochasticProcess>& ir_sim = boost::shared_ptr<StochasticProcess>(new OrnsteinUhlenbeckProcess(lamda_mll, sigma_mll, ir_simt0, mu_mll));

			// Simulate fx rates from the default date
			fx_simt1 = 0;

			// Simulate short rates from the default date
			ir_simt1 = 0;

			sim_EAD = 0;

			for (int z = 0; z < n_sim; z++)
			{
				BigInteger seed_fxsim = SeedGenerator::instance().get();
				MersenneTwisterUniformRng mersenneRng_fxsim(seed_fxsim);
				MersenneBoxMuller boxMullerRng_fxsim(mersenneRng_fxsim);
				RandomSequenceGenerator<MersenneBoxMuller> gsg_fxsim(timeSteps_sim, boxMullerRng_fxsim);
				PathGenerator<RandomSequenceGenerator<MersenneBoxMuller>> fxPathGenerator_sim(gbm_sim, dt_sim, timeSteps_sim, gsg_fxsim, false);
				const Path& samplePath_fxsim = fxPathGenerator_sim.next().value;

				// Calculate the simulated cashflow at the expiration
				sim_cf = Notional * Pos * (samplePath_fxsim[timeSteps_sim] - K);

				BigInteger seed_irsim = SeedGenerator::instance().get();
				MersenneTwisterUniformRng mersenneRng_irsim(seed_irsim);
				MersenneBoxMuller boxMullerRng_irsim(mersenneRng_irsim);
				RandomSequenceGenerator<MersenneBoxMuller> gsg_irsim(1, boxMullerRng_irsim);
				PathGenerator<RandomSequenceGenerator<MersenneBoxMuller>> irPathGenerator_sim(ir_sim, dt_sim, 1, gsg_irsim, false);
				const Path& samplePath_irsim = irPathGenerator_sim.next().value;

				// The unit of the interest rate in the Vasicek model is %, not decimal place
				sim_ir = samplePath_irsim[1] / 100;

				// Calculate stochastic discount factor from default date to maturity
				sim_df = exp(-1 * ((mu_mll - 0.5 * (sigma_mll / lamda_mll) * (sigma_mll / lamda_mll)) * dt_sim + (sim_ir - mu_mll + (sigma_mll / lamda_mll) * (sigma_mll / lamda_mll)) * (1 - exp(-lamda_mll * dt_sim)) / lamda_mll - 0.5 * (sigma_mll / lamda_mll) * (sigma_mll / lamda_mll) * 0.5 * (1 - exp(-2 * lamda_mll * dt_sim)) / lamda_mll));

				// Exposure at default
				sim_EAD = sim_EAD + sim_cf * sim_df;

				// Update the initial fx rate of the GBM process for the next simulation
				fx_simt1 = fx_simt1 + samplePath_fxsim[1];

				// Update the initial short rate of the interest rate process for the next simulation
				ir_simt1 = ir_simt1 + samplePath_irsim[1];
			}

			fx_simt0 = fx_simt1 / n_sim;
			ir_simt0 = ir_simt1 / n_sim;
			sim_NPV[i] = sim_EAD / n_sim;
			
			// Valuation date t0 to default date
			dt_df = dayCounter.yearFraction(t0, defaultDate);

			// Assume constant rate when dt_df is out of range
			if (dt_df < tenor_indexes[0])
			{
				r_df = d_ZeroRate_indexes[0];
			}
			else if (dt_df > tenor_indexes[31])
			{
				r_df = d_ZeroRate_indexes[31];
			}
			else
			{
				r_df = IntZero(dt_df);
			}

			// Discount factor from valuation date t0 to default date
			df_df = std::exp(-1 * r_df * dt_df);

			// Only positive EAD is taken into account of the CVA calculation
			if (sim_NPV[i] > 0)
			{
				sumNPV = sumNPV + sim_NPV[i] * LGD * PD[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}

	// Average of the CVA[j]s
	Real CVA = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA = CVA + CVA_sim[i];
	}

	CVA = CVA / n_sim;

	// Output the CVA
	std::cout << "CVA of the FX forward = " << CVA << std::endl;

	// ********************** //
	// End of CVA calculation //
	// ********************** //

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	system("pause");
	return 0;
}