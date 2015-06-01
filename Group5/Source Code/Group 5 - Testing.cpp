// HKUST MAFS 5220 Spring 2014-15
// Group 5 Project
// Project #8: Compute the unilateral CVA of a FX forward
// Arthors: CHEONG, Siu Sing (20279891); TSUI, Hin Kan (08109319); XU, Zhou Lei (20279889) 
// Description: Test the calculation of the unilateral CVA of a FX forward contract.
//              The computation logic of this program is exactly the same as the program calculating the unlateral CVA of a FX forwward.
//              This program includes 20 model testing cases.

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

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 113858.239031519
	Real diff = abs(FX_Value - 113858.239031519);

	if ((diff / 113858.239031519) <= 0.05)
		std::cout << "The output from the program is within +/- 5% of the result from the excel." << std::endl;
	else
		std::cout << "The output from the program exceeds +/- 5% of the result from the excel." << std::endl;

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
	
	// Check the parameters estimated by this program with those estimated by the excel spreadsheet
	// The Mean reversion parameter (¦Ë) estimated by the excel is 0.065345
	// The Equilibrium rate (¦Ì) estimated by the excel is 0.242915
	// The Volatility estimated by the excel is 0.004748
	if ((abs(lamda_mll / 0.065345 - 1) <= 0.05) && (abs(mu_mll / 0.242915 - 1) <= 0.05) && (abs(sigma_mll / 0.004748 - 1) <= 0.05))
		cout << "The parameters of the Vasicek Interest Model estimated by the program are within +/- 5% of those estimated by the excel." << endl;
	else
		cout << "The parameters of the Vasicek Interest Model estimated by the program exceed +/- 5% of those estimated by the excel." << endl;


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
	int n_sim = 100;

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

		// Simulate 1 day ahead fx rate
		fx_simt0 = 0;

		for (int z = 0; z < n_sim; z++)
		{
			BigInteger seed_fx = SeedGenerator::instance().get();
			MersenneTwisterUniformRng mersenneRng_fx(seed_fx);
			MersenneBoxMuller boxMullerRng_fx(mersenneRng_fx);
			RandomSequenceGenerator<MersenneBoxMuller> gsg_fx(1, boxMullerRng_fx);
			PathGenerator<RandomSequenceGenerator<MersenneBoxMuller>> fxPathGenerator_t0(gbm, dt_1d, 1, gsg_fx, false);
			const Path& fxpath_t0 = fxPathGenerator_t0.next().value;

			fx_simt0 = fx_simt0 + fxpath_t0[1];
		}

		fx_simt0 = fx_simt0 / n_sim;

		// Vasicek interest rate process
		r0 = d_ZeroRate_indexes[0];
		const boost::shared_ptr<StochasticProcess>& ir = boost::shared_ptr<StochasticProcess>(new OrnsteinUhlenbeckProcess(lamda_mll, sigma_mll, r0, mu_mll));

		// Simulate 1 day ahead short rate
		ir_simt0 = 0;

		for (int z = 0; z < n_sim; z++)
		{
			BigInteger seed_ir = SeedGenerator::instance().get();
			MersenneTwisterUniformRng mersenneRng_ir(seed_ir);
			MersenneBoxMuller boxMullerRng_ir(mersenneRng_ir);
			RandomSequenceGenerator<MersenneBoxMuller> gsg_ir(1, boxMullerRng_ir);
			PathGenerator<RandomSequenceGenerator<MersenneBoxMuller>> irPathGenerator_t0(ir, dt_1d, 1, gsg_ir, false);
			const Path& irpath_t0 = irPathGenerator_t0.next().value;

			ir_simt0 = ir_simt0 + irpath_t0[1];
		}

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


	// ************* //
	// Model testing //
	// ************* //

	// Test 1: K = 0

	Real K_test1 = 0;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test1 = FX_MTM(Notional, Pos, S0, df_foreign, K_test1, df);

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
				sim_cf = Notional * Pos * (samplePath_fxsim[timeSteps_sim] - K_test1);

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
	Real CVA_test1 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test1 = CVA_test1 + CVA_sim[i];
	}

	CVA_test1 = CVA_test1 / n_sim;

	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 781667.946445869
	Real fxvalue_test1 = 781667.946445869;
	Real diff_test1 = abs(FX_Value_test1 - fxvalue_test1);

	// Test 1 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test1 / fxvalue_test1) <= 0.05) && (((CVA_test1 >= CVA) && (Pos > 0)) || ((CVA_test1 <= CVA) && (Pos < 0))))
		std::cout << "Test 1: K = 0 passes. FX_Value_test1 = " << FX_Value_test1 << ", CVA_test1 = " << CVA_test1 << std::endl;
	else
		std::cout << "Test 1: K = 0 fails." << std::endl;

	// Test 2: S0 = 0

	Real S0_test2 = 0;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test2 = FX_MTM(Notional, Pos, S0_test2, df_foreign, K, df);

	// Simulate CVA for each path
	for (int j = 0; j < n_sim; j++)
	{
		sumNPV = 0;

		// 1 day
		dt_1d = dayCounter.yearFraction(t0, t0 + 1 * Days);

		typedef BoxMullerGaussianRng<MersenneTwisterUniformRng> MersenneBoxMuller;

		// Geometric Brownian Motion process
		const boost::shared_ptr<StochasticProcess>& gbm = boost::shared_ptr<StochasticProcess>(new GeometricBrownianMotionProcess(S0_test2, mu, sigma));

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
	Real CVA_test2 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test2 = CVA_test2 + CVA_sim[i];
	}

	CVA_test2 = CVA_test2 / n_sim;

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is -667809.70741435
	Real fxvalue_test2 = -667809.70741435;
	Real diff_test2 = abs(FX_Value_test2 - fxvalue_test2);

	// Test 2 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test2 / fxvalue_test2) <= 0.05) && (((CVA_test2 <= CVA) && (Pos > 0)) || ((CVA_test2 >= CVA) && (Pos < 0))))
		std::cout << "Test 2: S0 = 0 passes. FX_Value_test2 = " << FX_Value_test2 << ", CVA_test2 = " << CVA_test2 << std::endl;
	else
		std::cout << "Test 2: S0 = 0 fails." << std::endl;

	// Test 3: Notional Amount = 0

	Real Notional_test3 = 0;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test3 = FX_MTM(Notional_test3, Pos, S0, df_foreign, K, df);

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
				sim_cf = Notional_test3 * Pos * (samplePath_fxsim[timeSteps_sim] - K);

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
	Real CVA_test3 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test3 = CVA_test3 + CVA_sim[i];
	}

	CVA_test3 = CVA_test3 / n_sim;

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 0
	Real fxvalue_test3 = 0;
	Real diff_test3 = abs(FX_Value_test3 - fxvalue_test3);

	// Test 3 Result
	// Check the % difference of fair value calculation and the change of CVA
	if ((diff_test3 == fxvalue_test3) && (CVA_test3 == 0))
		std::cout << "Test 3: Notional Amount = 0 passes. FX_Value_test3 = " << FX_Value_test3 << ", CVA_test3 = " << CVA_test3 << std::endl;
	else
		std::cout << "Test 3: Notional Amount = 0 fails." << std::endl;

	// Test 4: S0 is increased by 30%

	Real S0_test4 = S0 * 1.3;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test4 = FX_MTM(Notional, Pos, S0_test4, df_foreign, K, df);

	// Simulate CVA for each path
	for (int j = 0; j < n_sim; j++)
	{
		sumNPV = 0;

		// 1 day
		dt_1d = dayCounter.yearFraction(t0, t0 + 1 * Days);

		typedef BoxMullerGaussianRng<MersenneTwisterUniformRng> MersenneBoxMuller;

		// Geometric Brownian Motion process
		const boost::shared_ptr<StochasticProcess>& gbm = boost::shared_ptr<StochasticProcess>(new GeometricBrownianMotionProcess(S0_test4, mu, sigma));

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
	Real CVA_test4 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test4 = CVA_test4 + CVA_sim[i];
	}

	CVA_test4 = CVA_test4 / n_sim;

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 348358.62296528
	Real fxvalue_test4 = 348358.62296528;
	Real diff_test4 = abs(FX_Value_test4 - fxvalue_test4);

	// Test 4 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test4 / fxvalue_test4) <= 0.05) && (((CVA_test4 >= CVA) && (Pos > 0)) || ((CVA_test4 <= CVA) && (Pos < 0))))
		std::cout << "Test 4: S0 is increased by 30% passes. FX_Value_test4 = " << FX_Value_test4 << ", CVA_test4 = " << CVA_test4 << std::endl;
	else
		std::cout << "Test 4: S0 is increased by 30% fails." << std::endl;

	// Test 5: S0 is decreased by 30%

	Real S0_test5 = S0 * 0.7;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test5 = FX_MTM(Notional, Pos, S0_test5, df_foreign, K, df);

	// Simulate CVA for each path
	for (int j = 0; j < n_sim; j++)
	{
		sumNPV = 0;

		// 1 day
		dt_1d = dayCounter.yearFraction(t0, t0 + 1 * Days);

		typedef BoxMullerGaussianRng<MersenneTwisterUniformRng> MersenneBoxMuller;

		// Geometric Brownian Motion process
		const boost::shared_ptr<StochasticProcess>& gbm = boost::shared_ptr<StochasticProcess>(new GeometricBrownianMotionProcess(S0_test5, mu, sigma));

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
	Real CVA_test5 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test5 = CVA_test5 + CVA_sim[i];
	}

	CVA_test5 = CVA_test5 / n_sim;

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is -120642.144902242
	Real fxvalue_test5 = -120642.144902242;
	Real diff_test5 = abs(FX_Value_test5 - fxvalue_test5);

	// Test 5 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test5 / fxvalue_test5) <= 0.05) && (((CVA_test5 <= CVA) && (Pos > 0)) || ((CVA_test5 >= CVA) && (Pos < 0))))
		std::cout << "Test 5: S0 is decreased by 30% passes. FX_Value_test5 = " << FX_Value_test5 << ", CVA_test5 = " << CVA_test5 << std::endl;
	else
		std::cout << "Test 5: S0 is decreased by 30% fails." << std::endl;

	// Test 6: Notional amount is increased by 30%

	Real Notional_test6 = Notional * 1.3;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test6 = FX_MTM(Notional_test6, Pos, S0, df_foreign, K, df);

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
				sim_cf = Notional_test6 * Pos * (samplePath_fxsim[timeSteps_sim] - K);

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
	Real CVA_test6 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test6 = CVA_test6 + CVA_sim[i];
	}

	CVA_test6 = CVA_test6 / n_sim;

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 148015.710740975
	Real fxvalue_test6 = 148015.710740975;
	Real diff_test6 = abs(FX_Value_test6 - fxvalue_test6);

	// Test 6 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test6 / fxvalue_test6) <= 0.05) && (((CVA_test6 >= CVA) && (Pos > 0)) || ((CVA_test6 >= CVA) && (Pos < 0))))
		std::cout << "Test 6: Notional amount is increased by 30% passes. FX_Value_test6 = " << FX_Value_test6 << ", CVA_test6 = " << CVA_test6 << std::endl;
	else
		std::cout << "Test 6: Notional amount is increased by 30% fails." << std::endl;



	// Test 7: Notional amount is decreased by 30%

	Real Notional_test7 = Notional * 0.7;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test7 = FX_MTM(Notional_test7, Pos, S0, df_foreign, K, df);

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
				sim_cf = Notional_test7 * Pos * (samplePath_fxsim[timeSteps_sim] - K);

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
	Real CVA_test7 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test7 = CVA_test7 + CVA_sim[i];
	}

	CVA_test7 = CVA_test7 / n_sim;

	// Check the output from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 79700.7673220636
	Real fxvalue_test7 = 79700.7673220636;
	Real diff_test7 = abs(FX_Value_test7 - fxvalue_test7);

	// Test 7 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test7 / fxvalue_test7) <= 0.05) && (((CVA_test7 <= CVA) && (Pos > 0)) || ((CVA_test7 <= CVA) && (Pos < 0))))
		std::cout << "Test 7: Notional amount is decreased by 30% passes. FX_Value_test7 = " << FX_Value_test7 << ", CVA_test7 = " << CVA_test7 << std::endl;
	else
		std::cout << "Test 7: Notional amount is decreased by 30% fails." << std::endl;

	// Test 8: Contract rate is increased by 30%

	Real K_test8 = K * 1.3;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test8 = FX_MTM(Notional, Pos, S0, df_foreign, K_test8, df);

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
				sim_cf = Notional * Pos * (samplePath_fxsim[timeSteps_sim] - K_test8);

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
	Real CVA_test8 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test8 = CVA_test8 + CVA_sim[i];
	}

	CVA_test8 = CVA_test8 / n_sim;

	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is -86484.6731927856
	Real fxvalue_test8 = -86484.6731927856;
	Real diff_test8 = abs(FX_Value_test8 - fxvalue_test8);

	// Test 8 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test8 / fxvalue_test8) <= 0.05) && (((CVA_test8 <= CVA) && (Pos > 0)) || ((CVA_test8 >= CVA) && (Pos < 0))))
		std::cout << "Test 8: Contract rate is increased by 30% passes. FX_Value_test8 = " << FX_Value_test8 << ", CVA_test8 = " << CVA_test8 << std::endl;
	else
		std::cout << "Test 8: Contract rate is increased by 30% fails." << std::endl;

	// Test 9: Contract rate is decreased by 30%

	Real K_test9 = K * 0.7;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test9 = FX_MTM(Notional, Pos, S0, df_foreign, K_test9, df);

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
				sim_cf = Notional * Pos * (samplePath_fxsim[timeSteps_sim] - K_test9);

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
	Real CVA_test9 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test9 = CVA_test9 + CVA_sim[i];
	}

	CVA_test9 = CVA_test9 / n_sim;

	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 314201.151255824
	Real fxvalue_test9 = 314201.151255824;
	Real diff_test9 = abs(FX_Value_test9 - fxvalue_test9);

	// Test 9 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test9 / fxvalue_test9) <= 0.05) && (((CVA_test9 >= CVA) && (Pos > 0)) || ((CVA_test9 <= CVA) && (Pos < 0))))
		std::cout << "Test 9: Contract rate is decreased by 30% passes. FX_Value_test9 = " << FX_Value_test9 << ", CVA_test9 = " << CVA_test9 << std::endl;
	else
		std::cout << "Test 9: Contract rate is decreased by 30% fails." << std::endl;

	// Test 10: Domestic zero rate curve shifts up by 200bps

	std::vector<Real> d_ZeroRate_indexes_test10(32);

	// Shift up the domestic zero rate curve by 200bps
	for (int i = 0; i < 32; i++)
	{
		d_ZeroRate_indexes_test10[i] = d_ZeroRate_indexes[i] + 0.02;
	}

	// Linear Interpolation of the domestic zero rates
	LinearInterpolation IntZero_test10(tenor_indexes.begin(), tenor_indexes.end(), d_ZeroRate_indexes_test10.begin());
	Real r_test10;

	// Assume constant rate when dt is out of range
	if (dt < tenor_indexes[0])
	{
		r_test10 = d_ZeroRate_indexes_test10[0];
	}
	else if (dt > tenor_indexes[31])
	{
		r_test10 = d_ZeroRate_indexes_test10[31];
	}
	else
	{
		r_test10 = IntZero_test10(dt);
	}

	// Discount factor of the domestic currency
	Real df_test10 = std::exp(-1 * r_test10*dt);

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test10 = FX_MTM(Notional, Pos, S0, df_foreign, K, df_test10);
	
	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 133577.242958822
	Real fxvalue_test10 = 133577.242958822;
	Real diff_test10 = abs(FX_Value_test10 - fxvalue_test10);

	// Test 10 Result
	// Check the % difference of fair value calculation
	// Increase in domestic interest rate can impact the interest rate forecast, PD estimation, EAD estimation etc.
	// As a result, it is hard to determine the movement of CVA for upward shift of domestic zero curve.
	if ((diff_test10 / fxvalue_test10) <= 0.05)
		std::cout << "Test 10: Domestic zero rate curve shifts up by 200bps passes. FX_Value_test10 = " << FX_Value_test10 << std::endl;
	else
		std::cout << "Test 10: Domestic zero rate curve shifts up by 200bps fails." << std::endl;

	// Test 11: Domestic zero rate curve shifts down by 200bps

	std::vector<Real> d_ZeroRate_indexes_test11(32);

	// Shift down the domestic zero rate curve by 200bps
	for (int i = 0; i < 32; i++)
	{
		d_ZeroRate_indexes_test11[i] = d_ZeroRate_indexes[i] - 0.02;
	}

	// Linear Interpolation of the domestic zero rates
	LinearInterpolation IntZero_test11(tenor_indexes.begin(), tenor_indexes.end(), d_ZeroRate_indexes_test11.begin());
	Real r_test11;

	// Assume constant rate when dt is out of range
	if (dt < tenor_indexes[0])
	{
		r_test11 = d_ZeroRate_indexes_test11[0];
	}
	else if (dt > tenor_indexes[31])
	{
		r_test11 = d_ZeroRate_indexes_test11[31];
	}
	else
	{
		r_test11 = IntZero_test11(dt);
	}

	// Discount factor of the domestic currency
	Real df_test11 = std::exp(-1 * r_test11*dt);

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test11 = FX_MTM(Notional, Pos, S0, df_foreign, K, df_test11);

	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 93539.2587216009
	Real fxvalue_test11 = 93539.2587216009;
	Real diff_test11 = abs(FX_Value_test11 - fxvalue_test11);

	// Test 11 Result
	// Check the % difference of fair value calculation
	// Decrease in domestic interest rate can impact the interest rate forecast, PD estimation, EAD estimation etc.
	// As a result, it is hard to determine the movement of CVA for downward shift of domestic zero curve.
	if ((diff_test11 / fxvalue_test11) <= 0.05)
		std::cout << "Test 11: Domestic zero rate curve shifts down by 200bps passes. FX_Value_test11 = " << FX_Value_test11 << std::endl;
	else
		std::cout << "Test 11: Domestic zero rate curve shifts down by 200bps fails. " << std::endl;

	// Test 12: Foreign zero rate curve shifts up by 200bps

	std::vector<Real> f_ZeroRate_indexes_test12(32);

	// Shift up the foreign zero rate curve by 200bps
	for (int i = 0; i < 32; i++)
	{
		f_ZeroRate_indexes_test12[i] = f_ZeroRate_indexes[i] + 0.02;
	}

	// Linear interpolation of the foreign zero rates
	LinearInterpolation IntZero_foreign_test12(tenor_indexes.begin(), tenor_indexes.end(), f_ZeroRate_indexes_test12.begin());
	Real r_foreign_test12;

	// Assume constant rate when dt is out of range
	if (dt < tenor_indexes[0])
	{
		r_foreign_test12 = f_ZeroRate_indexes_test12[0];
	}
	else if (dt > tenor_indexes[31])
	{
		r_foreign_test12 = f_ZeroRate_indexes_test12[31];
	}
	else
	{
		r_foreign_test12 = IntZero_foreign_test12(dt);
	}

	// Discount factor of the foreign currency
	Real df_foreign_test12 = std::exp(-1 * r_foreign_test12*dt);

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test12 = FX_MTM(Notional, Pos, S0, df_foreign_test12, K, df);

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
	Real CVA_test12 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test12 = CVA_test12 + CVA_sim[i];
	}

	CVA_test12 = CVA_test12 / n_sim;

	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 90777.2428514964 
	Real fxvalue_test12 = 90777.2428514964;
	Real diff_test12 = abs(FX_Value_test12 - fxvalue_test12);

	// Test 12 Result
	// Check the % difference of fair value calculation and the change of CVA
	// The CVA doesn't reply on the foreign interest rates, but the Change of CVA should be very small but not zero.
	// This is because the seeds are not fixed, it is random. A +/-5% tolerance level is assigned for the CVA change in test 12.
	if (((diff_test12 / fxvalue_test12) <= 0.05) && (abs(CVA_test12/CVA - 1) <= 0.05))
		std::cout << "Test 12: Foreign zero rate curve shifts up by 200bps passes. FX_Value_test12 = " << FX_Value_test12 << ", CVA_test12 = " << CVA_test12 << std::endl;
	else
		std::cout << "Test 12: Foreign zero rate curve shifts up by 200bps fails." << std::endl;

	// Test 13: Foreign zero rate curve shifts down by 200bps

	std::vector<Real> f_ZeroRate_indexes_test13(32);

	// Shift down the foreign zero rate curve by 200bps
	for (int i = 0; i < 32; i++)
	{
		f_ZeroRate_indexes_test13[i] = f_ZeroRate_indexes[i] - 0.02;
	}

	// Linear interpolation of the foreign zero rates
	LinearInterpolation IntZero_foreign_test13(tenor_indexes.begin(), tenor_indexes.end(), f_ZeroRate_indexes_test13.begin());
	Real r_foreign_test13;

	// Assume constant rate when dt is out of range
	if (dt < tenor_indexes[0])
	{
		r_foreign_test13 = f_ZeroRate_indexes_test13[0];
	}
	else if (dt > tenor_indexes[31])
	{
		r_foreign_test13 = f_ZeroRate_indexes_test13[31];
	}
	else
	{
		r_foreign_test13 = IntZero_foreign_test13(dt);
	}

	// Discount factor of the foreign currency
	Real df_foreign_test13 = std::exp(-1 * r_foreign_test13*dt);

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test13 = FX_MTM(Notional, Pos, S0, df_foreign_test13, K, df);

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
	Real CVA_test13 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test13 = CVA_test13 + CVA_sim[i];
	}

	CVA_test13 = CVA_test13 / n_sim;

	// Check the fair value of the FX forward from this program with the result from the excel spreadsheet
	// The fair value calculated by the excel spreadsheet is 137641.504588145
	Real fxvalue_test13 = 137641.504588145;
	Real diff_test13 = abs(FX_Value_test13 - fxvalue_test13);

	// Test 13 Result
	// Check the % difference of fair value calculation and the change of CVA
	// The CVA doesn't reply on the foreign interest rates, but the Change of CVA should be very small but not zero.
	// This is because the seeds are not fixed, it is random. A +/-5% tolerance level is assigned for the CVA change in test 13.
	if (((diff_test13 / fxvalue_test13) <= 0.05) && (abs(CVA_test13 / CVA - 1) <= 0.05))
		std::cout << "Test 13: Foreign zero rate curve shifts down by 200bps passes. FX_Value_test13 = " << FX_Value_test13 << ", CVA_test13 = " << CVA_test13 << std::endl;
	else
		std::cout << "Test 13: Foreign zero rate curve shifts down by 200bps fails." << std::endl;

	// Test 14: Recovery Rate = 0 assuming there is no change to the default rate
	
	Real LGD_test14 = 1;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test14 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

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
				sumNPV = sumNPV + sim_NPV[i] * LGD_test14 * PD[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test14 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test14 = CVA_test14 + CVA_sim[i];
	}

	CVA_test14 = CVA_test14 / n_sim;

	// The recovery rate should not cause any impact to the default-free fair value of the FX forward
	Real fxvalue_test14 = FX_Value;
	Real diff_test14 = abs(FX_Value_test14 - fxvalue_test14);

	// Test 14 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test14 / fxvalue_test14) == 0) && (((CVA_test14 >= CVA) && (Pos > 0)) || ((CVA_test14 >= CVA) && (Pos < 0))))
		std::cout << "Test 14: Recovery rate = 0 assuming there is no change to the default rate passes. FX_Value_test14 = " << FX_Value_test14 << ", CVA_test14 = " << CVA_test14 << std::endl;
	else
		std::cout << "Test 14: Recovery rate = 0 assuming there is no change to the default rate fails." << std::endl;

	// Test 15: Recovery Rate = 1 assuming there is no change to the default rate

	Real LGD_test15 = 0;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test15 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

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
				sumNPV = sumNPV + sim_NPV[i] * LGD_test15 * PD[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test15 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test15 = CVA_test15 + CVA_sim[i];
	}

	CVA_test15 = CVA_test15 / n_sim;

	// The recovery rate should not cause any impact to the default-free fair value of the FX forward
	Real fxvalue_test15 = FX_Value;
	Real diff_test15 = abs(FX_Value_test15 - fxvalue_test15);

	// Test 15 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test15 / fxvalue_test15) == 0) && (((CVA_test15 == 0) && (Pos > 0)) || ((CVA_test15 == 0) && (Pos < 0))))
		std::cout << "Test 15: Recovery Rate = 1 assuming there is no change to the default rate passes. FX_Value_test15 = " << FX_Value_test15 << ", CVA_test15 = " << CVA_test15 << std::endl;
	else
		std::cout << "Test 15: Recovery Rate = 1 assuming there is no change to the default rate fails." << std::endl;

	// Test 16: Recovery Rate is increased by 30% assuming there is no change to the default rate

	Real LGD_test16 = 1 - recovery_rate * 1.3;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test16 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

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
				sumNPV = sumNPV + sim_NPV[i] * LGD_test16 * PD[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test16 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test16 = CVA_test16 + CVA_sim[i];
	}

	CVA_test16 = CVA_test16 / n_sim;

	// The recovery rate should not cause any impact to the default-free fair value of the FX forward
	Real fxvalue_test16 = FX_Value;
	Real diff_test16 = abs(FX_Value_test16 - fxvalue_test16);

	// Test 16 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test16 / fxvalue_test16) == 0) && (((CVA_test16 <= CVA) && (Pos > 0)) || ((CVA_test16 <= CVA) && (Pos < 0))))
		std::cout << "Test 16: Recovery Rate is increased by 30% assuming there is no change to the default rate passes. FX_Value_test16 = " << FX_Value_test16 << ", CVA_test16 = " << CVA_test16 << std::endl;
	else
		std::cout << "Test 16: Recovery Rate is increased by 30% assuming there is no change to the default rate fails." << std::endl;

	// Test 17: Recovery Rate is decreased by 30% assuming there is no change to the default rate

	Real LGD_test17 = 1 - recovery_rate * 0.7;

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test17 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

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
				sumNPV = sumNPV + sim_NPV[i] * LGD_test17 * PD[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test17 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test17 = CVA_test17 + CVA_sim[i];
	}

	CVA_test17 = CVA_test17 / n_sim;

	// The recovery rate should not cause any impact to the default-free fair value of the FX forward
	Real fxvalue_test17 = FX_Value;
	Real diff_test17 = abs(FX_Value_test17 - fxvalue_test17);

	// Test 17 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test17 / fxvalue_test17) == 0) && (((CVA_test17 >= CVA) && (Pos > 0)) || ((CVA_test17 >= CVA) && (Pos < 0))))
		std::cout << "Test 17: Recovery Rate is decreased by 30% assuming there is no change to the default rate passes, FX_Value_test17 = " << FX_Value_test17 << ", CVA_test17 = " << CVA_test17 << std::endl;
	else
		std::cout << "Test 17: Recovery Rate is decreased by 30% assuming there is no change to the default rate fails." << std::endl;

	// Test 18: CDS spreads are increased by 50bps

	// CDS Quotes (1Y,2Y,3Y,4Y,5Y,7Y,10Y)
	std::vector<Real> quoted_spreads_test18(7);
	for (int i = 0; i < 7; i++)
	{
		quoted_spreads_test18[i] = quoted_spreads[i] + 0.005;
	}

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test18 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

	// DefaultProbabilityHelper
	std::vector<boost::shared_ptr<DefaultProbabilityHelper> > instruments_test18;

	for (int k = 0; k < 7; k++) {
		instruments_test18.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
			new SpreadCdsHelper(
			Handle<Quote>(boost::shared_ptr<Quote>(
			new SimpleQuote(quoted_spreads_test18[k]))),
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
	boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >	hazardRateStructure_test18(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(t0, instruments_test18, dayCounter));
	
	std::vector < Real > PD_test18(timeSteps + 1);

	// PD at t (i.e. PD(t-1 < t < t+1))
	for (Size i = 1; i <= timeSteps; i++)
	{
		todaysDate = t0 + i * Days;
		PD_test18[i] = hazardRateStructure_test18->hazardRate(todaysDate) / 365;
	}

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
				sumNPV = sumNPV + sim_NPV[i] * LGD * PD_test18[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test18 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test18 = CVA_test18 + CVA_sim[i];
	}

	CVA_test18 = CVA_test18 / n_sim;

	// CDS spreads should not impact the default-free fair value of the FX forward 
	Real fxvalue_test18 = FX_Value;
	Real diff_test18 = abs(FX_Value_test18 - fxvalue_test18);

	// Test 18 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test18 / fxvalue_test18) == 0) && (((CVA_test18 >= CVA) && (Pos > 0)) || ((CVA_test18 >= CVA) && (Pos < 0))))
		std::cout << "Test 18: CDS spreads are increased by 50bps passes. FX_Value_test18 = " << FX_Value_test18 << ", CVA_test18 = " << CVA_test18 << std::endl;
	else
		std::cout << "Test 18: CDS spreads are increased by 50bps fails." << std::endl;

	// Test 19: CDS spreads are decreased by 50bps

	// CDS Quotes (1Y,2Y,3Y,4Y,5Y,7Y,10Y)
	std::vector<Real> quoted_spreads_test19(7);
	for (int i = 0; i < 7; i++)
	{
		quoted_spreads_test19[i] = quoted_spreads[i] - 0.005;
	}

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test19 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

	// DefaultProbabilityHelper
	std::vector<boost::shared_ptr<DefaultProbabilityHelper> > instruments_test19;

	for (int k = 0; k < 7; k++) {
		instruments_test19.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
			new SpreadCdsHelper(
			Handle<Quote>(boost::shared_ptr<Quote>(
			new SimpleQuote(quoted_spreads_test19[k]))),
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
	boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >	hazardRateStructure_test19(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(t0, instruments_test19, dayCounter));

	std::vector < Real > PD_test19(timeSteps + 1);

	// PD at t (i.e. PD(t-1 < t < t+1))
	for (Size i = 1; i <= timeSteps; i++)
	{
		todaysDate = t0 + i * Days;
		PD_test19[i] = hazardRateStructure_test19->hazardRate(todaysDate) / 365;
	}

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
				sumNPV = sumNPV + sim_NPV[i] * LGD * PD_test19[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test19 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test19 = CVA_test19 + CVA_sim[i];
	}

	CVA_test19 = CVA_test19 / n_sim;

	// CDS spreads should not impact the default-free fair value of the FX forward 
	Real fxvalue_test19 = FX_Value;
	Real diff_test19 = abs(FX_Value_test19 - fxvalue_test19);

	// Test 19 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test19 / fxvalue_test19) == 0) && (((CVA_test19 <= CVA) && (Pos > 0)) || ((CVA_test19 <= CVA) && (Pos < 0))))
		std::cout << "Test 19: CDS spreads are decreased by 50bps passes. FX_Value_test19 = " << FX_Value_test19 << ", CVA_test19 = " << CVA_test19 << std::endl;
	else
		std::cout << "Test 19: CDS spreads are decreased by 50bps fails." << std::endl;

	// Test 20: CDS spreads are increased by 200bps

	// CDS Quotes (1Y,2Y,3Y,4Y,5Y,7Y,10Y)
	std::vector<Real> quoted_spreads_test20(7);

	for (int i = 0; i < 7; i++)
	{
		quoted_spreads_test20[i] = quoted_spreads[i] + 0.02;
	}

	// The fair value of the FX forward at the valuation date
	Real FX_Value_test20 = FX_MTM(Notional, Pos, S0, df_foreign, K, df);

	// DefaultProbabilityHelper
	std::vector<boost::shared_ptr<DefaultProbabilityHelper> > instruments_test20;

	for (int k = 0; k < 7; k++) {
		instruments_test20.push_back(boost::shared_ptr<DefaultProbabilityHelper>(
			new SpreadCdsHelper(
			Handle<Quote>(boost::shared_ptr<Quote>(
			new SimpleQuote(quoted_spreads_test20[k]))),
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
	boost::shared_ptr<PiecewiseDefaultCurve<HazardRate, BackwardFlat> >	hazardRateStructure_test20(new PiecewiseDefaultCurve<HazardRate, BackwardFlat>(t0, instruments_test20, dayCounter));

	std::vector < Real > PD_test20(timeSteps + 1);

	// PD at t (i.e. PD(t-1 < t < t+1))
	for (Size i = 1; i <= timeSteps; i++)
	{
		todaysDate = t0 + i * Days;
		PD_test20[i] = hazardRateStructure_test20->hazardRate(todaysDate) / 365;
	}

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
				sumNPV = sumNPV + sim_NPV[i] * LGD * PD_test20[i] * df_df;
			}

		}

		// CVA for path j
		CVA_sim[j] = sumNPV;

	}


	// Average of the CVA[j]s
	Real CVA_test20 = 0;

	for (int i = 0; i < n_sim; i++)
	{
		CVA_test20 = CVA_test20 + CVA_sim[i];
	}

	CVA_test20 = CVA_test20 / n_sim;

	// CDS spreads should not impact the default-free fair value of the FX forward 
	Real fxvalue_test20 = FX_Value;
	Real diff_test20 = abs(FX_Value_test20 - fxvalue_test20);

	// Test 20 Result
	// Check the % difference of fair value calculation and the change of CVA
	if (((diff_test20 / fxvalue_test20) == 0) && (((CVA_test20 >= CVA) && (Pos > 0)) || ((CVA_test20 >= CVA) && (Pos < 0))))
		std::cout << "Test 20: CDS spreads are increased by 200bps passes. FX_Value_test20 = " << FX_Value_test20 << ", CVA_test20 = " << CVA_test20 << std::endl;
	else
		std::cout << "Test 20: CDS spreads are increased by 200bps fails." << std::endl;


	// ******************** //
	// End of model testing //
	// ******************** //

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	system("pause");
	return 0;
}