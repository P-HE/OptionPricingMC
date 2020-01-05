//Part1 of the end Term Assignment. Some of the code is from the code sample written by Prof. RS

#include <iostream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
using namespace std;
double up_factor, uptick_prob, risk_free_rate, strike_price, barrier_price;
double initial_stock_price, expiration_time, volatility, R;
double dt;
int no_of_divisions, no_of_trials;
int no_of_time_steps;

unordered_map<int, double> price_table_call;
unordered_map<int, double> price_table_put;

int hash_pairing(int a, int b) {
	return a >= b ? a * a + a + b : a + b * b;
}

double get_uniform() {
	return (double)rand() / (double)RAND_MAX;
}

//It
void monte_carlo_simulation(int no_of_divisions, int no_of_trials) {
	double R = (risk_free_rate - 0.5 * pow(volatility, 2)) * dt;
	double SD = volatility * sqrt(dt);

	double put_option_price = 0.0;
	double call_option_price = 0.0;

	double call_option_price_2 = 0.0;
	double put_option_price_2 = 0.0;

	for (int i = 0; i < no_of_trials; i++) {
		//initialize to the starting point of the price path
		double S1 = initial_stock_price;
		double S2 = initial_stock_price;
		double S3 = initial_stock_price;
		double S4 = initial_stock_price;
		
		//Boolean to indicate if a curtain price path is killed
		bool o1 = 1;
		bool o2 = 1;
		bool o3 = 1;
		bool o4 = 1;

		for (int j = 0; j < no_of_divisions; j++) {
			double x = get_uniform();
			double y = get_uniform();
			
			//In case it generates a 0 and cause box-muller to geneate an INF
			while (x * y == 0) {
				x = get_uniform();
				y = get_uniform();
			}

			//box-muller
			double a = sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
			double b = sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);

			//S1 = S1 * exp(R + SD * a);
			//S2 = S2 * exp(R - SD * a);
			//S1 = S1 * exp(R + SD * b);
			//S4 = S4 * exp(R - SD * b);
			

			S1 = S1 * exp(R + SD * a);
			if (S1 <= barrier_price) {
				o1 = 0;
			}



			S2 = S2 * exp(R + SD * a);
			if (S2 <= barrier_price) {
				o2 = 0;
			}



			S3 = S3 * exp(R + SD * b);
			if (S3 <= barrier_price) {
				o3 = 0;
			}


			S4 = S4 * exp(R - SD * b);
			if (S4 <= barrier_price) {
				o4 = 0;
			}

			//cout << x << y << endl;
			//cout << a << b << endl;
		}

		//cout << S1 << S2 << S3 << S4 << endl;

		double p1, p2, p3, p4;
		double c1, c2, c3, c4;
		if (o1 == 1) {
			p1 = max(0.0, S1 - strike_price);
		}
		else {
			p1 = 0;
		}

		if (o2 == 1) {
			p2 = max(0.0, S2 - strike_price);
		}
		else {
			p2 = 0;
		}

		if (o3 == 1) {
			p3 = max(0.0, S3 - strike_price);
		}
		else {
			p3 = 0;
		}

		if (o4 == 1) {
			p4 = max(0.0, S4 - strike_price);
		}
		else {
			p4 = 0;
		}

		call_option_price += (p1 + p2 + p3 + p4) / 4.0 / no_of_trials;
		if (call_option_price == INFINITY) {
			cout << p1 << p2 << p3 << p4 << endl;
			cout << S1 << S2 << S3 << S4 << endl;
			break;
		}

		if (o1 == 1) {
			c1 = max(0.0, strike_price - S1);
		}
		else {
			c1 = 0;
		}

		if (o2 == 1) {
			c2 = max(0.0, strike_price - S2);
		}
		else {
			c2 = 0;
		}

		if (o3 == 1) {
			c3 = max(0.0, strike_price - S3);
		}
		else {
			c3 = 0;
		}

		if (o4 == 1) {
			c4 = max(0.0, strike_price - S4);
		}
		else {
			c4 = 0;
		}

		put_option_price += (c1 + c2 + c3 + c4) / 4.0 / no_of_trials;
		//cout << call_option_price << endl;
		//cout << put_option_price << endl;

		//Part 2
		double probability_stock_path_has_hit_barrier_by_the_time_it_got_here;
		double c1_2, c2_2, c3_2, c4_2;
		double p1_2, p2_2, p3_2, p4_2;
		
		probability_stock_path_has_hit_barrier_by_the_time_it_got_here =
			exp((-2.0 / (volatility * volatility * expiration_time)) *
				log(initial_stock_price / barrier_price) *
				log(S1 / barrier_price));
		c1_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
			max(0.0, (S1 - strike_price)));
		if (S1 > barrier_price) {
			p1_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
				max(0.0, (strike_price - S1)));
		}
		else {
			p1_2 = 0.0;
		}

		probability_stock_path_has_hit_barrier_by_the_time_it_got_here =
			exp((-2.0 / (volatility * volatility * expiration_time)) *
				log(initial_stock_price / barrier_price) *
				log(S2 / barrier_price));
		c2_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
			max(0.0, (S2 - strike_price)));
		if (S2 > barrier_price) {
			p2_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
				max(0.0, (strike_price - S2)));
		}
		else {
			p2_2 = 0.0;
		}

		probability_stock_path_has_hit_barrier_by_the_time_it_got_here =
			exp((-2.0 / (volatility * volatility * expiration_time)) *
				log(initial_stock_price / barrier_price) *
				log(S3 / barrier_price));
		c3_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
			max(0.0, (S3 - strike_price)));
		if (S3 > barrier_price) {
			p3_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
				max(0.0, (strike_price - S3)));
		}
		else {
			p3_2 = 0.0;
		}

		probability_stock_path_has_hit_barrier_by_the_time_it_got_here =
			exp((-2.0 / (volatility * volatility * expiration_time)) *
				log(initial_stock_price / barrier_price) *
				log(S4 / barrier_price));
		c4_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
			max(0.0, (S4 - strike_price)));
		if (S4 > barrier_price) {
			p4_2 = ((1 - probability_stock_path_has_hit_barrier_by_the_time_it_got_here) *
				max(0.0, (strike_price - S4)));
		}
		else {
			p4_2 = 0.0;
		}

		call_option_price_2 += (c1_2 + c2_2 + c3_2 + c4_2) / 4.0 / no_of_trials;
		put_option_price_2 += (p1_2 + p2_2 + p3_2 + p4_2) / 4.0 / no_of_trials;
	}
	call_option_price = exp(-risk_free_rate * expiration_time) * call_option_price;
	put_option_price = exp(-risk_free_rate * expiration_time) * put_option_price ;

	call_option_price_2 = exp(-risk_free_rate * expiration_time) * call_option_price_2;
	put_option_price_2 = exp(-risk_free_rate * expiration_time) * put_option_price_2;

	cout << "Call Price by simulation = " << call_option_price << endl;
	cout << "Put Price by simulation = " << put_option_price << endl;

	cout << "Call Price by Adjusted Method = " << call_option_price_2 << endl;
	cout << "Put Price by Adjusted Method = " << put_option_price_2 << endl;
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a * p);
	double b = c2 * exp((-z) * (z / 2.0));
	double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
};

double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time)
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return K * exp(-r * time) * N(-d2) - S * N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time)	  // time to maturity 
{
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
	double d2 = d1 - (sigma * time_sqrt);
	return S * N(d1) - K * exp(-r * time) * N(d2);
};

float closed_form_down_and_out_european_call_option()
{
	// I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
	float K = (2 * risk_free_rate) / (volatility * volatility);
	float A = option_price_call_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	float B = (barrier_price * barrier_price) / initial_stock_price;
	float C = pow(initial_stock_price / barrier_price, -(K - 1));
	float D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
	return (A - D * C);
}

float closed_form_down_and_in_european_put_option()
{
	// just making it easier by renaming the global variables locally
	float S = initial_stock_price;
	float r = risk_free_rate;
	float T = expiration_time;
	float sigma = volatility;
	float H = barrier_price;
	float X = strike_price;

	// Took these formulae from some online reference
	float lambda = (r + ((sigma * sigma) / 2)) / (sigma * sigma);
	float temp = 2 * lambda - 2.0;
	float x1 = (log(S / H) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
	float y = (log(H * H / (S * X)) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
	float y1 = (log(H / S) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
	return (-S * N(-x1) + X * exp(-r * T) * N(-x1 + sigma * sqrt(T)) +
		S * pow(H / S, 2 * lambda) * (N(y) - N(y1)) -
		X * exp(-r * T) * pow(H / S, temp) * (N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

float closed_form_down_and_out_european_put_option()
{
	float vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
		risk_free_rate, volatility, expiration_time);
	float put_down_in = closed_form_down_and_in_european_put_option();
	return (vanilla_put - put_down_in);
}

int main(int argc, char* argv[])
{

	sscanf(argv[1], "%lf", &expiration_time);
	sscanf(argv[2], "%lf", &risk_free_rate);
	sscanf(argv[3], "%lf", &volatility);
	sscanf(argv[4], "%lf", &initial_stock_price);
	sscanf(argv[5], "%lf", &strike_price);
	sscanf(argv[6], "%d", &no_of_trials);
	sscanf(argv[7], "%d", &no_of_divisions);
	sscanf(argv[8], "%lf", &barrier_price);

	up_factor = exp(volatility * sqrt(expiration_time / ((double)no_of_divisions)));
	R = exp(risk_free_rate * expiration_time / ((double)no_of_divisions));
	uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
	dt = expiration_time / no_of_divisions;

	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials (4 paths each trial) = " << no_of_trials << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << dt << endl;
	cout << "--------------------------------------" << endl;
	cout << "Price of an European Down and Out Call Option from Theory = " <<
		closed_form_down_and_out_european_call_option() << endl;
	cout << "Price of an European Down and Out Put Option from Theory = " <<
		closed_form_down_and_out_european_put_option() << endl;
	cout << "--------------------------------------" << endl;
	monte_carlo_simulation(no_of_divisions, no_of_trials);
}