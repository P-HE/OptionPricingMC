//Modified on the code example provided by Prof.RS
//Takes significantly longer time than the sample screenshot
//And generates slightly different number(with error < 0.005)

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <vector>
using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price;
float initial_stock_price, expiration_time, volatility, R, barrier_price;
int no_of_divisions = 500;
int no_of_sampling_instants;
int no_of_trials;
double dt;

float max(float a, float b) {
	return (b < a) ? a : b;
}

double get_uniform() {
	return (double)rand() / (double)RAND_MAX;
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

		bool o1 = 1;
		bool o2 = 1;
		bool o3 = 1;
		bool o4 = 1;

		int count = 0;
		int count1 = 0;

		for (int j = 0; j < no_of_divisions; j++) {
			double x = get_uniform();
			double y = get_uniform();

			while (x * y == 0) {
				x = get_uniform();
				y = get_uniform();
			}

			double a = sqrt(-2.0 * log(x)) * cos(6.283185307999998 * y);
			double b = sqrt(-2.0 * log(x)) * sin(6.283185307999998 * y);

			S1 = S1 * exp(R + SD * a);
			S2 = S2 * exp(R - SD * a);
			S3 = S1 * exp(R + SD * b);
			S4 = S4 * exp(R - SD * b);


			if (count % (no_of_divisions / no_of_sampling_instants) == 0) {
				if (S1 <= barrier_price) {
					o1 = 0;
				}
				if (S2 <= barrier_price) {
					o2 = 0;
				}
				if (S3 <= barrier_price) {
					o3 = 0;
				}
				if (S4 <= barrier_price) {
					o4 = 0;
				}
				//cout << o1 << o2 << o3 << o4 << endl;
			}
			count++;
			//cout << x << y << endl;
			//cout << a << b << endl;
		}
		//cout << S1 <<" "<< S2 <<" "<< S3 << " " << S4 << endl;

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

		call_option_price += (p1 + p2 + p3 + p4) / 4.0 ;

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

		//cout << put_option_price << endl;
		//put_option_price += (c1 + c2 + c3 + c4) / 4.0;
		put_option_price += (c1 + c2 + c3 + c4) / 4.0;

		double probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant;
		double c1_2, c2_2, c3_2, c4_2;
		double p1_2, p2_2, p3_2, p4_2;

		vector<double> mean_at_sampling_instant;
		vector<double> variance_at_sampling_instant;
		mean_at_sampling_instant.push_back(0.0);
		variance_at_sampling_instant.push_back(0.0);


		for (int j = 1; j <= no_of_sampling_instants; j++) {
			double mean = initial_stock_price +
				(((float)j) / ((float)no_of_sampling_instants) * (S1 - initial_stock_price));
			mean_at_sampling_instant.push_back(mean);
			double variance = (((float)j) / ((float)no_of_sampling_instants)) * expiration_time *
				(1.0 - ((float)j) / ((float)no_of_sampling_instants));
			variance_at_sampling_instant.push_back(variance);
		}

		probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
		for (int j = 1; j <= no_of_sampling_instants; j++)
			probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
			(1.0 - N((barrier_price - mean_at_sampling_instant[j]) / sqrt(variance_at_sampling_instant[j])));

		c1_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
			max(0.0, (S1 - strike_price)));
		if (S1 > barrier_price) {
			p1_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
				max(0.0, (strike_price - S1)));
		}
		else {
			p1_2 = 0.0;
		}


		for (int j = 1; j <= no_of_sampling_instants; j++) {
			double mean = initial_stock_price +
				(((float)j) / ((float)no_of_sampling_instants) * (S2 - initial_stock_price));
			mean_at_sampling_instant.push_back(mean);
			double variance = (((float)j) / ((float)no_of_sampling_instants)) * expiration_time *
				(1.0 - ((float)j) / ((float)no_of_sampling_instants));
			variance_at_sampling_instant.push_back(variance);
		}

		probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
		for (int j = 1; j <= no_of_sampling_instants; j++)
			probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
			(1.0 - N((barrier_price - mean_at_sampling_instant[j]) / sqrt(variance_at_sampling_instant[j])));

		c2_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
			max(0.0, (S2 - strike_price)));
		if (S2 > barrier_price) {
			p2_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
				max(0.0, (strike_price - S2)));
		}
		else {
			p2_2 = 0.0;
		}


		for (int j = 1; j <= no_of_sampling_instants; j++) {
			double mean = initial_stock_price +
				(((float)j) / ((float)no_of_sampling_instants) * (S3 - initial_stock_price));
			mean_at_sampling_instant.push_back(mean);
			double variance = (((float)j) / ((float)no_of_sampling_instants)) * expiration_time *
				(1.0 - ((float)j) / ((float)no_of_sampling_instants));
			variance_at_sampling_instant.push_back(variance);
		}

		probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
		for (int j = 1; j <= no_of_sampling_instants; j++)
			probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
			(1.0 - N((barrier_price - mean_at_sampling_instant[j]) / sqrt(variance_at_sampling_instant[j])));


		c3_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
			max(0.0, (S3 - strike_price)));
		if (S3 > barrier_price) {
			p3_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
				max(0.0, (strike_price - S3)));
		}
		else {
			p3_2 = 0.0;
		}


		for (int j = 1; j <= no_of_sampling_instants; j++) {
			double mean = initial_stock_price +
				(((float)j) / ((float)no_of_sampling_instants) * (S4 - initial_stock_price));
			mean_at_sampling_instant.push_back(mean);
			double variance = (((float)j) / ((float)no_of_sampling_instants)) * expiration_time *
				(1.0 - ((float)j) / ((float)no_of_sampling_instants));
			variance_at_sampling_instant.push_back(variance);
		}

		probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant = 1.0;
		for (int j = 1; j <= no_of_sampling_instants; j++)
			probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *=
			(1.0 - N((barrier_price - mean_at_sampling_instant[j]) / sqrt(variance_at_sampling_instant[j])));

		c4_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
			max(0.0, (S4 - strike_price)));
		if (S4 > barrier_price) {
			p4_2 = (probability_stock_path_has_not_hit_barrier_at_some_discrete_sampling_instant *
				max(0.0, (strike_price - S4)));
		}
		else {
			p4_2 = 0.0;
		}

		call_option_price_2 += (c1_2 + c2_2 + c3_2 + c4_2) / 4.0 / no_of_trials;
		put_option_price_2 += (p1_2 + p2_2 + p3_2 + p4_2) / 4.0 / no_of_trials;
	}

	call_option_price = exp(-risk_free_rate * expiration_time) * call_option_price / no_of_trials;
	put_option_price = exp(-risk_free_rate * expiration_time) * put_option_price / no_of_trials;

	call_option_price_2 = exp(-risk_free_rate * expiration_time) * call_option_price_2;
	put_option_price_2 = exp(-risk_free_rate * expiration_time) * put_option_price_2;

	cout << "Call Price by simulation = " << call_option_price << endl;
	cout << "Put Price by simulation = " << put_option_price << endl;

	cout << "Call Price by Brownian-Bridge Adjusted Method = " << call_option_price_2 << endl;
	cout << "Put Price by Brownian-Bridge Adjusted Method = " << put_option_price_2 << endl;
}

int main(int argc, char* argv[])
{

	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%f", &risk_free_rate);
	sscanf(argv[3], "%f", &volatility);
	sscanf(argv[4], "%f", &initial_stock_price);
	sscanf(argv[5], "%f", &strike_price);
	sscanf(argv[6], "%d", &no_of_trials);
	sscanf(argv[7], "%d", &no_of_sampling_instants);
	sscanf(argv[8], "%f", &barrier_price);

	cout << "European Down and Out Discrete Barrier Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of (uniformly spaced) Discrete Barrier-Samples = " << no_of_sampling_instants << endl;
	cout << "Number of Trails (Each Trial has 4 Price Paths) = " << no_of_trials << endl;
	dt = expiration_time / no_of_divisions;
	cout << dt << endl;
	cout << "--------------------------------------" << endl;
	monte_carlo_simulation(no_of_divisions, no_of_trials);
	//monte_carlo_simulation(no_of_divisions, 2);
}