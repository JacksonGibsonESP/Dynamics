#include <iostream>
#include <random>
#include <chrono>
using namespace std;

static void G(double * x, double * e, int n) {
	static long long int seed = chrono::system_clock::now().time_since_epoch().count();
	static default_random_engine generator(seed);
	static vector<uniform_real_distribution<double>*> distribution;
	static bool firstStart = true;
	if (firstStart) {
		firstStart = false;
		for (int i = 0; i < n; i++) {
			distribution.push_back(new uniform_real_distribution<double>(0, 2 * e[i]));
		}
	}
	for (int i = 0; i < n; i++) {
		x[i] = (*distribution[i])(generator);
	}
}

static double my_rand() {
	static long long int seed = chrono::system_clock::now().time_since_epoch().count();
	static default_random_engine generator(seed);
	static uniform_real_distribution<double> distribution(0, 1);
	return distribution(generator);
}


static double T(int i) {
	return exp(8 - i * 0.05);
}


void annealing(double fn(double x[]), int n, double start[], double xmin[]) {
	int iteration = 0;
	double t = T(iteration);
	double * x = new double[n];
	double * temp = new double[n];
	G(x, start, n);
	while (t > 0.01) {
		G(temp, start, n);
		double dE = fn(temp) - fn(x);
		if (dE <= 0) {
			memcpy(x, temp, n * sizeof(double));
		}
		else {
			double p = exp(-dE / t);
			if (my_rand() <= p) {
				memcpy(x, temp, n * sizeof(double));
			}
		}
		//		for (int i = 0; i < n; i++) {
		//			std::cout << x[i] << " ";
		//		}
		//		std::cout << fn(x) << " " << t << "\n";
		++iteration;
		t = T(iteration);
	}
	memcpy(xmin, x, n * sizeof(double));
	//std::cout << iteration << "\n";
}
