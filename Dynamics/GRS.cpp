#include <iostream>
#include <random>
#include <chrono>
#include <cstring>
#include "fitting.h"

static int my_rand_disc(int n) {
	static long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::default_random_engine generator(seed);
	static std::uniform_int_distribution<int> distribution(0, n - 1);
	return distribution(generator);
}

static double my_rand() {
	static long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::default_random_engine generator(seed);
	static std::uniform_real_distribution<double> distribution(-0.5, 0.5);
	return distribution(generator);
}

static double parameter_tweak(double x, double g) {
	double t = my_rand() * g * x;
	return x + t;
}

void GRS(double fn(double x[], CalcData &calcData), int n, double start[], double xmin[], CalcData &calcData) {
	//memcpy(xmin, start, sizeof(double) * n);
	for (int i = 0; i < n; i++) {
		xmin[i] = start[i] * 2;
	}
	#pragma omp parallel
	{
		CalcData calcDataLocal = calcData;
		int iteration = 0;
		int count = 0;
		int successes = 0;
		double g = 0.1;
		double *x_curr = new double[n];
		memcpy(x_curr, start, sizeof(double) * n);
		double *x_new = new double[n];
		double m_curr = fn(x_curr, calcDataLocal);
		while (iteration < 100 && g > 0.000001) {
			int choice = my_rand_disc(n);
			double v_new = parameter_tweak(x_curr[choice], g);
			memcpy(x_new, x_curr, sizeof(double) * n);
			x_new[choice] = v_new;
			double m_new = fn(x_new, calcDataLocal);
			if (m_new < m_curr) {
				double * tmp = x_new;
				x_new = x_curr;
				x_curr = tmp;
				m_curr = m_new;
				++successes;
			}
			if (count == 100) {
				count = 0;
				if (successes < 5) {
					g /= 10;
				}
				successes = 0;
			}
			//std::cout << iteration << " " << g << " " << x_curr[0] << " " << x_curr[1] << " " << m_curr << "\n";
			++iteration;
			++count;
		}
		#pragma omp critical
		{
			if (fn(x_curr, calcDataLocal) < fn(xmin, calcDataLocal)) {
				memcpy(xmin, x_curr, sizeof(double) * n);
			}
		}
		//std::cout << iteration << "\n";
	}
}
