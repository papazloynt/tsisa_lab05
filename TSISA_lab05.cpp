//Variant 16
#include <random>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

struct neuron {
	double x;
	double y;
};

double lin_func(const double c, const double d, const double x) {
	return  c * x + d;
}

double random(const double a, const double b) {
	return a + rand() * 1. / RAND_MAX * (b - a);
}

double rand_er(const double a, const double b, const double A) {
	return A * random(a, b);
}

std::vector<neuron> fill(const double c, const double d, const double A, const double a, const double b, const size_t N) {
	double delta = (b - a) / double(N + 1);
	std::vector<neuron> n_(N);
	double current_x = a;
	for (auto& n : n_) {
		current_x += delta;
		n.x = current_x;
		n.y = lin_func(c, d, n.x) + rand_er(-0.5, 0.5, A);
	}
	return n_;
}


double sum_of_squares(const double c, const double d, std::vector<neuron> n_) {
	double s = 0.0;
	for (auto& n : n_) {
		double x = n.x;
		double t = n.y;
		double	y = c * x + d;
		s += (y - t) * (y - t);
		return s;
	}
}




double comp(const std::vector<neuron>& n_, const double c, const double d) {
	double sum = 0.;
	for (auto n : n_) {
		sum += pow(n.y - (c * n.x + d), 2);
	}
	return sum;
}


double random_search(const std::vector<neuron>& n, double c) {
	const double D_MIN = 4;
	const double D_MAX = 6;
	const size_t iterations = 50;

	double d = D_MIN;
	for (size_t i = 0; i < iterations; ++i) {
		double new_d = random(D_MIN, D_MAX);
		if (comp(n, c, new_d) < comp(n, c, d)) {
			d = new_d;
		}
	}
	return d;
}




double dichotomy(double c_min, double c_max, std::vector<neuron>& n) {
	double eps = 0.01;
	double delta = 0.001;
	while ((c_max - c_min) > eps) {
		double c1 = 0.5 * (c_min + c_max) - delta;
		double c2 = 0.5 * (c_min + c_max) + delta;
		double d1 = random_search(n, c1);
		double d2 = random_search(n, c2);
		if (sum_of_squares(c2, d2, n) > sum_of_squares(c1, d1, n)) {
			c_max = c2;
		}
		else {
			c_min = c1;
		}
	}
	return (c_min + c_max) / 2;
}


void print(const std::vector<neuron>& n_) {
	for (const auto& n : n_) {
		std::cout << "| " << std::setw(10) << n.x
			<< " | " << std::setw(10) << n.y << " |" << std::endl;
	}
}

int main() {
	const double a = 0;
	const double b = 4;
	const double c = -5;
	const double d = 5;
	const size_t N = 16;
	const double A = 2;

	std::vector<neuron> n = fill(c, d, A, a, b, N);

	double res_1 = dichotomy(-6, -4, n);

	double res_2 = random_search(n, res_1);

	print(n);

	std::cout << "Result :" << std::endl;
	std::cout << "w1 = " << res_1 << "  ;  " << "w0 = " << res_2 << std::endl;

	return 0;
}