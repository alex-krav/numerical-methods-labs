#include "math_func.h"

using namespace std;

double cos_func(double x) {
	return cos(x) / (x + 1);
}

double cos_func_der1(double x) {
	return -sin(x) / (x + 1) - cos(x) / pow((x + 1), 2);
}

double cos_func_der2(double x) {
	return -cos(x) / (x + 1) + 2 * sin(x) / pow((x + 1), 2) + 2 * cos(x) / pow((x + 1), 3);
}

vector<double> find_min_max(double (*func)(double), double a, double b) {
	double x, f, step = (b - a) / 1000;
	vector<double> min_max(2, 0);

	for (x = a; x <= b - step; x += step) {
		f = func(x);
		if (f < min_max[0])
			min_max[0] = f;
		if (f > min_max[1])
			min_max[1] = f;
	}

	return min_max;
}

vector<double> find_min_max_abs(double (*func)(double), double a, double b) {
	double x, f, step = (b - a) / 1000;
	vector<double> min_max(2, 0);

	for (x = a; x <= b - step; x += step) {
		f = abs(func(x));
		if (f < min_max[0])
			min_max[0] = f;
		if (f > min_max[1])
			min_max[1] = f;
	}

	return min_max;
}

vector<double> find_min_max_abs_pow(double (*func)(double), double a, double b, int power) {
	double x, f, step = (b - a) / 1000;
	vector<double> min_max(2, 0);

	for (x = a; x <= b - step; x += step) {
		f = abs(pow(func(x), power));
		if (f < min_max[0])
			min_max[0] = f;
		if (f > min_max[1])
			min_max[1] = f;
	}

	return min_max;
}

int factorial(int n) {
	return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}