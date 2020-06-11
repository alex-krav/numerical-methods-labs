#ifndef MATH_FUNCTIONS_H_INCLUDED
#define MATH_FUNCTIONS_H_INCLUDED

#include <cmath>
#include <vector>

using namespace std;

int factorial(int);
double cos_func(double);
double cos_func_der1(double);
double cos_func_der2(double);
vector<double> find_min_max(double (*func)(double), double, double);
vector<double> find_min_max_abs(double (*func)(double), double, double);
vector<double> find_min_max_abs_pow(double (*func)(double), double, double, int);

#endif
