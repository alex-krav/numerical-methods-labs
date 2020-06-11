#include "io_func.h"
#include "matrix_func.h"
#include <cmath>

using namespace std;

vector<double> find_roots(string, string, double(*func)(vector<double>, double, double, double,int&),
	vector<double>, vector<vector<double>>, double);
double get_root_dich(vector<double>, double, double, double, int&);
double get_root_chord(vector<double>, double, double, double, int&);
double get_root_tan(vector<double>, double, double, double, int&);
double get_func_val(vector<double>, double);
double get_func_val(double(*func)(double), double);
vector<double> derivative1(vector<double>);
vector<double> derivative2(vector<double>);

int main5()
{
	string in_intervals = "input_lab5_intervals.txt",  // 0.263 1.3
		in_polynomials = "input_lab5_polynomials.txt", 
		out = "output_lab5.txt";
	vector<vector<double>> intervals;
	vector<double> polyVec;

	double eps = 0.0001;
	read_vec_from_file(in_intervals, out, &intervals);
	read_vec_from_file(in_polynomials, out, &polyVec);

	find_roots(out, "dichotomy roots:", get_root_dich, polyVec, intervals, eps);
	find_roots(out, "chord roots:", get_root_chord, polyVec, intervals, eps);
	find_roots(out, "tangent roots:", get_root_tan, polyVec, intervals, eps);

	return 0;
}

vector<double> find_roots(string out, string mes, 
	double(*func)(vector<double>,double,double,double,int&), 
	vector<double> polyVec, vector<vector<double>> intervals, double eps) {
	vector<double> roots;
	int size = intervals.end() - intervals.begin(), counter = 0;

	for (int i = 0; i < size; ++i) {
		roots.push_back(func(polyVec, intervals[i][0], intervals[i][1], eps, counter));
	}
	print_vec_to_cout(mes, roots);
	save_vec_to_file(out, mes, roots);

	char iter_str[100];
	sprintf(iter_str, "%d iterations", counter);
	cout << iter_str << endl;
	save_str_to_file(out, iter_str);

	return roots;
}


double get_root_tan(vector<double> polyVec, double a, double b, double eps, int& counter) {
	double x0, xp, xn, r=1;
	double der1 = get_func_val(derivative1(polyVec), a);
	double der2 = get_func_val(derivative2(polyVec), a);
	
	if (der1 * der2 >= 0)
		x0 = b;
	else 
		x0 = a;
	xp = x0;

	while (abs(r) > eps) {
		xn = xp - get_func_val(polyVec, xp) / get_func_val(derivative1(polyVec), xp);
		r = xn - xp;
		xp = xn;
		++counter;
	}

	return xp;
}


double get_root_chord(vector<double> polyVec, double a, double b, double eps, int& counter) {
	double x0, c, xp, xn, r=1;
	double der1 = get_func_val(derivative1(polyVec), a);
	double der2 = get_func_val(derivative2(polyVec), a);
	
	if (der1 * der2 < 0) {
		x0 = b;
		c = a;
	}
	else {
		x0 = a;
		c = b;
	}
	xp = x0;

	while (abs(r) > eps) {
		xn = xp - get_func_val(polyVec, xp) * (c - xp) / (get_func_val(polyVec, c) - get_func_val(polyVec, xp));
		r = xn - xp;
		xp = xn;
		++counter;
	}

	return xp;
}

double get_root_dich(vector<double> polyVec, double a, double b, double eps, int& counter) {
	double m;

	while (abs(b - a) > eps) {
		m = (a + b) / 2;
		if (get_func_val(polyVec, m) < 0)
			a = m;
		else
			b = m;
		++counter;
	}

	return (a + b) / 2;
}

double get_func_val(vector<double> polyVec, double x) {
	double val = 0;
	int size = polyVec.end() - polyVec.begin();

	for (int i = 0; i < size; ++i) {
		val += polyVec[i] * pow(x, i);
	}

	return val;
}

double get_func_val(double(*func)(double), double x) {
	return func(x);
}

vector<double> derivative1(vector<double> polyVec) {
	vector<double> der;
	int size = polyVec.end() - polyVec.begin();

	for (int i = 1; i < size; ++i) {
		der.push_back(polyVec[i] * i);
	}

	return der;
}

vector<double> derivative2(vector<double> polyVec) {
	vector<double> der;
	int size = polyVec.end() - polyVec.begin();

	for (int i = 2; i < size; ++i) {
		der.push_back(polyVec[i] * i);
	}

	return der;
}