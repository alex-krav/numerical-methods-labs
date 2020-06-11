#include "io_func.h"
#include "matrix_func.h"
#include "math_func.h"
#include <map>

using namespace std;

double find_integral_trapeze(string, double (*func)(double), vector<double>, double);
double find_integral_gauss(string, double (*func)(double), vector<double>, double);
int find_n(string, double (*func)(double), int, int, double);
int find_m(string, double (*func)(double), int, int, double);
vector<vector<double>> get_wx_coef(int);

int main6()
{
	string in = "input_lab6.txt", out = "output_lab6.txt";
	vector<double> intervals;
	double trapeze, gauss, eps = 0.0001;

	read_vec_from_file(in, out, &intervals);
	find_integral_trapeze(out, cos_func, intervals, eps);
	find_integral_gauss(out, cos_func, intervals, eps);

	return 0;
}

double find_integral_trapeze(string out, double (*func)(double), vector<double> intervals, double eps) {
	double sum = 0, a = intervals[0], b = intervals[1], fa, fb;
	int n = find_n(out, cos_func_der2, a, b, eps); n = 100;
	double h = (b - a) / n;
	double x = a + h;
	char mes[100];
	
	while (x < b) {
		sum += func(x);
		x += h;
	}

	fa = func(a);
	fb = func(b);
	sum = (h / 2) * (fa + fb + 2 * sum);

	sprintf(mes, "\ntrapeze integral: %0.6f", sum);
	cout << mes << endl;
	save_str_to_file(out, mes);

	return sum;
}

double find_integral_gauss(string out, double (*func)(double), vector<double> intervals, double eps) {
	int m = find_m(out, cos_func, intervals[0], intervals[1], eps); m = 4;
	vector<vector<double>> wx_coef = get_wx_coef(m);
	int size = wx_coef.end() - wx_coef.begin();
	double sum = 0, x;
	char mes[100];

	for (int i = 0; i < size; ++i) {
		x = (intervals[1] + intervals[0]) / 2 
			+ ((intervals[1] - intervals[0]) / 2) * wx_coef[i][1];
		sum += wx_coef[i][0] * func(x);
	}
	sum = sum * (intervals[1] - intervals[0]) / 2;

	sprintf(mes, "\ngauss integral: %0.6f", sum);
	cout << mes << endl;
	save_str_to_file(out, mes);

	return sum;
}

int find_n(string out, double (*func)(double), int a, int b, double eps) {
	int n = 2;
	double error; 
	char mes[100];

	while (1) {
		error = (pow((b - a), 3) / (12 * pow(n, 2))) 
			* find_min_max_abs(func, a, b)[1];
		if (error <= eps)
			break;
		else
			++n;
	}
	
	sprintf(mes, "\ntrapeze error: %0.6f", error);
	cout << mes << endl; 
	save_str_to_file(out, mes);

	return n;
}

int find_m(string out, double (*func)(double), int a, int b, double eps) {
	int m = 2;
	double error = 1;
	char mes[100];

	while (1) {
		error = ((pow(factorial(m),4)*pow((b-a),(m-1)))
			/ ((2*m+1)*pow(factorial(2*m),3)))
			* find_min_max_abs_pow(func, a, b, 2*m)[1];
		if (error <= eps)
			break;
		else
			++m;
	}

	sprintf(mes, "\ngauss error: %0.6f", error);
	cout << mes << endl;;
	save_str_to_file(out, mes);

	return m;
}

vector<vector<double>> get_wx_coef(int m) {
	map<int, vector<vector<double>>> wx_coef;
	vector<vector<double>> vec;

	vec.push_back(vector<double>{1, -0.577350269});
	vec.push_back(vector<double>{1, 0.577350269});
	wx_coef[2] = vec; vec.clear();

	vec.push_back(vector<double>{0.555555556, -0.774596669});
	vec.push_back(vector<double>{0.888888889, 0});
	vec.push_back(vector<double>{0.555555556, 0.774596669});
	wx_coef[3] = vec; vec.clear();

	vec.push_back(vector<double>{0.347854845, -0.861136312});
	vec.push_back(vector<double>{0.652145155, -0.339981044});
	vec.push_back(vector<double>{0.652145155, 0.339981044});
	vec.push_back(vector<double>{0.347854845, 0.861136312});
	wx_coef[4] = vec; vec.clear();

	vec.push_back(vector<double>{0.236926885, -0.906179846});
	vec.push_back(vector<double>{0.478628670, -0.538469310});
	vec.push_back(vector<double>{0.568888889, 0});
	vec.push_back(vector<double>{0.478628670, 0.538469310});
	vec.push_back(vector<double>{0.236926885, 0.906179846});
	wx_coef[5] = vec; vec.clear();

	vec.push_back(vector<double>{0.171324492, -0.932469514});
	vec.push_back(vector<double>{0.360761573, -0.661209386});
	vec.push_back(vector<double>{0.467913935, -0.238619186});
	vec.push_back(vector<double>{0.467913935, 0.238619186});
	vec.push_back(vector<double>{0.360761573, 0.661209386});
	vec.push_back(vector<double>{0.171324492, 0.932469514});
	wx_coef[6] = vec; vec.clear();

	return wx_coef.at(m);
}
