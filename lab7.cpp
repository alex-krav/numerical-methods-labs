#include "io_func.h"
#include "matrix_func.h"
#include "math_func.h"

using namespace std;

vector<vector<double>> runge_cutta_method(string out, double (*func)(double, double), double a, double b, double h, double y0, int n);
vector<vector<double>> adams_method(string out, double (*func)(double, double), double a, double b, double h, double y0);
double ode_func(double x, double y);

int main7() {
	// input: a, b, h, y0 - 0 1 0.1 0
	string in = "input_lab7.txt", out = "output_lab7.txt";
	vector<double> input;
	vector<vector<double>> roots;

	read_vec_from_file(in, out, &input);
	roots = adams_method(out, ode_func, input[0], input[1], input[2], input[3]);
	print_vec_to_cout("x,y:", roots);
	save_vec_to_file(out, "x,y:", roots);

	return 0;
}

vector<vector<double>> runge_cutta_method(string out, double (*func)(double, double), double a, double b, double h, double y0, int n) {
	vector<vector<double>> coef;
	vector<double> row{ 2,0 };
	coef.push_back(row);
	coef[0][0] = a;
	coef[0][1] = y0;

	double k1=1, k2=1, k3=1, k4=1, dy, tau, eps, xn = a+h;
	int i = 0;
	char iter_mes[100];

	string mes = "\nrunge-cutta method:";
	cout << mes << endl;
	save_str_to_file(out, mes);
	while (i < n) {
		++i;
		row.clear(); row.push_back(0); row.push_back(0); coef.push_back(row);
		coef[i][0] = coef[i - 1][0] + h;
		k1 = h * func(coef[i][0], coef[i][1]);
		k2 = h * func(coef[i][0] + h/2, coef[i][1] + k1/2);
		k3 = h * func(coef[i][0] + h/2, coef[i][1] + k2/2);
		k4 = h * func(coef[i][0] + h, coef[i][1] + k3);
		dy = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		coef[i][1] = coef[i - 1][1] + dy;
		xn += h;

		tau = abs((k2-k3) / (k1-k2));
		eps = (pow(coef[i][1], h) - pow(coef[i][1], h / 2)) / (pow(2, 4) - 1);
		sprintf(iter_mes, "eps=%f, tau=%f", eps, tau);
		cout << iter_mes << endl;
		save_str_to_file(out, iter_mes);
	}
	
	return coef;
}

vector<vector<double>> adams_method(string out, double (*func)(double, double), double a, double b, double h, double y0) {
	vector<vector<double>> coef = runge_cutta_method(out, func, a, b, h, y0, 4);
	print_vec_to_cout("x,y:", coef);
	vector<double> row;
	double f1, f2, f3, f4, eps, xn = coef[4][1];
	int i = 4;
	char iter_mes[100];

	string mes = "\nadams method:";
	cout << mes << endl;
	save_str_to_file(out, mes);
	while (xn <= b) {
		++i;
		row.clear(); row.push_back(0); row.push_back(0); coef.push_back(row);
		f1 = func(coef[i - 1][0], coef[i - 1][1]);
		f2 = func(coef[i - 2][0], coef[i - 2][1]);
		f3 = func(coef[i - 3][0], coef[i - 3][1]);
		f4 = func(coef[i - 4][0], coef[i - 4][1]);
		coef[i][0] = coef[i - 1][0] + h;
		coef[i][1] = coef[i - 1][1] + h / 24 * (55*f1 - 59*f2 + 37*f3 - 9*f4);
		xn += h;

		eps = (pow(coef[i][1], h) - pow(coef[i][1], h/2)) / (pow(2,4)-1);
		sprintf(iter_mes, "eps=%f", eps);
		cout << iter_mes << endl;
		save_str_to_file(out, iter_mes);
	}
	
	return coef;
}

double ode_func(double x, double y) {
	return 1 + 1.8 * y * sin(x) - pow(y, 2);
}
