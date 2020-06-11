#include "io_func.h"
#include "matrix_func.h"
#include <cmath>

using namespace std;

double myFunc(double);
double polyFunc(double);

double calc_diff(vector<double>, int, int&, double (*func)(double), vector<double>*);
vector<double> find_newton_polynom_coef(string, vector<double>, double (*func)(double));
string build_newton_polynom_str(string, vector<double>, vector<double>);

vector<vector<double>> find_spline_coef(string, vector<double>, vector<double>);
vector<double> find_spline_a_coef(vector<double>);
vector<double> find_spline_c_coef(vector<double>, vector<double>);
vector<double> find_spline_d_coef(vector<double>, vector<double>);
vector<double> find_spline_b_coef(vector<double>, vector<double>, vector<double>);
vector<double> progonka(vector<vector<double>>, vector<double>);

vector<double> find_error_coef(string, vector<double>, vector<double>, vector<double>);
vector<double> calc_polynom(string, vector<double>, vector<double>);

int main()
{
	string in = "input_lab4.txt", out = "output_lab4.txt";
	vector<double> xVec, yVec, poly_coef;
	read_vec_from_file(in, out, &xVec);
	
	poly_coef = find_newton_polynom_coef(out, xVec, myFunc);
	build_newton_polynom_str(out, xVec, poly_coef);

	yVec = calc_values(xVec, myFunc);
	find_spline_coef(out, xVec, yVec);

	print_vec_to_cout("y values:", yVec);
	save_vec_to_file(out, "y values:", yVec);
	find_error_coef(out, xVec, yVec, poly_coef);

	return 0;
}

string get_x_str(double val) {
	char x_str[100];
	if (val > 0)
		sprintf(x_str, " + (x - %.0f)", val);
	else if (val < 0)
		sprintf(x_str, " + (x + %.0f)", val*(-1));
	else
		sprintf(x_str, " + x");
	return x_str;
}

string get_y_str(double val) {
	char y_str[100];
	if (val == 0)
		sprintf(y_str, "0");
	else
		sprintf(y_str, "%f", val);
	return y_str;
}

string build_newton_polynom_str(string out, vector<double> x, vector<double> y) {
	string str = get_y_str(y[0]) + get_x_str(x[0]);
	int size = x.end() - x.begin();

	for (int i = 1; i < size; ++i) {
		str += " * [" + get_y_str(y[i]) + get_x_str(x[i]);
	}
	for (int i = 1; i < size; ++i) {
		str += "]";
	}

	cout << "\nnewton polynom equation:" << endl;
	cout << str << endl;
	save_str_to_file(out, "\nnewton polynom equation:\n");
	save_str_to_file(out, str + "\n");
	 
	return str;
}

vector<double> find_newton_polynom_coef(string out, vector<double> x, double (*func)(double)) {
	vector<double> coef;
	int counter = 1;
	calc_diff(x, (x.end()-x.begin()), counter, myFunc, &coef);

	print_vec_to_cout("newton polynom coef:", coef);
	save_vec_to_file(out, "newton polynom coef:", coef);
	return coef;
}

vector<double> find_error_coef(string out, vector<double> x, vector<double> y, vector<double> poly_coef) {
	vector<double> error_coef, polynoms = calc_polynom(out, x, poly_coef);
	int size = x.end() - x.begin();
	for (int i = 0; i < size; ++i) {
		error_coef.push_back(abs(polynoms[i]-y[i]));
	}
	print_vec_to_cout("error coef:", error_coef);
	save_vec_to_file(out, "error coef:", error_coef);
	return error_coef;
}

vector<double> calc_polynom(string out, vector<double> xVec, vector<double> poly_coef) {
	vector<double> polynoms;
	double mult_val, cur_val;
	int size = xVec.end() - xVec.begin();

	for (int j = 0; j < size; ++j) {
		mult_val = 1;
		for (int i = size - 1; i > -1; --i) {
			cur_val = (xVec[j] - xVec[i]) * mult_val + poly_coef[i];
			mult_val = cur_val;
		}
		polynoms.push_back(mult_val);
	}
	print_vec_to_cout("polynom values:", polynoms);
	save_vec_to_file(out, "polynom values:", polynoms);

	return polynoms;
}

double calc_diff(vector<double> x, int n, int& counter, double (*func)(double), vector<double> *coef) {
	double val;
	if (n == 1) {
		val =  func(*x.begin());
		if (counter == 1)
			coef->push_back(val);
		++counter;
		return val;
	}
	else {
		double val = (calc_diff({ x.begin(), x.end() - 1 }, n - 1, counter, func, coef)
			- calc_diff({x.begin()+1, x.end()}, n-1, counter, func, coef))
			/ (x.front() - x.back());
		++counter;
		if (pow(2,n) == counter)
			coef->push_back(val);
		return val;
	}
}

double myFunc(double x) {
	return sin(0.5 * x) + cbrt(x);
}

double polyFunc(double x) {
	return -2.496698 + (x + 4) * (0.197653 + (x + 2) * (0.213261 + x * (-0.035543 + (x - 2) * (x - 4))));
}

vector<vector<double>> find_spline_coef(string out, vector<double> x, vector<double> y) {
	vector<vector<double>> spline_coef;

	vector<double> aCoef = find_spline_a_coef(y);
	spline_coef.push_back(aCoef);
	print_vec_to_cout("a coef:", aCoef);
	save_vec_to_file(out, "a coef:", aCoef);

	vector<double> cCoef = find_spline_c_coef(x, y);
	spline_coef.push_back(cCoef);
	print_vec_to_cout("c coef:", cCoef);
	save_vec_to_file(out, "c coef:", cCoef);

	vector<double> dCoef = find_spline_d_coef(x, cCoef);
	spline_coef.push_back(dCoef);
	print_vec_to_cout("d coef:", dCoef);
	save_vec_to_file(out, "d coef:", dCoef);

	vector<double> bCoef = find_spline_b_coef(x, y, cCoef);
	spline_coef.push_back(bCoef);
	print_vec_to_cout("b coef:", bCoef);
	save_vec_to_file(out, "b coef:", bCoef);

	return spline_coef;
}

vector<double> find_spline_a_coef(vector<double> yVec) {
	vector<double> aCoef;
	int size = yVec.end() - yVec.begin();
	for (int i = 1; i < size; ++i) {
		aCoef.push_back(yVec[i - 1]);
	}
	return aCoef;
}

vector<double> find_spline_d_coef(vector<double> x, vector<double> c) {
	vector<double> d; double val;
	int size = c.end() - c.begin();

	for (int i = 0; i < size-1; ++i) {
		val = (c[i + 1] - c[i]) / (3 * (x[i+1] - x[i]));
		d.push_back(val);
	}
	val = -c[size - 1] / (3 * (x[size - 1] - x[size - 2]));
	d.push_back(val);

	return d;
}

vector<double> find_spline_b_coef(vector<double> x, vector<double> y, vector<double> c) {
	vector<double> b; double val;
	int size = c.end() - c.begin();

	for (int i = 0; i < size-1; ++i) {
		val = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (x[i+1] - x[i]) * (c[i + 1] + 2 * c[i]) / 3;
		b.push_back(val);
	}
	val = (y[size - 1] - y[size - 2]) / (x[size - 1] - x[size - 2]) - (2 * c[size - 1] * (x[size - 1] - x[size - 2])) / 3;
	b.push_back(val);

	return b;
}

vector<double> find_spline_c_coef(vector<double> x, vector<double> y) {
	vector<double> dVec, cCoef;
	int size = y.end() - y.begin();
	vector<vector<double>> cMatrix = init_with_zeros(size - 1);
	int cSize = size - 1;

	for (int i = 1; i <= cSize; ++i) {
		if (i == 1) {
			cMatrix[0][0] = 1;
			dVec.push_back(0);
		}
		else if (i == cSize) {
			cMatrix[cSize - 1][cSize - 1] = 1;
			dVec.push_back(0);
		}
		else {
			cMatrix[i - 1][i - 2] = x[i - 1] - x[i - 2];
			cMatrix[i - 1][i-1] = 2 * (x[i-1]-x[i-2]+x[i]-x[i-1]);
			cMatrix[i - 1][i] = x[i] - x[i - 1];

			dVec.push_back( 3 * ((y[i]-y[i-1])/(x[i]-x[i-1]) - (y[i-1]-y[i-2])/(x[i-1]-x[i-2])) );
		}
	}
	cCoef = progonka(cMatrix, dVec);

	return cCoef;
}

vector<double> progonka(vector<vector<double>> cVec, vector<double> dVec) {
	vector<double> roots = copy_vec(dVec), a, b;
	int size = dVec.end() - dVec.begin();

	for (int i = 0; i < size-1; ++i) {
		if (i == 0) 
			a.push_back(-cVec[i][i+1]/cVec[i][i]);
		else 
			a.push_back(-cVec[i][i+1] / (cVec[i][i-1]*a[i-1]+cVec[i][i]));
	}

	for (int i = 0; i < size; ++i) {
		if (i == 0)
			b.push_back(dVec[i]/cVec[i][i]);
		else
			b.push_back((dVec[i]-cVec[i][i-1]*b[i-1])/(cVec[i][i-1]*a[i-1]+cVec[i][i]));
	}

	roots[size - 1] = b[size - 1];
	for (int i = size - 2; i > -1; --i) {
		roots[i] = a[i] * roots[i + 1] + b[i];
	}

	return roots;
}