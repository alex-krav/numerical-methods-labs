#include "io_func.h"
#include "matrix_func.h"
#include <cmath>

using namespace std;

vector<vector<double>> matrix_c(string, vector<vector<double>>);
vector<double> vector_d(string, vector<vector<double>>, vector<double>);
vector<double> vector_x_iteration(string, vector<vector<double>>, vector<double>, vector<vector<double>>, vector<double>, double);

bool x_less_m(vector<double>, vector<double>, double);
vector<double> multiply_matrix(vector<vector<double>>, vector<double>);
vector<double> add_matrix(vector<double>, vector<double>);

vector<double> vector_x_Seidel(string, vector<vector<double>>, vector<double>, double);
vector<double> Seidel_iteration(vector<vector<double>>, vector<double>, vector<double>);

int main2()
{
	string in = "input_lab2.txt", out = "output_lab2_iter.txt";
	vector< vector<double> > fullVec, aVec, cVec;
	vector<double> bVec, dVec, roots, mVec;
	double m = 0.000001; int variant = 2;

	read_vec_from_file(in, &fullVec, &aVec, &bVec);
	print_vec_to_cout("input matrix:", fullVec);
	save_vec_to_file(out, "input matrix:", fullVec);

	if (variant % 2 == 0) {
		cout << "\nSolving with simple iteration method" << endl;
		save_str_to_file(out, "\nSolving with simple iteration method");

		cVec = matrix_c(out, aVec);
		dVec = vector_d(out, aVec, bVec);

		roots = vector_x_iteration(out, aVec, bVec, cVec, dVec, m);
		print_vec_to_cout("roots:", roots);
		save_vec_to_file(out, "roots:", roots);

		mVec = mismatch_vector(bVec, aVec, roots);
		print_vec_to_cout("mismatch vector:", mVec);
		save_vec_to_file(out, "mismatch vector:", mVec);
	}
	else {
		cout << "\nSolving with Seidel method" << endl;
		save_str_to_file(out, "\nSolving with Seidel method");

		roots = vector_x_Seidel(out, aVec, bVec, m);
		print_vec_to_cout("roots:", roots);
		save_vec_to_file(out, "roots:", roots);

		mVec = mismatch_vector(bVec, aVec, roots);
		print_vec_to_cout("mismatch vector:", mVec);
		save_vec_to_file(out, "mismatch vector:", mVec);
	}

	return 0;
}

vector<vector<double>> matrix_c(string out, vector<vector<double>> aVec) {
	vector<vector<double>> aCopy = copy_vec(aVec);
	int aSize = aVec.end() - aVec.begin();

	for (int i = 0; i < aSize; ++i) {
		for (int j = 0; j < aSize; ++j) {
			if (i == j)
				aCopy[i][j] = 0;
			else
				aCopy[i][j] = -1 * aVec[i][j] / aVec[i][i];
		}
	}
	print_vec_to_cout("matrix C:", aCopy);
	save_vec_to_file(out, "matrix C:", aCopy);

	return aCopy;
}

vector<double> vector_d(string out, vector<vector<double>> aVec, vector<double> bVec) {
	vector<double> bCopy = copy_vec(bVec);
	int bSize = bVec.end() - bVec.begin();

	for (int i = 0; i < bSize; ++i) {
		bCopy[i] = bVec[i] / aVec[i][i];
	}
	print_vec_to_cout("vector d:", bCopy);
	save_vec_to_file(out, "vector d:", bCopy);

	return bCopy;
}

vector<double> vector_x_iteration(string out, vector<vector<double>> aVec, vector<double> bVec,
						vector<vector<double>> cVec, vector<double> dVec, double m) {
	vector<double> prev_roots = copy_vec(dVec), next_roots, mVec;
	int dSize = dVec.end() - dVec.begin();
	int iterNum = 0;

	while (true) {
		++iterNum;	
		char iter_mes[100];
		sprintf(iter_mes, "mismatch vector for iteration %d:", iterNum);
		mVec = mismatch_vector(bVec, aVec, prev_roots);
		print_vec_to_cout(iter_mes, mVec);
		save_vec_to_file(out, iter_mes, mVec);

     	// x(k+1) = c*x(k) + d
		vector<double> next_roots = add_matrix(multiply_matrix(cVec, prev_roots), dVec);
		if (x_less_m(next_roots, prev_roots, m)) {
			return next_roots;
		}
		else 
			prev_roots = copy_vec(next_roots);
	}
}

vector<double> vector_x_Seidel(string out, vector<vector<double>> aVec, vector<double> bVec, double m) {
	vector<double> prev_roots = vector_d(out, aVec, bVec), next_roots, mVec;
	int iterNum = 0;

	while (true) {
		char iter_mes[100];
		sprintf(iter_mes, "mismatch vector for iteration %d:", ++iterNum);
		mVec = mismatch_vector(bVec, aVec, prev_roots);
		print_vec_to_cout(iter_mes, mVec);
		save_vec_to_file(out, iter_mes, mVec);

		vector<double> next_roots = Seidel_iteration(aVec, bVec, prev_roots);
		if (x_less_m(next_roots, prev_roots, m)) {
			return next_roots;
		}
		else
			prev_roots = copy_vec(next_roots);
	}
}

bool x_less_m(vector<double> next_roots, vector<double> prev_roots, double m) {
	int size = next_roots.end() - next_roots.begin();
	double maxDiff = -1, diff;

	for (int i = 0; i < size; ++i) {
		diff = abs(next_roots[i] - prev_roots[i]);
		if (diff > maxDiff)
			maxDiff = diff;
	}
	return (maxDiff < m) ? true : false;
}

vector<double> multiply_matrix(vector<vector<double>> cVec, vector<double> xVec) {
	vector<double> cxVec = copy_vec(xVec);
	int xSize = xVec.end() - xVec.begin();
	double sum;

	for (int i = 0; i < xSize; ++i) {
		sum = 0;
		for (int j = 0; j < xSize; ++j) {
			sum += cVec[i][j] * xVec[j];
		}
		cxVec[i] = sum;
	}

	return cxVec;
}

vector<double> add_matrix(vector<double> cxVec, vector<double> dVec) {
	vector<double> xVec = copy_vec(cxVec);
	int xSize = xVec.end() - xVec.begin();

	for (int i = 0; i < xSize; ++i) {
		xVec[i] = cxVec[i] + dVec[i];
	}

	return xVec;
}

vector<double> Seidel_iteration(vector<vector<double>> aVec, vector<double> bVec, vector<double> prev_roots) {
	vector<double> roots = copy_vec(prev_roots);
	int xSize = prev_roots.end() - prev_roots.begin();
	double sum1, sum2;

	for (int i = 0; i < xSize; ++i) {
		sum1 = sum2 = 0;
		for (int j = 0; j < i; ++j) {
			sum1 += ((aVec[i][j] / aVec[i][i]) * roots[j]);
		}
		for (int j = i+1; j < xSize; ++j) {
			sum2 += ((aVec[i][j] / aVec[i][i]) * prev_roots[j]);
		}
		roots[i] = -1 * sum1 - sum2 + bVec[i]/aVec[i][i];
	}

	return roots;
}