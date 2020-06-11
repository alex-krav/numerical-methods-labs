#include "io_func.h"
#include "matrix_func.h"

using namespace std;

vector<vector<double>> frobenius_normal_form(string, vector<vector<double>>);
vector < vector<double>> identity_matrix(int);
vector<vector<double>> find_m_vec(string, string, vector<vector<double>>, int, int);
vector<vector<double>> find_mI_vec(string, string, vector<vector<double>>, int, int);
vector<vector<double>> multiply_vec(vector<vector<double>>, vector<vector<double>>);

int main3()
{
	string in = "input_lab3.txt", out = "output_lab3.txt";
	vector< vector<double> > aVec, pVec;

	read_vec_from_file(in, out, &aVec);
	
	pVec = frobenius_normal_form(out, aVec);

	return 0;
}

vector<vector<double>> frobenius_normal_form(string out, vector<vector<double>> aVec) {

	vector<vector<double>> pVec = copy_vec(aVec);
	vector<vector<double>> matrix, matrix1;
	int m = aVec.end() - aVec.begin();

	for (int i = 1; i <= m-1; ++i) {
		char iter_mes[100];
		sprintf(iter_mes, "M(%d) matrix:", i);
		matrix = find_m_vec(out, iter_mes, pVec, m, i);
		sprintf(iter_mes, "M-1(%d) matrix:", i);
		matrix1 = find_mI_vec(out, iter_mes, pVec, m, i);
		pVec = multiply_vec(multiply_vec(matrix1, pVec), matrix);
	}
	print_vec_to_cout("frobenius normal form:", pVec);
	save_vec_to_file(out, "frobenius normal form:", pVec);

	return pVec;
}

vector < vector<double>> identity_matrix(int size) {
	vector < vector<double>> matrix;

	for (int i = 0; i < size; ++i) {
		vector<double> row;
		for (int j = 0; j < size; ++j) {
			if (i == j)
				row.push_back(1);
			else
				row.push_back(0);
		}
		matrix.push_back(row);
	}

	return matrix;
}

vector<vector<double>> find_m_vec(string out, string mes, vector<vector<double>> aVec, int m, int i) {

	vector<vector<double>> mVec = identity_matrix(m);
	
	for (int j = 1; j <= m; ++j) {
		if (j == m-i) 
			mVec[m-i-1][j-1] = 1 / aVec[m-i][m-i-1];
		else 
			mVec[m-i-1][j-1] = -aVec[m-i][j-1] / aVec[m-i][m-i-1];
	}
	print_vec_to_cout(mes, mVec);
	save_vec_to_file(out, mes, mVec);

	return mVec;
}

vector<vector<double>> find_mI_vec(string out, string mes, vector<vector<double>> aVec, int m, int i) {
	vector<vector<double>> mIVec = identity_matrix(m);

	for (int j = 1; j <= m; ++j) {
		mIVec[m-i-1][j-1] = aVec[m-i][j-1];
	}
	print_vec_to_cout(mes, mIVec);
	save_vec_to_file(out, mes, mIVec);

	return mIVec;
}

vector<vector<double>> multiply_vec(vector<vector<double>> a, vector<vector<double>> b) {
	vector<vector<double>> m = copy_vec(a);
	int size = a.end() - a.begin();
	double val;

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			val = 0;
			for (int k = 0; k < size; ++k) {
				val += a[i][k] * b[k][j];
			}
			m[i][j] = val;
		}
	}

	return m;
}