#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "matrix_func.h"

using namespace std;

bool is_symetric(vector<vector<double>> vec) {
	vector<vector<double>> tVec = transposed_matrix(vec);
	return vec_equal(vec, tVec);
}

bool vec_equal(vector<vector<double>> first, vector<vector<double>> second) {
	for (vector< vector<double> >::iterator fRow = first.begin(), sRow = second.begin();
		fRow != first.end() && sRow != second.end();
		++fRow, ++sRow) {
		for (vector<double>::iterator fCol = fRow->begin(), sCol = sRow->begin();
			fCol != fRow->end() && sCol != sRow->end();
			++fCol, ++sCol) {
			if (*fCol != *sCol)
				return false;
		}
	}

	return true;
}

vector<vector<double>> transposed_matrix(vector<vector<double>> vec) {
	vector<vector<double>> tVec;
	int size = vec.end() - vec.begin();

	for (int i = 0; i < size; ++i) {
		vector<double> tRow;
		for (int j = 0; j < size; ++j) {
			tRow.push_back(vec[j][i]);
		}
		tVec.push_back(tRow);
	}

	return tVec;
}

vector<vector<double>> copy_vec(vector<vector<double>> vec) {
	vector<vector<double>> copy;

	for (vector< vector<double> >::iterator row = vec.begin(); row != vec.end(); ++row) {
		vector<double> copyRow;
		for (vector<double>::iterator col = row->begin(); col != row->end(); ++col) {
			copyRow.push_back(*col);
		}
		copy.push_back(copyRow);
	}

	return copy;
}

vector<double> copy_vec(vector<double> vec) {
	vector<double> copy;

	for (vector<double>::iterator col = vec.begin(); col != vec.end(); ++col) {
		copy.push_back(*col);
	}

	return copy;
}

// r = b – Ax
vector<double> mismatch_vector(vector<double> bVec, vector<vector<double>> aVec, vector<double> roots) {
	vector<double> rVec;
	int size = bVec.end() - bVec.begin();
	for (int i = 0; i < size; ++i) {
		double aVal = 0;
		for (int j = 0; j < size; ++j) {
			aVal += (aVec[i][j] * roots[j]);
		}
		rVec.push_back(bVec[i] - aVal);
	}

	return rVec;
}

vector<double> calc_values(vector<double> xVec, double (*func)(double)) {
	vector<double> yVec;
	for (vector<double>::iterator elem = xVec.begin(); elem != xVec.end(); ++elem) {
		yVec.push_back(func(*elem));
	}
	return yVec;
}

vector<vector<double>> init_with_zeros(int n) {
	vector<vector<double>> zero_matrix;

	for (int i = 0; i < n; ++i) {
		vector<double> zero_row;
		for (int j = 0; j < n; ++j) {
			zero_row.push_back(0);
		}
		zero_matrix.push_back(zero_row);
	}

	return zero_matrix;
}