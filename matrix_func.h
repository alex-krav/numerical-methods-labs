#ifndef MATRIX_FUNCTIONS_H_INCLUDED
#define MATRIX_FUNCTIONS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

bool is_symetric(vector<vector<double>>);
bool vec_equal(vector<vector<double>>, vector<vector<double>>);
vector<vector<double>> transposed_matrix(vector<vector<double>>);
vector<vector<double>> copy_vec(vector<vector<double>>);
vector<double> copy_vec(vector<double>);
vector<double> mismatch_vector(vector<double>, vector<vector<double>>, vector<double>);
vector<double> calc_values(vector<double>, double (*func)(double));
vector<vector<double>> init_with_zeros(int);

#endif
