#include "io_func.h"
#include "matrix_func.h"

using namespace std;

vector<vector<double>> gauss_direct_way(string, vector<vector<double>>, vector<vector<double>>);
vector<double> gauss_return_way(vector<vector<double>>);

vector<vector<double>> square_tVec(vector<vector<double>>);
vector<double> square_yVec(vector<double>, vector<vector<double>>);
vector<double> square_xVec(vector<double>, vector<vector<double>>);

int main1()
{
	string in = "input_lab1_gauss.txt", out = "output_lab1_gauss.txt";
	vector< vector<double> > fullVec, aVec, taVec, gauss_direct, squareTVec;
	vector<double> bVec, roots, mVec, squareYVec;

	read_vec_from_file(in, &fullVec, &aVec, &bVec);
	print_vec_to_cout("input matrix:", fullVec);
	save_vec_to_file(out, "input matrix:", fullVec);

	if (is_symetric(aVec)) {
		cout << "\nSolving with square roots method" << endl;
		save_str_to_file(out, "\nSolving with square roots method");

		squareTVec = square_tVec(aVec);
		print_vec_to_cout("T matrix:", squareTVec);
		save_vec_to_file(out, "T matrix:", squareTVec);

		squareYVec = square_yVec(bVec, squareTVec);
		print_vec_to_cout("Y matrix:", squareYVec);
		save_vec_to_file(out, "Y matrix:", squareYVec);

		roots = square_xVec(squareYVec, squareTVec);
		print_vec_to_cout("roots:", roots);
		save_vec_to_file(out, "roots:", roots);

		mVec = mismatch_vector(bVec, aVec, roots);
		print_vec_to_cout("mismatch vector:", mVec);
		save_vec_to_file(out, "mismatch vector:", mVec);
	}
	else {
		cout << "\nSolving with Gauss method" << endl;
		save_str_to_file(out, "\nSolving with Gauss method");
		gauss_direct = gauss_direct_way(out, fullVec, aVec);

		roots = gauss_return_way(gauss_direct);
		print_vec_to_cout("roots vector:", roots);
		save_vec_to_file(out, "roots vector:", roots);

		mVec = mismatch_vector(bVec, aVec, roots);
		print_vec_to_cout("mismatch vector:", mVec);
		save_vec_to_file(out, "mismatch vector:", mVec);
	}

	return 0;
}

vector<vector<double>> gauss_direct_way(string out, vector<vector<double>> fullVec, vector<vector<double>> aVec) {
	vector<vector<double>> fullCopy = copy_vec(fullVec);
	vector<vector<double>> aCopy = copy_vec(aVec);
	int fullSize = fullVec[0].end() - fullVec[0].begin();
	int aSize = aVec.end() - aVec.begin();
	int iter_num = 0;

	for (int i = 0; i < aSize; ++i) {
		for (int j = i+1; j < aSize; ++j) {
			double m = fullCopy[j][i];
			double n = fullCopy[i][i];

			for (int k = i; k < fullSize; ++k) {
				fullCopy[j][k] -= fullCopy[i][k] * m / n;
				if (k < aSize)
					aCopy[j][k] = fullCopy[j][k];
			}
			char iter_mes[100];
			sprintf(iter_mes, "iteration %d:", ++iter_num);
			print_vec_to_cout(iter_mes, aCopy);
			save_vec_to_file(out, iter_mes, aCopy);
		}
	}
	print_vec_to_cout("return stroke matrix:", fullCopy);
	save_vec_to_file(out, "return stroke matrix:", fullCopy);

	return fullCopy;
}

vector<double> gauss_return_way(vector<vector<double>> vec) {
	int size = vec.end() - vec.begin();
	vector<double> roots(size);

	for (int i = size-1; i >= 0; --i) {
		double sum = 0;
		for (int j = size-1; j > i; --j) {
			sum += vec[i][j] * roots[j];
		}
		roots[i] = (vec[i][size] - sum) / vec[i][i];
	}

	return roots;
}

vector<vector<double>> square_tVec(vector<vector<double>> vec) {
	vector<vector<double>> tVec = copy_vec(vec);
	int size = vec.end() - vec.begin();

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {	
			if (i == 0 && j == 0)
				tVec[0][0] = sqrt(vec[0][0]);
			else if (i == 0 && j > 0)
				tVec[0][j] = vec[0][j] / tVec[0][0];
			else if (i == j && (0 < i < size)) {
				double sum = 0;
				for (int k = 0; k <= i - 1; ++k)
					sum += (tVec[k][i] * tVec[k][i]);
				tVec[i][i] = sqrt(vec[i][i] - sum);
			}
			else if (i < j) {
				double sum = 0;
				for (int k = 0; k <= i - 1; ++k)
					sum += (tVec[k][i] * tVec[k][j]);
				tVec[i][j] = (vec[i][j] - sum) / tVec[i][i];
			}
			else if (i > j)
				tVec[i][j] = 0;
		}
	}

	return tVec;
}

vector<double> square_yVec(vector<double> bVec, vector<vector<double>> tVec) {
	vector<double> yVec;
	int size = bVec.end() - bVec.begin();

	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int k = 0; k <= i - 1; ++k)
			sum += tVec[k][i] * yVec[k];
		double val = (bVec[i] - sum) / tVec[i][i];
		yVec.push_back(val);
	}

	return yVec;
}

vector<double> square_xVec(vector<double> yVec, vector<vector<double>> tVec) {
	vector<double> roots = copy_vec(yVec);
	int size = yVec.end() - yVec.begin();

	for (int i = size-1; i >= 0; --i) {
		double sum = 0;
		for (int k = i + 1; k < size; ++k) //i=0,k=1,2
			sum += tVec[i][k] * roots[k];
		double val = (yVec[i] - sum) / tVec[i][i];
		roots[i] = val;
	}

	return roots;
}
