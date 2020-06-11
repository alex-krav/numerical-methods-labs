#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "io_func.h"

using namespace std;

void read_vec_from_file(string in, string out, vector<double>* fullVec) {
	ifstream input(in);
	if (input.is_open())
	{
		int n;
		input >> n;
		double elem;

		for (int j = 0; j < n; ++j) {
			input >> elem;
			fullVec->push_back(elem);
		}
		input.close();

		print_vec_to_cout("input vector:", *fullVec);
		save_vec_to_file(out, "input vector:", *fullVec);
	}
	else cout << "Unable to open input file: " << in << endl;
}

void read_vec_from_file(string in, string out, vector<vector<double>>* fullVec) {
	ifstream input(in);
	if (input.is_open())
	{
		int n, m;
		input >> n >> m;
		double elem;

		for (int i = 0; i < n; ++i) {
			vector<double> fullRow;
			for (int j = 0; j < m; ++j) {
				input >> elem;
				fullRow.push_back(elem);
			}
			fullVec->push_back(fullRow);
		}

		input.close();

		string mes = (n == 1) ? "input vector:" : "input matrix:";
		print_vec_to_cout(mes, *fullVec);
		save_vec_to_file(out, mes, *fullVec);
	}
	else cout << "Unable to open input file: " << in << endl;
}

void read_vec_from_file(string filename, vector<vector<double>>* fullVec, vector<vector<double>>* aVec, vector<double>* bVec) {
	ifstream input(filename);
	if (input.is_open())
	{
		int n, m;
		input >> n >> m;
		double elem;

		for (int i = 0; i < n; ++i) {
			vector<double> fullRow, aRow;
			for (int j = 0; j < m; ++j) {
				input >> elem;
				fullRow.push_back(elem);
				if (j < m - 1)
					aRow.push_back(elem);
				if (j == (m - 1))
					bVec->push_back(elem);
			}
			fullVec->push_back(fullRow);
			aVec->push_back(aRow);
		}

		input.close();
	}
	else cout << "Unable to open input file: " << filename << endl;
}

void print_vec_to_cout(string mes, vector<vector<double>> vec) {
	cout << endl;
	cout << mes << endl;
	for (vector< vector<double> >::iterator row = vec.begin(); row != vec.end(); ++row) {
		for (vector<double>::iterator col = row->begin(); col != row->end(); ++col) {
			printf("%12.6f", *col);
		}
		cout << endl;
	}
}

void print_vec_to_cout(string mes, vector<double> vec) {
	cout << endl;
	cout << mes << endl;
	for (vector<double>::iterator it = vec.begin(); it != vec.end(); ++it) {
		printf("%12.6f\n", *it);
	}
}

void save_str_to_file(string out, string mes) {
	ofstream output(out, ios_base::app);
	if (output.is_open())
	{
		output << mes << endl;
		output.close();
	}
	else cout << "Unable to open output file " << out << endl;;
}

void save_vec_to_file(string out, string mes, vector<double> vec) {
	ofstream output(out, ios_base::app);
	char buffer[100];
	if (output.is_open())
	{
		output << endl;
		output << mes << endl;
		for (vector<double>::iterator it = vec.begin(); it != vec.end(); ++it) {
			snprintf(buffer, sizeof(buffer), "%12.6f", *it);
			output << buffer << endl;
		}
		output.close();
	}
	else cout << "Unable to open output file " << out << endl;
}

void save_vec_to_file(string out, string mes, vector<vector<double>> vec) {
	ofstream  output(out, ios_base::app);
	char buffer[100];
	if (output.is_open())
	{
		output << endl;
		output << mes << endl;
		for (vector< vector<double> >::iterator row = vec.begin(); row != vec.end(); ++row) {
			for (vector<double>::iterator col = row->begin(); col != row->end(); ++col) {
				snprintf(buffer, sizeof(buffer), "%12.6f", *col);
				output << buffer;
			}
			output << endl;
		}
		output.close();
	}
	else cout << "Unable to open output file " << out << endl;
}

void print_file(string filename) {
	string line;
	ifstream input(filename);
	if (input.is_open())
	{
		while (input) {
			getline(input, line);
			cout << line << endl;
		}

		input.close();
	}
	else cout << "Unable to open input file " << filename << endl;
}