#ifndef IO_FUNCTIONS_H_INCLUDED
#define IO_FUNCTIONS_H_INCLUDED

#include <string>
#include <vector>

using namespace std;

void print_vec_to_cout(string, vector<vector<double>>);
void print_vec_to_cout(string, vector<double>);

void read_vec_from_file(string, vector<vector<double>>*, vector<vector<double>>*, vector<double>*);
void read_vec_from_file(string, string, vector<vector<double>>*);
void read_vec_from_file(string, string, vector<double>*);
void save_vec_to_file(string, string, vector<double>);
void save_vec_to_file(string, string, vector<vector<double>>);
void save_str_to_file(string, string);
void print_file(string);

#endif
