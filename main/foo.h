#pragma once
#include "Header.h"

int subtractIndex(int, int);
int addIndex(int, int);
double dirac(double);

//Functions 
vector<double> f_ph(vector<double>&, vector<double>&, vector<double>&, MatrixPH, IRREP, PHtwo);
vector<double> f_mg_alpha(vector<double>&, vector<double>&, vector<double>&, IRREP, MatrixMG);
vector<double> f_mg_beta(vector<double>&, vector<double>&,vector<double>&,IRREP, MatrixMG);

//RKfour
void RKfour(vector<double>&, vector<double>&, vector<double>&,
	MatrixPH, IRREP, PHtwo, MatrixMG, MatrixMG, double h);

void readFiles(vector<double>&, string);

MatrixE createSmallerM_Fone(vector<double>, vector<int>, vector<int>, vector<vector<int>>,int);
MatrixE createSmallerM_Ftwo(vector<double>, vector<int>, vector<int>, vector<vector<int>>,int);
MatrixE createSmallerM_Fthree(vector<double>, vector<int>, vector<int>, vector<vector<int>>,int);

Vector3i get_vector(int, int);
vector<int> vec_flip(vector<int>);
void move_inside_BZ(Vector3i&);
int get_third(int, int, int);

Total gen_total();

int mEnergyOne(int, int, int, vector<double>);

int mEnergyTwo(int, int, int, vector<double>);

int mEnergyThree(int, int, int, vector<double>);

bool IsZero(int, vector<vector<int>>,int);

string test_of_matching(MatrixE);

string test_of_matching_TWO(MatrixE);

void find_the_largest(MatrixE);


