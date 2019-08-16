#pragma once
#include "Header.h"

int subtractIndex(int, int);
int addIndex(int, int);
double dirac(double);

//f_ph(phonon,mg_alpha,mg_beta,matrixAdd, A,B,C,theta1,theta2,theta3,theta4,theta5,theta6)
vector<double> f_ph(vector<double>&, vector<double>&, vector<double>&, MatrixPH, IRREP, PHtwo);

//mg_alpha(phonon,mg_alpha,mg_beta,matrixAdd, alpha1,alpha2)
vector<double> f_mg_alpha(vector<double>&, vector<double>&, IRREP, MatrixMG);

//mg_beta(phonon_mg_alpha,mg_beta,matrixAdd,beta1,beta2)
vector<double> f_mg_beta(vector<double>&, vector<double>&,
IRREP, MatrixMG);

//RKfour
void RKfour(vector<double>&, vector<double>&, vector<double>&,
	MatrixPH, IRREP, PHtwo, MatrixMG, MatrixMG, double h);

void readFiles(vector<double>&, string);

MatrixE createSmallerM_Fone(vector<double>, vector<int>, vector<int>, vector<vector<int>>);
MatrixE createSmallerM_Ftwo(vector<double>, vector<int>, vector<int>, vector<vector<int>>);
MatrixE createSmallerM_Fthree(vector<double>, vector<int>, vector<int>, vector<vector<int>>);

Vector3i get_vector(int, int);
vector<int> vec_flip(vector<int>);
void move_inside_BZ(Vector3i&);
int get_third(int, int, int);

Total gen_total();

int mEnergyOne(int, int, int, vector<double>);

int mEnergyTwo(int, int, int, vector<double>);

int mEnergyThree(int, int, int, vector<double>);

bool IsZero(int, vector<vector<int>>);


