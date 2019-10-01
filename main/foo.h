#pragma once
#include "Header.h"

//Functions in Small.cpp
int subtractIndex(int, int);
int addIndex(int, int);
double dirac(double); //input = E, output = 1/sqrt(2pi)*e^-(1/2(E^2/smearing^2)) 
void readFiles(vector<double>&, string);

//Functions that keep multiplicative factors in the interaction that
// are larger than some threshold called compare
MatrixE createSmallerM_Fone(vector<double>, vector<int>, vector<int>, vector<int>);
MatrixE createSmallerM_Ftwo(vector<double>, vector<int>, vector<int>, vector<int>);
MatrixE createSmallerM_Fthree(vector<double>, vector<int>, vector<int>, vector<int>);


Vector3i get_vector(int); //Generates a vector from an index
vector<int> vec_flip(vector<int>); // Flips two indices of a vector (a,b,c) = (a,c,b)
void move_inside_BZ(Vector3i&); // moves a vector inside the first BZ by adding the reciprocal vector
int get_third(int, int); // Given k1,k2 in the ph-ph it generates k3

//Functions that find the branch of k3 in ph-ph interaction that
// gives the largest dirac delta
int mEnergyOne(int, int, int, vector<double>);
int mEnergyTwo(int, int, int, vector<double>);
int mEnergyThree(int, int, int, vector<double>);

bool IsZero(int, vector<int>);

//Functions used for testing
string test_of_matching(MatrixE);
string test_of_matching_TWO(MatrixE);
void find_the_largest(MatrixE);


//Function in irrep.cpp
Total gen_total();

//Functions in RHS.cpp
vector<double> f_ph(vector<double>&, vector<double>&, vector<double>&, MatrixPH, IRREP, PHtwo);
vector<double> f_mg_alpha(vector<double>&, vector<double>&, vector<double>&, IRREP, MatrixMG);
vector<double> f_mg_beta(vector<double>&, vector<double>&,vector<double>&,IRREP, MatrixMG);

//Function in RKfour.cpp
void RKfour(vector<double>&, vector<double>&, vector<double>&,
	MatrixPH, IRREP, PHtwo, MatrixMG, MatrixMG, double h);

