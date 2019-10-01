#pragma once
#define _USE_MATH_DEFINES
#include <fstream> 
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <iomanip>
#include<algorithm>
#include <Eigen/Dense> //Make sure to load this library

using namespace std;
using namespace Eigen;

/*
//Special declaration of reduction
//Important in order to avoid race conditions in RHS.cpp
#pragma omp declare reduction(vec_double_plus : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                     initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
*/

//Constant params used in simulation
// Parameters that can be changed: size_q,compare,smearing

constexpr double k_B = 1.38064852 * 1e-23; 
constexpr double bohrM = 9.274009 * 1e-24;// J/T
constexpr double g = 2.0; 
constexpr double H_ext = 0; //T
constexpr double M = g * bohrM * H_ext;
constexpr double compare = 1e-8; // Threshold for multiplicative factors, introduces small errors
// but saves computational time

constexpr int size_q = 4; //Resolution of lattice, uniform sampling of first BZ
constexpr int qpoints = size_q * size_q * size_q;
constexpr int size_k = size_q;
constexpr int kpoints = size_k * size_k * size_k;
constexpr int branches = 12;
constexpr int size_ph = kpoints * branches;
constexpr int size_add = kpoints * qpoints;
constexpr int size_ALL = size_ph * qpoints;

constexpr double smear = 1e-23; //J

constexpr double h_bar = 1.054571817 * 1e-34; // J*s
constexpr double preFactor = (1 / h_bar) * (1 / smear) * 1e-15 * 1e-38; //1e-15 gives femtosecond units, 1e-38 is units from input



struct MatrixE {
	vector<double> m;
	vector<int> kb;
	vector<int> qqP;
};

struct MatrixPH {
	MatrixE F1nz;
	MatrixE A;
	MatrixE B;
	MatrixE C;
	MatrixE D;
};

struct MatrixMG {
	MatrixE F1anz;
	MatrixE A;
	MatrixE B;
};
struct IRREP {
	int C;
	vector<int> irrep;
	vector<vector<int>> RED;
};
struct nzTRIP {
	vector<double> theta;
	vector<Vector4i> index;
};
struct PHtwo {
	nzTRIP One;
	nzTRIP Two;
	nzTRIP Three;
};
struct Total {
	IRREP irrep;
	vector<Vector4i> triplets;
};