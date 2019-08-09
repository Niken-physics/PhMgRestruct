#include "Header.h"
#include "foo.h"
#include <omp.h>
#define _USE_MATH_DEFINES
#include <fstream> 
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <iomanip>
#include<algorithm>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;



vector<double> f_ph(vector<double>& phonon, vector<double>& mg_alpha, vector<double>& mg_beta, MatrixPH MPH, IRREP irrep, PHtwo phTWO) {
	vector<double> RHS(size_ph);
	//std::cout << "HELLO FROM f_ph, I AM WORKING HARD TODAY :)" << std::endl;

	MatrixE tmp = MPH.A;
#pragma omp parallel for
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}

	tmp = MPH.B;
#pragma omp parallel for
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}

	tmp = MPH.C;
#pragma omp parallel for
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}
	tmp = MPH.D;


#pragma omp parallel for
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}
#pragma omp parallel for
	for (size_t count = 0; count <phTWO.One.theta.size(); count++)
	{
		RHS[phTWO.One.index[count][0]] += phTWO.One.theta[count] * phonon[phTWO.One.index[count][0]]; //nonsense but correct later
	}
#pragma omp parallel for
	for (size_t count = 0; count < phTWO.Two.theta.size(); count++)
	{
		RHS[phTWO.Two.index[count][0]] += phTWO.Two.theta[count] * phonon[phTWO.Two.index[count][0]]; //nonsense but correct later
	}

	for (size_t b = 0; b < branches; b++)
	{
		phonon[b * kpoints] = 0;
	}

#pragma omp parallel for
	for (size_t count = 0; count < irrep.irrep.size(); count++)
	{
		size_t i = irrep.irrep[count];
		for (size_t b = 0; b < branches; b++)
		{
			for (auto&& ind : irrep.RED[count]) {
				RHS[ind + b * irrep.C] = RHS[i + b * irrep.C];
			}
		}
	}
	return RHS;
}


//fcn used to update in RK-4 for magnons alpha
vector<double> f_mg_alpha(vector<double>& phonon, vector<double>& mg_alpha, IRREP irrep, MatrixMG MGA){
vector<double> RHS(qpoints);
//std::cout << "HELLO FROM f_mg_alpha, I AM WORKING HARD TODAY :)" << std::endl;
MatrixE tmp = MGA.A;
#pragma omp parallel for
for (size_t count = 0; count < tmp.m.size(); count++)
{
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	RHS[j] -= tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
}

tmp = MGA.B;
#pragma omp parallel for
for (size_t count = 0; count < tmp.m.size(); count++)
{
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	RHS[j] -= tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
}

RHS[0] = 0;

#pragma omp parallel for
for (size_t count = 0; count < irrep.irrep.size(); count++)
{
	size_t i = irrep.irrep[count];
	for (auto&& ind : irrep.RED[count])
	{
		RHS[ind] = RHS[i];
	}
}
	return RHS;
}

//fcn used to update in RK-4 for magnons beta

vector<double> f_mg_beta(vector<double>& phonon, vector<double>& mg_alpha, IRREP irrep, MatrixMG MGB) {
	vector<double> RHS(qpoints);
	//std::cout << "HELLO FROM f_mg_alpha, I AM WORKING HARD TODAY :)" << std::endl;
	MatrixE tmp = MGB.A;
#pragma omp parallel for
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[j] -= tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
	}

	tmp = MGB.B;
#pragma omp parallel for
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[j] -= tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
	}

	RHS[0] = 0;

#pragma omp parallel for
	for (size_t count = 0; count < irrep.irrep.size(); count++)
	{
		size_t i = irrep.irrep[count];
		for (auto&& ind : irrep.RED[count])
		{
			RHS[ind] = RHS[i];
		}
	}
	return RHS;
}
