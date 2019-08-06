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
#include "Header.h"
#include "foo.h"
using namespace Eigen;
using namespace std;



vector<double> f_ph(vector<double>& phonon, vector<double>& mg_alpha, vector<double>& mg_beta, MatrixPH MPH, IRREP irrep, PHtwo phTWO) {
	vector<double> RHS(size_ph);

	MatrixE tmp = MPH.A;
#pragma omp parallel for
	for (auto&& element : tmp.m) 
	{
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}

	tmp = MPH.B;
#pragma omp parallel for
	for (auto&& element : tmp.m) 
	{
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
	}

	tmp = MPH.C;
#pragma omp parallel for
	for (auto&& element : tmp.m) 
	{
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
	}
	tmp = MPH.D;
#pragma omp parallel for
	for (auto&& element : tmp.m) 
	{
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j]);
	}

	for (auto&& t : phTWO.One.triplets)
	{
		auto count = &t - &phTWO.One.triplets[0];
		RHS[t[0] + irrep.C * t[3]] += phTWO.One.theta[count] * phonon[t[0]]; //nonsense but correct later
	}
#pragma omp parallel for
	for (auto&& t : phTWO.Two.triplets)
	{
		auto count = &t - &phTWO.Two.triplets[0];
		RHS[t[0] + irrep.C * t[3]] += phTWO.Two.theta[count] * phonon[t[0]]; //nonsense but correct later
	}

#pragma omp parallel for
	for (auto&& i : irrep.IRREP)
	{
		for (size_t b = 0; b < branches; b++)
		{
			for (auto&& ind : irrep.RED[&i - &irrep.IRREP[0]]) {
				RHS[ind + b * irrep.C] == RHS[i + b * irrep.C];
			}
		}
	}
	return RHS;
}


//fcn used to update in RK-4 for magnons alpha
vector<double> f_mg_alpha(vector<double>& phonon, vector<double>& mg_alpha, IRREP irrep, MatrixMG MGA){
vector<double> RHS(qpoints);
MatrixE tmp = MGA.A;
#pragma omp parallel for
for (auto&& element : tmp.m)
{
	auto count = &element - &tmp.m[0];
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
}
tmp = MGA.B;
#pragma omp parallel for
for (auto&& element : tmp.m)
{
	auto count = &element - &tmp.m[0];
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
}
#pragma omp parallel for
for (auto&& i : irrep.IRREP)
{
	for (size_t b = 0; b < branches; b++)
	{
		for (auto&& ind : irrep.RED[&i - &irrep.IRREP[0]]) {
			RHS[ind + b * irrep.C] == RHS[i + b * irrep.C];
		}
	}
}
	return RHS;
}



vector<double> f_mg_alpha(vector<double>& phonon, vector<double>& mg_alpha, IRREP irrep, MatrixMG MGB) {
	vector<double> RHS(qpoints);
	MatrixE tmp = MGB.A;
#pragma omp parallel for
	for (auto&& element : tmp.m)
	{
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}
	tmp = MGB.B;
#pragma omp parallel for
	for (auto&& element : tmp.m)
	{
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - mg_alpha[i] * mg_alpha[j]);
	}
#pragma omp parallel for
	for (auto&& i : irrep.IRREP)
	{
		for (size_t b = 0; b < branches; b++)
		{
			for (auto&& ind : irrep.RED[&i - &irrep.IRREP[0]]) {
				RHS[ind + b * irrep.C] == RHS[i + b * irrep.C];
			}
		}
	}
	return RHS;
}
//fcn used to update in RK-4 for magnons beta
