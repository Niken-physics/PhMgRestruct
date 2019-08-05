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
using std::vector;

extern const double k_b;
extern const double bohrM;
extern const double g;
extern const double H_ext;
extern const double M;

extern const int size_q;
extern const int qpoints;
extern const int size_k;
extern const int kpoints;
extern const int branches;
extern const int size_ph;
extern const int size_Theta;
extern const int size_A;
extern const int size_add;

extern const double smear;

//fcns used to update in RK-4 for phonons
//f_ph, f_mg_alpha f_mg_beta are all below
/*vector<double> f_ph_One(vector<double>& phonon, vector<double>& mg_alpha,
	vector<double>& mg_beta, vector<int> matrixA, vector<int> matrixS, vector<int> Aaphq,
	vector<int> Abphq, vector<int> Baphq, vector<int> Bbphq,
	vector<double> Aaph, vector<double> Baph, vector<double> Abph, vector<double> Bbph, vector<double> IRREP, int C, vector<vector<int>> RED) {
	vector<double> RHS(size_ph);

#pragma omp parallel for
	for (auto&& k : IRREP)
	{
		int a = &k - &IRREP[0];
		for (size_t b = 0; b < branches; b++)
		{
			int i = k + C * b;
			for (size_t tmp = a + 1; Aaph[tmp] != 0 || Baph[tmp] != 0 || Abph[tmp] != 0 || Bbph[tmp] != 0; tmp += C)
			{
				//getting the q-index ...
				int jAa = Aaphq[tmp];
				int jAb = Abphq[tmp];
				int jBa = Baphq[tmp];
				int jBb = Bbphq[tmp];
				int sSAa = matrixS[a + C * jAa];
				int sSAb = matrixS[a + C * jAb];
				int sABa = matrixA[a + C * jBa];
				int sABb = matrixA[a + C * jBb];

				RHS[i] += Aaph[tmp] * (phonon[i] * mg_alpha[jAa] + mg_alpha[jAa] * mg_alpha[sSAa] + mg_alpha[jAa] - phonon[i] * mg_alpha[sSAa]) +
					Abph[tmp] * (phonon[i] * mg_beta[jAb] + mg_beta[jAb] * mg_beta[sSAb] + mg_beta[jAb] - phonon[i] * mg_beta[sSAb]) +
					Baph[tmp] * (phonon[i] * mg_alpha[sABa] + mg_alpha[jBa] * mg_alpha[sABa] + mg_alpha[sABa] - phonon[i] * mg_alpha[jBa]) +
					Bbph[tmp] * (phonon[i] * mg_beta[sABb] + mg_beta[jBb] * mg_beta[sABb] + mg_beta[sABb] - phonon[i] * mg_beta[jBb]);


			}

			for (auto&& ind : RED[a])
			{
				RHS[ind+b*C] == RHS[i];
			}
		}
	}
	return RHS;
}*/
vector<double> f_ph_One(State state, matrixPH, IRREP  iRREP) {
	vector<double> RHS(size_ph);

	MatrixE tmp = MatrixPH.A;
#pragma omp parallel for
	for (auto&& element : tmp.m) {
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.q[count];
		int jPrim = tmp.q[count];
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j])
	}

	tmp = MatrixPH.B;
#pragma omp parallel for
	for (auto&& element : tmp.m) {
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.q[count];
		int jPrim = tmp.q[count];
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j])
	}

	tmp = MatrixPH.C;
#pragma omp parallel for
	for (auto&& element : tmp.m) {
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.q[count];
		int jPrim = tmp.q[count];
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j])
	}
	tmp = MatrixPH.D;
#pragma omp parallel for
	for (auto&& element : tmp.m) {
		auto count = &element - &tmp.m[0];
		int i = tmp.kb[count];
		int j = tmp.q[count];
		int jPrim = tmp.q[count];
		RHS[i] += element * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[j])
	}
#pragma omp parallel for
	for (auto&& i : IRREP)
	{
		for (size_t b = 0; b < branches; b++)
		{
			for (auto&& ind : RED[&i - &IRREP[0]]) {
				RHS[ind + b * C] == RHS[i + b * C];
			}
		}
	}
}

vector<double> f_ph_Two(vector<double> & phonon, IRREP irrep)


vector<vector<int>> RED, vector<Vector4i> TripOne, vector<Vector4i> TripTwo, vector<double> ThetaOne, vector<double> ThetaTwo, vector<int> degen, int C) {
vector<double> RHS(size_ph);
#pragma omp parallel  for
for (auto&& t : TripOne)
{
	auto count = &t - &TripOne[0];
	RHS[t[0] + C * t[3]] += ThetaOne[count] * phonon[t[0]]; //nonsense but correct later
}
#pragma omp parallel for
for (auto&& t : TripTwo)
{
	auto count = &t - &TripTwo[0];
	RHS[t[0] + C * t[3]] += ThetaTwo[count] * phonon[t[0]]; //nonsense but correct later
}
#pragma omp parallel for
for (auto&& i : IRREP)
{
	for (size_t b = 0; b < branches; b++)
	{
		for (auto&& ind : RED[&i - &IRREP[0]]) {
			RHS[ind + b * C] == RHS[i + b * C];
		}
	}
}

return RHS;
}


//fcn used to update in RK-4 for magnons alpha
vector<double> f_mg_alpha(vector<double> & phonon, vector<double> & mg_alpha,
	vector<double> & mg_beta, vector<int> matrixA, vector<int> matrixS, vector<int> Aaphq,
	vector<int> Abphq, vector<int> Baphq, vector<int> Bbphq,
	vector<double> Aaph, vector<double> Baph, vector<double> Abph, vector<double> Bbph, vector<double> IRREP, int C, vector<vector<int>> RED) {
	vector<double> RHS(size_ph);
#pragma omp parallel for
	for (auto&& k : IRREP)
	{
		int a = &k - &IRREP[0];
		for (size_t b = 0; b < branches; b++)
		{
			int i = k + C * b;
			for (size_t tmp = a + 1; Aaph[tmp] != 0 || Baph[tmp] != 0 || Abph[tmp] != 0 || Bbph[tmp] != 0; tmp += C)
			{
				//getting the q-index ... 
				int jAa = Aaphq[tmp];
				int jAb = Abphq[tmp];
				int jBa = Baphq[tmp];
				int jBb = Bbphq[tmp];
				int sSAa = matrixS[a + C * jAa];
				int sSAb = matrixS[a + C * jAb];
				int sABa = matrixA[a + C * jBa];
				int sABb = matrixA[a + C * jBb];

				RHS[i] += Aaph[tmp] * (phonon[i] * mg_alpha[jAa] + mg_alpha[jAa] * mg_alpha[sSAa] + mg_alpha[jAa] - phonon[i] * mg_alpha[sSAa]) +
					Abph[tmp] * (phonon[i] * mg_beta[jAb] + mg_beta[jAb] * mg_beta[sSAb] + mg_beta[jAb] - phonon[i] * mg_beta[sSAb]) +
					Baph[tmp] * (phonon[i] * mg_alpha[sABa] + mg_alpha[jBa] * mg_alpha[sABa] + mg_alpha[sABa] - phonon[i] * mg_alpha[jBa]) +
					Bbph[tmp] * (phonon[i] * mg_beta[sABb] + mg_beta[jBb] * mg_beta[sABb] + mg_beta[sABb] - phonon[i] * mg_beta[jBb]);


			}

			for (auto&& ind : RED[a])
			{
				RHS[ind + b * C] == RHS[i];
			}
		}
	}
	return RHS;
}
//fcn used to update in RK-4 for magnons beta
