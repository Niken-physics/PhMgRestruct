#include "Header.h"
#include "foo.h"

vector<double> f_ph(vector<double>& phonon, vector<double>& mg_alpha, vector<double>& mg_beta, MatrixPH MPH, IRREP irrep, PHtwo phTWO) {
	vector<double> RHS(size_ph);
	//std::cout << "HELLO FROM f_ph, I AM WORKING HARD TODAY :)" << std::endl;
	MatrixE tmp = MPH.F1nz;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (mg_alpha[j]*mg_beta[jPrim]-phonon[i]*mg_alpha[j]-phonon[i]*mg_beta[jPrim]-phonon[i]);
	}

	tmp = MPH.A;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[jPrim]);
	}

	tmp = MPH.B;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_alpha[jPrim] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[jPrim] - phonon[i] * mg_alpha[j]);
	}

	tmp = MPH.C;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_beta[j] + mg_beta[jPrim] * mg_beta[j] + mg_beta[j] - phonon[i] * mg_beta[jPrim]);
	}
	tmp = MPH.D;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[i] += tmp.m[count] * (phonon[i] * mg_beta[jPrim] + mg_beta[jPrim] * mg_beta[j] + mg_beta[jPrim] - phonon[i] * mg_beta[j]);
	}
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	/*for (int count = 0; count <phTWO.One.theta.size(); count++)
	{
		int k1 = phTWO.One.index[count][0];
		int k2 = phTWO.One.index[count][1];
		int k3 = phTWO.One.index[count][2];
		RHS[phTWO.One.index[count][0]] += phTWO.One.theta[count] * 
			(phonon[k2]*phonon[k3]+phonon[k3]+phonon[k1]*phonon[k3]-phonon[k1]*phonon[k2]); 
	}
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (int count = 0; count < phTWO.Two.theta.size(); count++)
	{
		int k1 = phTWO.Two.index[count][0];
		int k2 = phTWO.Two.index[count][1];
		int k3 = phTWO.Two.index[count][2];
		RHS[phTWO.Two.index[count][0]] += phTWO.Two.theta[count] *
			(phonon[k3] * phonon[k2] + phonon[k2] + phonon[k1] * phonon[k2] - phonon[k1] * phonon[k3]);
	}
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (int count = 0; count < phTWO.Three.theta.size(); count++)
	{
		int k1 = phTWO.Three.index[count][0];
		int k2 = phTWO.Three.index[count][1];
		int k3 = phTWO.Three.index[count][2];
		RHS[phTWO.Two.index[count][0]] += phTWO.Three.theta[count] *
			(phonon[k3] * phonon[k2] - phonon[k1] - phonon[k1] * phonon[k2] - phonon[k1] * phonon[k3]);
	}
	*/

	for (size_t b = 0; b < branches; b++)
	{
		phonon[b * kpoints] = 0;
	}

//#pragma omp parallel for //reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < irrep.irrep.size(); count++)
	{
		int i = irrep.irrep[count];
		for (int b = 0; b < branches; b++)
		{
			for (auto&& ind : irrep.RED[count]) {
				RHS[ind + b * kpoints] = RHS[i + b * kpoints];
			}
		}
	}

	return RHS;
}


//fcn used to update in RK-4 for magnons alpha
vector<double> f_mg_alpha(vector<double>& phonon, vector<double>& mg_alpha, vector<double>& mg_beta, IRREP irrep, MatrixMG MGA){
vector<double> RHS(qpoints);
MatrixE tmp = MGA.F1anz;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
for (size_t count = 0; count < tmp.m.size(); count++)
{
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	RHS[j] -= tmp.m[count] * (mg_alpha[j] * mg_beta[jPrim] - phonon[i] * mg_alpha[j] - phonon[i] * mg_beta[jPrim] - phonon[i]);
}
    tmp = MGA.A;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
for (size_t count = 0; count < tmp.m.size(); count++)
{
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	/*if (j == 0 || j == 10 || j == 34 || j == 40 || jPrim == 0 || jPrim == 10 || jPrim == 34 || jPrim == 40) {
		std::cout << "You are having an inapropiate interaction 1" << std::endl;
		std::cout << i << "  " << j << " " << jPrim << std::endl;
	}*/
	RHS[j] -= tmp.m[count] * (phonon[i] * mg_alpha[j] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[j] - phonon[i] * mg_alpha[jPrim]);
}

tmp = MGA.B;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
for (size_t count = 0; count < tmp.m.size(); count++)
{
	int i = tmp.kb[count];
	int j = tmp.qqP[count] % qpoints;
	int jPrim = tmp.qqP[count] / qpoints;
	/*if (j == 0 || j == 10 || j == 34 || j == 40 || jPrim == 0 || jPrim == 10 || jPrim == 34 || jPrim == 40) {
		std::cout << "You are having an inapropiate interaction 2" << std::endl;
		std::cout << i << "  " << j << " " << jPrim << std::endl;
	}*/
	/*if (j == 23) {
		std::cout << i << "  " << j << " " << jPrim <<"  " << tmp.m[count] << std::endl;
	}*/
	//std::cout << i << " " << j << " " << jPrim << std::endl;
	RHS[j] -= tmp.m[count] * (phonon[i] * mg_alpha[jPrim] + mg_alpha[jPrim] * mg_alpha[j] + mg_alpha[jPrim] - phonon[i] * mg_alpha[j]);
}

//#pragma omp parallel for //reduction(vec_double_plus:RHS)
for (size_t count = 0; count < irrep.irrep.size(); count++)
{
	int i = irrep.irrep[count];
	for (auto&& ind : irrep.RED[count])
	{
		RHS[ind] = RHS[i];
	}
}
	return RHS;
}

//fcn used to update in RK-4 for magnons beta

vector<double> f_mg_beta(vector<double>& phonon, vector<double>& mg_alpha, vector<double>& mg_beta, IRREP irrep, MatrixMG MGB) {
	vector<double> RHS(qpoints);
	MatrixE tmp = MGB.F1anz;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[j] -= tmp.m[count] * (mg_beta[j] * mg_alpha[jPrim] - phonon[i] * mg_beta[j] - phonon[i] * mg_alpha[jPrim] - phonon[i]);
	}

	tmp = MGB.A;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[j] -= tmp.m[count] * (phonon[i] * mg_beta[j] + mg_beta[jPrim] * mg_beta[j] + mg_beta[j] - phonon[i] * mg_beta[jPrim]);
	}

	tmp = MGB.B;
//#pragma omp parallel for reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < tmp.m.size(); count++)
	{
		int i = tmp.kb[count];
		int j = tmp.qqP[count] % qpoints;
		int jPrim = tmp.qqP[count] / qpoints;
		RHS[j] -= tmp.m[count] * (phonon[i] * mg_beta[jPrim] + mg_beta[jPrim] * mg_beta[j] + mg_beta[jPrim] - phonon[i] * mg_beta[j]);
	}

//#pragma omp parallel for //reduction(vec_double_plus:RHS)
	for (size_t count = 0; count < irrep.irrep.size(); count++)
	{
		int i = irrep.irrep[count];
		for (auto&& ind : irrep.RED[count])
		{
			RHS[ind] = RHS[i];
		}
	}
	return RHS;
}
