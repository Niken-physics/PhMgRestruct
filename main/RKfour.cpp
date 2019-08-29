#include "foo.h"
#include "Header.h"

void RKfour(vector<double>& phonon, vector<double>& mg_alpha, vector<double>& mg_beta,
	        MatrixPH MPH, IRREP irrep, PHtwo phTWO, MatrixMG MGA, MatrixMG MGB,double h) {
	//std::cout << "HELLO FROM RK4, I AM WORKING HARD TODAY :)" << std::endl;
	//create vectors to keep track of input into next step
	vector<double> dyn1, dyn2, dyn3;
	vector<double> k1ph = f_ph(phonon, mg_alpha, mg_beta, MPH, irrep, phTWO);
	vector<double> k1mg_a = f_mg_alpha(phonon, mg_alpha, mg_beta, irrep, MGA);
	vector<double> k1mg_b = f_mg_beta(phonon,mg_alpha ,mg_beta, irrep, MGB);

	for (int i = 0; i < size_ph; i++) {
		dyn1.push_back(phonon[i] + (h / 2.0) * k1ph[i]);
	}
	for (int j = 0; j < qpoints; j++) {
		dyn2.push_back(mg_alpha[j] + (h / 2.0) * k1mg_a[j]);
		dyn3.push_back(mg_beta[j] + (h / 2.0) * k1mg_b[j]);
	}

	vector<double> k2ph = f_ph(dyn1, dyn2, dyn3, MPH, irrep, phTWO);
	vector<double> k2mg_a = f_mg_alpha(dyn1, dyn2,dyn3, irrep, MGA);
	vector<double> k2mg_b = f_mg_beta(dyn1, dyn2, dyn3, irrep, MGB);

	for (int i = 0; i < size_ph; i++) {
		dyn1[i] = phonon[i] + h / 2.0 * k2ph[i];
	}
	for (int j = 0; j < qpoints; j++) {
		dyn2[j] = mg_alpha[j] + h / 2.0 * k2mg_a[j];
		dyn3[j] = mg_beta[j] + h / 2.0 * k2mg_b[j];
	}

	vector<double> k3ph = f_ph(dyn1, dyn2, dyn3, MPH, irrep, phTWO);
	vector<double> k3mg_a = f_mg_alpha(dyn1, dyn2, dyn3, irrep, MGA);
	vector<double> k3mg_b = f_mg_beta(dyn1, dyn2, dyn3, irrep, MGB);

	for (int i = 0; i < size_ph; i++) {
		dyn1[i] = phonon[i] + h * k3ph[i];
	}
	for (int j = 0; j < qpoints; j++) {
		dyn2[j] = mg_alpha[j] + h * k3mg_a[j];
		dyn3[j] = mg_beta[j] + h * k3mg_b[j];
	}

	vector<double> k4ph = f_ph(dyn1, dyn2, dyn3, MPH, irrep, phTWO);
	vector<double> k4mg_a = f_mg_alpha(dyn1, dyn2, dyn3, irrep, MGA);
	vector<double> k4mg_b = f_mg_beta(dyn1, dyn2, dyn3, irrep, MGB);

	for (int i = 0; i < size_ph; i++) {
		phonon[i] += h / 6.0 * (k1ph[i] + 2 * k2ph[i] + 2 * k3ph[i] + k4ph[i]);
	}

	dyn2 = mg_alpha;
	dyn3 = mg_beta;
	for (int j = 0; j < qpoints; j++) {
		mg_alpha[j] += ((h / 6.0) * (k1mg_a[j] + 2 * k2mg_a[j] + 2 * k3mg_a[j] + k4mg_a[j]));
		mg_beta[j] += h / 6.0 * (k1mg_b[j] + 2 * k2mg_b[j] + 2 * k3mg_b[j] + k4mg_b[j]);
	}


	for (int i = 0; i < irrep.irrep.size(); i++) {
		int q = irrep.irrep[i];
		std::cout << "i: " << q << "   Delta mg: " << mg_alpha[q] - dyn2[q] << std::endl;
	}
}