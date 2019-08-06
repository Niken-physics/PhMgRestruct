//spar längsta raden
// stoppa när vi kommer till 0:a.
#define _USE_MATH_DEFINES
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <iomanip>
#include<algorithm>
#include <fstream> 
#include <cmath>
#include "Header.h"
#include "foo.h"

using namespace std;
//declare variables

int main(int argc, char const* argv[]) {
	//calculate the add matrix
	vector<int> matrixAdd;
	std::cout << "runs" << std::endl;

	for (int i = 0; i < size_add; i++) {
		matrixAdd.push_back(add(i % kpoints, i / kpoints, size_q));
	}
	vector<double> phonon, A(size_A), B(size_A), mg_alpha, mg_beta;

	//readfiles w_ph,w_mgBeta, w_mgAlpha, F1, F2, theta

	//w_ph
	vector<double> w_ph;
	readFiles(w_ph, "phonon.txt");

	//w_mg_alpha
	vector<double> w_mg_alpha;
	readFiles(w_mg_alpha, "magnon.txt");

	//w_mg_beta
	vector<double> w_mg_beta;
	readFiles(w_mg_beta, "magnon.txt");

	//F1,F2
	/*vector<double> F1;
	readFiles(F1, "F1four.txt");*/

	vector<double> F2;
	readFiles(F2, "F2four.txt");

	std::cout << w_mg_alpha.size() << std::endl;
	std::cout << phonon.size() << std::endl;


	// Get correct energies for w_mg_alpha and w_mg_beta
	for (int j = 0; j < qpoints; j++) {
		w_mg_alpha[j] += M;
		w_mg_beta[j] -= M;
		/*std::cout << w_mg_alpha[j] << std::endl;
		std::cout << w_mg_beta[j] << std::endl;*/
	}
	std::cout << "runs1" << std::endl;

	//calculate A,B
	for (int i = 0; i < size_ph; i++) {
		for (int j = 0; j < qpoints; j++) {
			int temp = j * size_ph + i;
			int temp_i = i % kpoints;

			double E3 = w_ph[i] - w_mg_alpha[j] + w_mg_alpha[matrixAdd[j + temp_i * kpoints]];
			double E4 = w_ph[i] + w_mg_alpha[j] - w_mg_alpha[matrixAdd[j + temp_i * kpoints]];
			double E5 = w_ph[i] - w_mg_beta[j] + w_mg_beta[matrixAdd[j + temp_i * kpoints]];
			double E6 = w_ph[i] + w_mg_beta[j] - w_mg_beta[matrixAdd[j + temp_i * kpoints]];
			A[temp] = 4.0 * M_PI * 1000.0 / 1.054 * F2[temp] * (dirac(E3));//This factor changes if you change params above
			B[temp] = 4.0 * M_PI * 1000.0 / 1.054 * F2[temp] * (dirac(E5));//This factor changes if you change params above
		}
	}

	// Calculate initial distributions;

	double T = 300.0;

	for (int i = 0; i < size_ph; i++) {
		if (i % 64 == 0) {
			phonon.push_back(0);
		}
		else {
			phonon.push_back(1 / (exp(w_ph[i] / k_B / T) - 1));
		}
	}
	T = 301.0;

	//Excite some phonon modes

	/*for(int i=0;i<branches;i++){
		phonon[i*64+1]+=10;
		phonon[i*64+5]+=10;
		phonon[i*64+17]+=10;
	}*/

	for (int j = 0; j < qpoints; j++) {
		mg_alpha.push_back(1 / (exp(w_mg_alpha[j] / k_B / T) - 1));
		mg_beta.push_back(1 / (exp(w_mg_beta[j] / k_B / T) - 1));
		if (j == 0) {
			mg_beta[j] = 0;
		}
		if (mg_alpha[j] >= 1) {
			mg_alpha[j] = 0;
		}
	}

	std::cout << "Runs this far" << std::endl;

	ofstream myfile;
	//myfile.open("AvTemp.txt");
	ofstream myfileOne;
	myfileOne.open("TotalEnergy.txt");
	//ofstream myfileTwo;
	//myfileTwo.open("tempMaxMin.txt");
	ofstream myfileThree;
	myfileThree.open("Conservation.txt");

	int time_max = 11;
	double h = 1.0;

	for (int t = 0; t < time_max; t++) {
		double E_ph = 0;
		double E_mg_alpha = 0;
		double E_mg_beta = 0;
		double phTot = 0;
		double mg_alphaTot = 0;
		double mg_betaTot = 0;
		for (auto&& i:irrep.IRREP)
		{
			for(size_t b=0; b<branches; b++){
				E_ph += phonon[i+kpoints*b] * w_ph[i+kpoints*b];
				phTot += phonon[i+kpoints*b];
			}
		}
		for (auto&& j:irrep.IRREP) 
		{
				E_mg_alpha += mg_alpha[j] * w_mg_alpha[j];
				E_mg_beta += mg_beta[j] * w_mg_beta[j];
				mg_alphaTot += mg_alpha[j];
				mg_betaTot += mg_beta[j];
		}
		RKfour(phonon, mg_alpha, mg_beta, matrixAdd, A, B, h);
	}
	//myfile.close();
	myfileOne.close();
	//myfileTwo.close();
	myfileThree.close();
	cout << "Done" << endl;
	return 0;
}
