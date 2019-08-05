//spar längsta raden
// stoppa när vi kommer till 0:a.


#define _USE_MATH_DEFINES
#include <fstream> 
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <iomanip>
#include<algorithm>

using namespace std;
//declare variables

const double k_B = 1.38064852 * pow(10, -23); // 
const double bohrM = 9.274009 * pow(10, -24); //  J/T
const double g = 2.0;
const double H_ext = 0;
const double M = g * bohrM * H_ext;

const int size_q = 4;
const int qpoints = size_q * size_q * size_q;
const int size_k = 4;
const int kpoints = size_k * size_k * size_k;
const int branches = 12;
const int size_ph = kpoints * branches;
const int size_Theta = size_ph * size_ph;
const int size_A = size_ph * qpoints;
const int size_add = kpoints * qpoints;

const double smear = pow(10, -22);

//declare functions
int add(int, int, int);

double dirac(double);

//f_ph(phonon,mg_alpha,mg_beta,matrixAdd, A,B,C,theta1,theta2,theta3,theta4,theta5,theta6)
vector<double> f_ph_One(vector<double>&, vector<double>&,
	vector<double>&, vector<int>,
	vector<double>, vector<double>);

vector<double> f_ph_Two(vector<double>&, vector<double>&,
	vector<double>&, vector<int>,
	vector<double>, vector<double>);

//mg_alpha(phonon,mg_alpha,mg_beta,matrixAdd, alpha1,alpha2)
vector<double> f_mg_alpha(vector<double>&, vector<double>&,
	vector<double>&, vector<int>, vector<double>, vector<double>);

//mg_beta(phonon_mg_alpha,mg_beta,matrixAdd,beta1,beta2)
vector<double> f_mg_beta(vector<double>&, vector<double>&,
	vector<double>&, vector<int>, vector<double>, vector<double>);

//RKfour
void RKfour(vector<double>&, vector<double>&,
	vector<double>&, vector<int>,
	vector<double>, vector<double>, double);

void readFiles(vector<double>&, string);

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
	myfile.open("AvTemp.txt");
	ofstream myfileOne;
	myfileOne.open("TotalEnergy.txt");
	ofstream myfileTwo;
	myfileTwo.open("tempMaxMin.txt");
	ofstream myfileThree;
	myfileThree.open("Conservation.txt");

	int time_max = 11;
	double h = 1.0;

	for (int t = 0; t < time_max; t++) {
		RKfour(phonon, mg_alpha, mg_beta, matrixAdd, A, B, h);
	}
	myfile.close();
	myfileOne.close();
	myfileTwo.close();
	myfileThree.close();
	std::cout << "Done" << std::endl;
	return 0;
}
