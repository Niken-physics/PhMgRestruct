//spar längsta raden
// stoppa när vi kommer till 0:a.
#include "Header.h"
#include "foo.h"

int main(int argc, char const* argv[]) {

	std::cout << "The program runs" << std::endl;
	Total total = gen_total();
	IRREP irrep = total.irrep;
	vector<Vector4i> triplets = total.triplets;

	//calculate the add matrix
	vector<int> matrixAdd;
	vector<int> matrixSub;
	

	for (int i = 0; i < size_add; i++) {
		matrixSub.push_back(subtractIndex(i % kpoints, i / kpoints));
	}
	for (int i = 0; i < size_add; i++) {
		matrixAdd.push_back(addIndex(i % kpoints, i / kpoints));
	}
	vector<double> phonon(size_ph), mg_alpha(qpoints), mg_beta(qpoints);

	//readfiles w_ph,w_mg_beta, w_mg_alpha, F2[k,q+k], F2[k,q+k], theta[k,k']

	//w_ph
	vector<double> w_ph;
	readFiles(w_ph, "phonon.txt");

	//w_mg_alpha
	vector<double> w_mg_alpha;
	readFiles(w_mg_alpha, "magnon.txt");

	//w_mg_beta
	vector<double> w_mg_beta;
	readFiles(w_mg_beta, "magnon.txt");

	//F2[k,q] "F2one",F2[k,q+k] "F2two"
	vector<double> F1one;
	readFiles(F1one, "F1four.txt");

	vector<double> F2one;
	readFiles(F2one, "F2four.txt");

	vector<double> F2two;
	readFiles(F2two, "F2four.txt");

	vector<double> theta;
	//readFiles(theta, "Theta.txt");
	for (size_t i = 0; i < (branches * branches * triplets.size()); i++)
	{
		theta.push_back(0.0001);

	}


	//std::cout << w_mg_alpha.size() << std::endl;
	//std::cout << phonon.size() << std::endl;


	// Get correct energies for w_mg_alpha and w_mg_beta
	for (int j = 0; j < qpoints; j++) {
		w_mg_alpha[j] += M;
		w_mg_beta[j] -= M;
	}

	//calculate A,B,C,D,E,F
	vector<double> phF1, ph1, ph2, ph3, ph4, mgF1a, mgF1b, mga1, mga2, mgb1, mgb2;
	vector<int> phIndex, mgIndex;
	for (auto&& k:irrep.irrep) 
	{
		for(size_t b=0;b<branches;b++)
		{ 
			for (int j = 0; j < qpoints; j++) 
			{
				int i = kpoints * b + k;
				int temp = j * size_ph + i;
				double F1 = w_ph[i] - w_mg_alpha[j] - w_mg_beta[matrixSub[k + j * kpoints]];
				double E1 = w_ph[i] - w_mg_alpha[j] + w_mg_alpha[matrixSub[j + k * kpoints]];
				double E2 = w_ph[i] + w_mg_alpha[j] - w_mg_alpha[matrixAdd[j + k * kpoints]];
				double E3 = w_ph[i] - w_mg_beta[j] + w_mg_beta[matrixSub[j + k * kpoints]];
				double E4 = w_ph[i] + w_mg_beta[j] - w_mg_beta[matrixAdd[j + k * kpoints]];
				phF1.push_back(4 * M_PI * preFactor * F1one[temp] * dirac(E1));
				ph1.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E1));
				ph2.push_back(2.0 * M_PI * preFactor * F2two[temp] * dirac(E2));
				ph3.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E3));
				ph4.push_back (2.0 * M_PI * preFactor * F2two[temp] * dirac(E4));
				phIndex.push_back(i);
			}
		}
	}

	for (auto&& q : irrep.irrep)
	{
		for (size_t i = 0; i < size_ph; i++)
		{
			int k = i % kpoints;
			int temp = size_ph * q + i;
			double F1a = w_ph[i] - w_mg_alpha[q] - w_mg_beta[matrixSub[q + k * kpoints]];
			double F1b = w_ph[i] - w_mg_beta[q] - w_mg_alpha[matrixSub[q + k * kpoints]];
			double E1 = w_ph[i] - w_mg_alpha[q] + w_mg_alpha[matrixSub[q + k * kpoints]];
			double E2 = w_ph[i] + w_mg_alpha[q] - w_mg_alpha[matrixSub[q + k * kpoints]];
			double E3 = w_ph[i] - w_mg_beta[q] + w_mg_beta[matrixSub[q + k * kpoints]];
			double E4 = w_ph[i] + w_mg_beta[q] - w_mg_beta[matrixSub[q + k * kpoints]];

			mgF1a.push_back(4.0 * M_PI * preFactor * F1one[temp] * dirac(F1a));
			mgF1b.push_back(4.0 * M_PI * preFactor * F1one[temp] * dirac(F1b));
			mga1.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E1));
			mga2.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E2));
			mgb1.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E3));
			mgb2.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E4));
			mgIndex.push_back(i + q * size_ph);

		}
	}


	MatrixE F1nz,A, B, C, D, F1anz,F1bnz, Aa, Ba, Ab, Bb;
	F1nz = createSmallerM_Fone(phF1, phIndex, matrixSub, irrep.RED);
	A = createSmallerM_Ftwo(ph1, phIndex, matrixSub,irrep.RED);
	B = createSmallerM_Ftwo(ph2, phIndex, matrixAdd, irrep.RED);
	C = createSmallerM_Ftwo(ph3, phIndex, matrixSub, irrep.RED);
	D = createSmallerM_Ftwo(ph4, phIndex, matrixAdd, irrep.RED);
	F1anz = createSmallerM_Fone(mgF1a, mgIndex, matrixSub, irrep.RED);
	Aa = createSmallerM_Ftwo(mga1, mgIndex, matrixSub, irrep.RED);
	Ba = createSmallerM_Fthree(mga2, mgIndex, matrixSub, irrep.RED);
	F1bnz = createSmallerM_Fone(mgF1b, mgIndex, matrixSub, irrep.RED);
	Ab = createSmallerM_Ftwo(mgb1, mgIndex, matrixSub, irrep.RED);
	Bb = createSmallerM_Fthree(mgb2, mgIndex, matrixSub, irrep.RED);

	const MatrixPH MPH = { F1nz, A, B, C, D };
	const MatrixMG MGA = {F1anz, Aa,Ba };
	const MatrixMG MGB = {F1anz, Ab,Bb };

	int sizeTRIP = triplets.size();
	

	nzTRIP One;
	nzTRIP Two;
	nzTRIP Three;
	for(auto&& i:triplets)
	{

		auto count = &i - &triplets[0];
		for (size_t b = 0; b < branches; b++)
		{
			for (size_t b1 = 0; b1 < branches; b1++)
			{
				if (IsZero(i[0], irrep.RED)) {
					continue;
				}
				if (IsZero(i[1], irrep.RED)) {
					continue;
				}
				int kPPOne = mEnergyOne(i[0] + b*kpoints, i[1] + b1*kpoints, i[2], w_ph);
				int kPPTwo = mEnergyTwo(i[0] + b*kpoints, i[1] + b1*kpoints, i[2], w_ph);
				int kPPThree = mEnergyThree(i[0] + b*kpoints, i[1] + b1*kpoints, i[2], w_ph);

				
				if (IsZero(kPPOne, irrep.RED)==0) {
					double E1 = w_ph[i[0] + b * kpoints] + w_ph[i[1] + b1 * kpoints] - w_ph[kPPOne];
					double one = theta[count + sizeTRIP * b + sizeTRIP * branches * b1] * dirac(E1);
					if (one > compare) {
						One.theta.push_back(one);
						Vector4i v;
						v << i[0] + b*kpoints, i[1] + b1*kpoints, kPPOne, i[3];
						One.index.push_back(v);
					}
				}

				if (IsZero(kPPTwo, irrep.RED)==0) {
					double E2 = w_ph[i[0] + b * kpoints] - w_ph[i[1] + b1 * kpoints] + w_ph[kPPTwo];
					double two = theta[count + sizeTRIP * b + sizeTRIP * branches * b1] * dirac(E2);
					if (two > compare) {
						Two.theta.push_back(two);
						Vector4i v;
						v << i[0] + b*kpoints, i[1] + b1*kpoints, kPPTwo, i[3];
						Two.index.push_back(v);
					}
				}

				if (IsZero(kPPThree, irrep.RED) == 0) {
					double E3 = w_ph[i[0] + b * kpoints] - w_ph[i[1] + b1 * kpoints] - w_ph[kPPThree];
					double three = theta[count + sizeTRIP * b + sizeTRIP * branches * b1] * dirac(E3);
					if (three > compare) {
						Three.theta.push_back(three);
						Vector4i v;
						v << i[0] + b * kpoints, i[1] + b1 * kpoints, kPPThree, i[3];
						Three.index.push_back(v);
					}
				}

			}
		}
	}

	const PHtwo phTWO = { One,Two,Three };

	// Calculate initial distributions;

	double T = 300.0;

	for (auto&& i : irrep.irrep)
	{
		for (size_t b = 0; b < branches; b++)
		{
			phonon[i + b * kpoints] = (1 / (exp(w_ph[i + b * kpoints] / k_B / T) - 1));
		}
	}
	for (size_t b = 0; b < branches; b++) 
	{
		phonon[b * kpoints] = 0;
	}

	for (auto&& i : irrep.irrep)
	{
		for (size_t b = 0; b < branches; b++)
		{
			for (auto&& ind : irrep.RED[&i - &irrep.irrep[0]])
			{
				phonon[ind + b * kpoints] = phonon[i + b * irrep.C];
			}
		}
	}
	T = 301.0;

	//Excite some phonon modes

	/*for(int i=0;i<branches;i++){
		phonon[i*64+1]+=10;
		phonon[i*64+5]+=10;
		phonon[i*64+17]+=10;
	}*/

	for (auto&& q:irrep.irrep) {
		mg_alpha[q]=(1 / (exp(w_mg_alpha[q] / k_B / T) - 1));
		mg_beta[q] = (1 / (exp(w_mg_beta[q] / k_B / T) - 1));
	}
	mg_alpha[0] = 0;
	mg_beta[0] = 0;

	for (auto&& i : irrep.irrep)
	{
		for (auto&& ind : irrep.RED[&i - &irrep.irrep[0]])
		{
			mg_alpha[ind] = mg_alpha[i];
			mg_beta[ind] = mg_beta[i];
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

	int time_max = 11; //number of timesteps
	constexpr double h = 1.0; //choose 1 for fs


	for (int t = 0; t < time_max; t++) 
	{
		double E_ph = 0;
		double E_mg_alpha = 0;
		double E_mg_beta = 0;
		double phTot = 0;
		double mg_alphaTot = 0;
		double mg_betaTot = 0;
		for (auto&& i:irrep.irrep)
		{
			for(size_t b=0; b<branches; b++)
			{
				E_ph += phonon[i+kpoints*b] * w_ph[i+kpoints*b];
				phTot += phonon[i+kpoints*b];
			}
		}
		for (auto&& j:irrep.irrep) 
		{
				E_mg_alpha += mg_alpha[j] * w_mg_alpha[j];
				E_mg_beta += mg_beta[j] * w_mg_beta[j];
				mg_alphaTot += mg_alpha[j];
				mg_betaTot += mg_beta[j];
		}
		if (myfileOne.is_open()) {
			myfileOne << std::setprecision(10) << E_ph << "    " << E_mg_alpha << "   " << E_mg_beta
				<< '\n';
		}
		if (myfileThree.is_open()) {
			myfileThree << std::setprecision(10) << phTot << "    " << mg_alphaTot << "   " << mg_betaTot
				<< '\n';
		}
		std::cout << t << std::endl;

		RKfour(phonon, mg_alpha, mg_beta, MPH, irrep, phTWO, MGA, MGB,h);
	}


	//myfile.close();
	myfileOne.close();
	//myfileTwo.close();
	myfileThree.close();
	std::cout << "Done" << std::endl;
	return 0;
}
