#include "Header.h"
#include "foo.h"

int main(int argc, char const* argv[]) {
	std::cout << "The program runs" << std::endl;
	//#pragma omp parallel
	//{std::cout << "Number of threads used: " <<omp_get_num_threads() << std::endl;
	//}
	Total total = gen_total(); //gen_total can be found in the file irrep.cpp
	IRREP irrep = total.irrep; //for more detail on IRREP, see the header file
	
	// In triplets there are 3 indices and 1 degeneracy stored
	// triplets[0] = b1k1, triplets[1] = b2k2, triplets[2]=k3, triplets[3]=degeneracy of triplet
	vector<Vector4i> triplets = total.triplets; 

	vector<int> degen; // Vector to save the degeneracy of the irreducible points
	// Note that RED does only include points that are equivalent to the
	// irreducible point, but not the irreducible point itself, hence the +1
	for (int count = 0; count < irrep.C; count++)
	{
		degen.push_back(irrep.RED[count].size()+1);
	}

	//calculate the add/subtract matrix 
	vector<int> matrixAdd; // q+k = q' Accessed in element q+k*kpoints
	vector<int> matrixSub; // q-k =q' Accessed in element q+k*kpoints
	// k-q =q' is accessed in element k+q*kpoints
	                              
	for (int i = 0; i < size_add; i++) {
		matrixSub.push_back(subtractIndex(i % kpoints, i / kpoints));
	}
	for (int i = 0; i < size_add; i++) {
		matrixAdd.push_back(addIndex(i % kpoints, i / kpoints));
	}
	vector<double> phonon(size_ph), mg_alpha(qpoints), mg_beta(qpoints);

	//Reading w_ph,w_mg_beta, w_mg_alpha, ,(F1[k,q])^2, (F2[k,q])^2, (F2[k,q+k])^2, (theta[k,k'])^2
	// Units of w_ph, w_mg_ are J
	// Units of F1^2,F2^2 and theta^2 are J^2*10^(-38)

	vector<double> w_ph, w_mg_alpha, w_mg_beta, F1, F2one, F2two, F,theta;
	readFiles(w_ph, "phonon-freq_4x4x4"); //change name depending on size 
	readFiles(w_mg_alpha, "magnon-freq_4x4x4"); //change name depending on size
	w_mg_beta = w_mg_alpha;
	readFiles(theta, "One-column-gkq-4x4x4-IRRED.dat"); //change name depending on size

	//F1[k,q] "F1" F2[k,q] "F2one",F2[k,q+k] "F2two"
	readFiles(F, "f_functions_s4"); //change name depending on size

	for(size_t i=0;i<F.size();i++)
	{
		if (i % 3 == 0) {
			F1.push_back(F[i]);
		}
		else if (i % 3 == 1) {
			F2one.push_back(F[i]);
		}
		else {
			F2two.push_back(F[i]);
		}
	}

	/*
	// Uncomment this and 3 other pieces of code to perform the index
	// matching test, comment reading of F above
	for(size_t q=0;q<qpoints;q++)
	{
		for (size_t k = 0; k < kpoints;k++)
		{
			for (size_t b = 0; b < branches; b++)
			{
				F1one.push_back(q * size_ph + k + b * kpoints);
				F2one.push_back(q * size_ph + k + b * kpoints);
				F2two.push_back(q * size_ph + k + b * kpoints);
			}
		}
	}
	*/

	// Get correct energies for w_mg_alpha and w_mg_beta
	// under an external magnetic field B_ext
	for (int j = 0; j < qpoints; j++) {
		w_mg_alpha[j] += M;
		w_mg_beta[j] -= M;
	}

	//We do not allow interactions with the magnons or phonons of 0 energy
	//They are in the indices k=(0,0,0) for phonons
	// and for q=(0,0,0), q=(1/2,1/2,0) as well as the equivalent points for magnons
	// The following creates a vector of q=(1/2,1/2,0) and its equivalent points

	int index_zero_for_mg = -1;
	for (size_t i = 0; i < irrep.irrep.size(); i++)
	{
		index_zero_for_mg += 1;
		if (irrep.irrep[i] == size_q * size_q / 2 + 2) {
			break;
		}
	}
	vector<int> of_limits_mg;
	of_limits_mg.push_back(size_q * size_q / 2 + 2);
	for (size_t i = 0; i < irrep.RED[index_zero_for_mg].size(); i++) {
		of_limits_mg.push_back(irrep.RED[index_zero_for_mg][i]);
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
				/*
				It is assumed that F1one,F2one,F2two are all ordered as qkb, where b changes the fastest,
				w_ph is has order bk where k changest the fastest.
				*/
				int i = kpoints * b + k;
				int tmp_i = b + branches * k;
				int temp = j * size_ph + tmp_i;
				double F1energy = w_ph[i] - w_mg_alpha[j] - w_mg_beta[matrixSub[k + j * kpoints]];
				double E1 = w_ph[i] - w_mg_alpha[j] + w_mg_alpha[matrixSub[j + k * kpoints]];
				double E2 = w_ph[i] + w_mg_alpha[j] - w_mg_alpha[matrixAdd[j + k * kpoints]];
				double E3 = w_ph[i] - w_mg_beta[j] + w_mg_beta[matrixSub[j + k * kpoints]];
				double E4 = w_ph[i] + w_mg_beta[j] - w_mg_beta[matrixAdd[j + k * kpoints]];
				phF1.push_back(4.0 * M_PI * preFactor * F1[temp] * dirac(F1energy));
				ph1.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E1));
				ph2.push_back(2.0 * M_PI * preFactor * F2two[temp] * dirac(E2));
				ph3.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E3));
				ph4.push_back(2.0 * M_PI * preFactor * F2two[temp] * dirac(E4));


				/*
				//To perform the index matching test, uncomment this and comment the
				// code right above it
				phF1.push_back(F1one[temp]);
				ph1.push_back(F2one[temp]);
				ph2.push_back(F2two[temp]);
				ph3.push_back(F2one[temp]);
				ph4.push_back(F2two[temp]);
				*/
				phIndex.push_back(i+j*size_ph);
			}
		}
	}

	for (auto&& q : irrep.irrep)
	{
		for (int i = 0; i < size_ph; i++)
		{
			/*
				It is assumed that F1one,F2one,F2two are all ordered as qkb, where b changes the fastest,
				w_ph is has order bk where k changest the fastest.
			*/
			int k = i % kpoints;
			int b = i / kpoints;
			int tmp_i = b + branches*k;
			int temp = q * size_ph + tmp_i;
			double F1a = w_ph[i] - w_mg_alpha[q] - w_mg_beta[matrixSub[k + q * kpoints]];
			double F1b = w_ph[i] - w_mg_beta[q] - w_mg_alpha[matrixSub[k + q * kpoints]];
			double E1 = w_ph[i] - w_mg_alpha[q] + w_mg_alpha[matrixSub[q + k * kpoints]];
			double E2 = w_ph[i] + w_mg_alpha[q] - w_mg_alpha[matrixSub[q + k * kpoints]];
			double E3 = w_ph[i] - w_mg_beta[q] + w_mg_beta[matrixSub[q + k * kpoints]];
			double E4 = w_ph[i] + w_mg_beta[q] - w_mg_beta[matrixSub[q + k * kpoints]];

			mgF1a.push_back(4.0 * M_PI * preFactor * F1[temp] * dirac(F1a));
			mgF1b.push_back(4.0 * M_PI * preFactor * F1[temp] * dirac(F1b));
			mga1.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E1));
			mga2.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E2));
			mgb1.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E3));
			mgb2.push_back(2.0 * M_PI * preFactor * F2one[temp] * dirac(E4));
	
			/*
			mgF1a.push_back(F1one[temp]);
			mgF1b.push_back(F1one[temp]);
			mga1.push_back(F2one[temp]);
			mga2.push_back(F2one[temp]);
			mgb1.push_back(F2one[temp]);
			mgb2.push_back(F2one[temp]);
			*/
			mgIndex.push_back(q*size_ph+b*kpoints+k);

		}
	}

	MatrixE F1nz,A, B, C, D, F1anz,F1bnz, Aa, Ba, Ab, Bb;
	F1nz = createSmallerM_Fone(phF1, phIndex, matrixSub, irrep.RED,index_zero_for_mg);
	A = createSmallerM_Ftwo(ph1, phIndex, matrixSub,irrep.RED, index_zero_for_mg);
	B = createSmallerM_Ftwo(ph2, phIndex, matrixAdd, irrep.RED, index_zero_for_mg);
	C = createSmallerM_Ftwo(ph3, phIndex, matrixSub, irrep.RED, index_zero_for_mg);
	D = createSmallerM_Ftwo(ph4, phIndex, matrixAdd, irrep.RED, index_zero_for_mg);
	F1anz = createSmallerM_Fone(mgF1a, mgIndex, matrixSub, irrep.RED, index_zero_for_mg);
	Aa = createSmallerM_Ftwo(mga1, mgIndex, matrixSub, irrep.RED, index_zero_for_mg);
	Ba = createSmallerM_Fthree(mga2, mgIndex, matrixSub, irrep.RED, index_zero_for_mg);
	F1bnz = createSmallerM_Fone(mgF1b, mgIndex, matrixSub, irrep.RED, index_zero_for_mg);
	Ab = createSmallerM_Ftwo(mgb1, mgIndex, matrixSub, irrep.RED, index_zero_for_mg);
	Bb = createSmallerM_Fthree(mgb2, mgIndex, matrixSub, irrep.RED, index_zero_for_mg);
	/*
	//This is the final thing to uncomment in order to see the matching of indices
	// Please note that the Ba and Bb saves the index -k instead of k
	// so there you have to another testing function
	std::cout << test_of_matching(F1nz) << std::endl;
	std::cout << test_of_matching(A) << std::endl;
    std::cout << test_of_matching(B) << std::endl;
	std::cout << test_of_matching(C) << std::endl;
	std::cout << test_of_matching(D) << std::endl;	
	std::cout << test_of_matching(F1anz) << std::endl;
	std::cout << test_of_matching(Aa) << std::endl;
	std::cout << test_of_matching_TWO(Ba) << std::endl;
	std::cout << test_of_matching(F1bnz) << std::endl;
	std::cout << test_of_matching(Ba) << std::endl;
	std::cout << test_of_matching_TWO(Bb) << std::endl;
	*/


	/*
	// The following is a test to get a sense of the size of the largest multiplicative factor
	// in the rate equations. 
	find_the_largest(Aa);
	find_the_largest(Ba);
	find_the_largest(Ab);
	find_the_largest(Bb);
	*/
	
	const MatrixPH MPH = { F1nz, A, B, C, D };
	const MatrixMG MGA = {F1anz, Aa,Ba };
	const MatrixMG MGB = {F1bnz, Ab,Bb };

	int sizeTRIP = triplets.size();

	nzTRIP One;
	nzTRIP Two;
	nzTRIP Three;

	//The following is done to set up things to be read in the correct way
	int index_theta = -1;
	vector<vector<Vector4i>> store_triplets_with_certain_index;
	int tmp_index_storage = -1;
	for (auto&& i : triplets) {
		if (i[0] == tmp_index_storage) {
			store_triplets_with_certain_index.back().push_back(i);
		}
		else {
			store_triplets_with_certain_index.push_back({i});
		}
		tmp_index_storage = i[0];
	}

	for(int b1 = 0;b1<branches;b1++)
	{
		for (auto&& r:irrep.irrep)
		{
			int count = &r - &irrep.irrep[0];
			for (int b2 = 0; b2 < branches; b2++)
			{
				for (auto&& i : store_triplets_with_certain_index[count])
				{
					index_theta++;
					if (i[0] ==0) {
						continue;
					}
					if (i[1] == 0) {
						continue;
					}
					int kPPOne = mEnergyOne(i[0] + b1 * kpoints, i[1] + b2 * kpoints, i[2], w_ph);
					int kPPTwo = mEnergyTwo(i[0] + b1 * kpoints, i[1] + b2 * kpoints, i[2], w_ph);
					int kPPThree = mEnergyThree(i[0] + b1 * kpoints, i[1] + b2 * kpoints, i[2], w_ph);


					if (kPPOne == 0) {
						double E1 = w_ph[i[0] + b1 * kpoints] + w_ph[i[1] + b2 * kpoints] - w_ph[kPPOne];
						double one = theta[index_theta] * dirac(E1);
						if (one > compare) {
							One.theta.push_back(one);
							Vector4i v;
							v << i[0] + b1 * kpoints, i[1] + b2 * kpoints, kPPOne, i[3];
							One.index.push_back(v);
						}
					}

					if (kPPTwo == 0) {
						double E2 = w_ph[i[0] + b1 * kpoints] - w_ph[i[1] + b2 * kpoints] + w_ph[kPPTwo];
						double two = theta[index_theta] * dirac(E2);
						if (two > compare) {
							Two.theta.push_back(two);
							Vector4i v;
							v << i[0] + b1 * kpoints, i[1] + b2 * kpoints, kPPTwo, i[3];
							Two.index.push_back(v);
						}
					}

					if (kPPThree == 0) {
						double E3 = w_ph[i[0] + b1 * kpoints] - w_ph[i[1] + b2 * kpoints] - w_ph[kPPThree];
						double three = theta[index_theta] * dirac(E3);
						if (three > compare) {
							Three.theta.push_back(three);
							Vector4i v;
							v << i[0] + b1 * kpoints, i[1] + b2 * kpoints, kPPThree, i[3];
							Three.index.push_back(v);
						}
					}
					
				}

			}
		}
	}
	const PHtwo phTWO = { One,Two,Three };

	// Calculate initial distributions;

	double T = 300.0; //phonon equilibrium temperature

	for (auto&& i : irrep.irrep)
	{
		for (int b = 0; b < branches; b++)
		{
			phonon[i + b * kpoints] = (1 / (exp(w_ph[i + b * kpoints] / k_B / T) - 1));
		}
	}

	// Enforce the 0 modes to have 0 phonons
	for (int b = 0; b < branches; b++) 
	{
		phonon[b * kpoints] = 0;
	}

	for (auto&& i : irrep.irrep)
	{
		for (int b = 0; b < branches; b++)
		{
			for (auto&& ind : irrep.RED[&i - &irrep.irrep[0]])
			{
				phonon[ind + b * kpoints] = phonon[i + b * kpoints];
			}
		}
	}
	
	//Excite some phonon modes
	//Should probab

	/*for(int i=0;i<branches;i++){
		phonon[i*64+1]+=10;
		phonon[i*64+5]+=10;
		phonon[i*64+17]+=10;
	}*/
	T = 300.0; //magnon equilibrium temperature
	for (auto&& q:irrep.irrep) {
		mg_alpha[q]=(1 / (exp(w_mg_alpha[q] / k_B / T) - 1));
		mg_beta[q] = (1 / (exp(w_mg_beta[q] / k_B / T) - 1));
	}
	mg_alpha[0] = 0;
	mg_beta[0] = 0;
	mg_alpha[size_q * size_q / 2 + 2] = 0;
	mg_beta[size_q * size_q / 2 + 2] = 0;
	for (auto&& i : irrep.irrep)
	{
		for (auto&& ind : irrep.RED[&i - &irrep.irrep[0]])
		{
			mg_alpha[ind] = mg_alpha[i];
			mg_beta[ind] = mg_beta[i];
		}
	}

	ofstream myfileOne, myfileThree;
	myfileOne.open("TotalEnergy.txt");
	myfileThree.open("Conservation.txt");

	int time_max = 100; //number of timesteps
	constexpr double h = 0.1; //choose 1 for fs

	//double START, END;
	//START = omp_get_wtime();
	for (int t = 0; t < time_max; t++) 
	{
		double E_ph = 0;
		double E_mg_alpha = 0;
		double E_mg_beta = 0;
		double phTot = 0;
		double mg_alphaTot = 0;
		double mg_betaTot = 0;
		for (int count = 0; count < irrep.C; count++) 
		{
			int i = irrep.irrep[count];
			for(int b=0; b<branches; b++)
			{
				E_ph += phonon[i+kpoints*b] * w_ph[i+kpoints*b]*degen[count];
				phTot += phonon[i+kpoints*b]*degen[count];
			}
		}
		for (int count = 0; count < irrep.C; count++)
		{
			int j = irrep.irrep[count];
			E_mg_alpha += mg_alpha[j] * w_mg_alpha[j]*degen[count];
			E_mg_beta += mg_beta[j] * w_mg_beta[j] * degen[count];
			mg_alphaTot += mg_alpha[j] * degen[count];
			mg_betaTot += mg_beta[j] * degen[count];
		}
		if (myfileOne.is_open()) {
			myfileOne << std::setprecision(10) << E_ph << "    " << E_mg_alpha << "   " << E_mg_beta
				<< '\n';
		}
		if (myfileThree.is_open()) {
			myfileThree << std::setprecision(10) << phTot << "    " << mg_alphaTot << "   " << mg_betaTot
				<< '\n';
		}
		RKfour(phonon, mg_alpha, mg_beta, MPH, irrep, phTWO, MGA, MGB,h);
	}
	//END = omp_get_wtime();
	//std::cout << "Loop time (s):" << END-START <<  std::endl;

	myfileOne.close();
	myfileThree.close();
	std::cout << "Done" << std::endl;
	return 0;
}
