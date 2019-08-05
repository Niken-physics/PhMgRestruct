double E_ph = 0;
double E_mg_alpha = 0;
double E_mg_beta = 0;
double phTot = 0;
double mg_alphaTot = 0;
double mg_betaTot = 0;
for (int i = 0; i < size_ph; i++) {
	if (i % 64 != 0) {
		E_ph += phonon[i] * w_ph[i];
		phTot += phonon[i];
	}
}
for (int j = 0; j < qpoints; j++) {
	if (j == 0 || j == 10 || j == 34 || j == 40) {
	}
	else {
		E_mg_alpha += mg_alpha[j] * w_mg_alpha[j];
		E_mg_beta += mg_beta[j] * w_mg_beta[j];
		mg_alphaTot += mg_alpha[j];
		mg_betaTot += mg_beta[j];
	}
}
double TavgF = 0.0;
double TavgMA = 0.0;
double TavgMB = 0.0;
vector<double> phonon_T;
for (int i = 0; i < size_ph; i++) {
	if (i % 64 != 0) {
		TavgF += w_ph[i] / k_B / log((1.0 / phonon[i]) + 1.0);
		phonon_T.push_back(w_ph[i] / k_B / log((1.0 / phonon[i]) + 1.0));
	}
}
vector<double> magnon_T;
for (int j = 0; j < qpoints; j++) {
	if (j == 0 || j == 10 || j == 34 || j == 40) {
	}
	else {
		TavgMA += w_mg_alpha[j] / k_B / log((1.0 / mg_alpha[j]) + 1.0);
		TavgMB += w_mg_beta[j] / k_B / log((1.0 / mg_beta[j]) + 1.0);
		magnon_T.push_back(w_mg_beta[j] / k_B / log((1.0 / mg_beta[j]) + 1.0));
	}
}
// add a comment;
double min = *min_element(magnon_T.begin(), magnon_T.end());
//std::cout<<"Min value: "<<min<<std::endl;
double max = *max_element(magnon_T.begin(), magnon_T.end());
//std::cout<<"Max value: "<<max<<std::endl;
double min1 = *min_element(phonon_T.begin(), phonon_T.end());
//std::cout<<"Min value: "<<min1<<std::endl;
double max1 = *max_element(phonon_T.begin(), phonon_T.end());
//std::cout<<"Max value: "<<max1<<std::endl;


TavgF = TavgF / (size_ph - 12);
TavgMA = TavgMA / (kpoints - 4);
TavgMB = TavgMB / (kpoints - 4);

//print_output to txt_file: E_ph, E_mg_alpha, E_mg_beta
double E_tot = (E_ph + E_mg_alpha + E_mg_beta);

if (myfile.is_open()) {
	myfile << std::setprecision(10) << TavgF << "    " << TavgMA
		<< "    " << TavgMB
		<< '\n';
}
if (myfileOne.is_open()) {
	myfileOne << std::setprecision(10) << E_ph << "    " << E_mg_alpha << "   " << E_mg_beta
		<< "    " << E_tot
		<< '\n';
}
if (myfileTwo.is_open()) {
	myfileTwo << std::setprecision(10) << max << "    " << min << "   " << max1
		<< "    " << min1
		<< '\n';
}
if (myfileThree.is_open()) {
	myfileThree << std::setprecision(10) << phTot << "    " << mg_alphaTot << "   " << mg_betaTot
		<< '\n';
}
