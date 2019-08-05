void RKfour(vector<double>& phonon, vector<double>& mg_alpha,
	vector<double>& mg_beta, vector<int> matrixAdd, vector<double> A, vector<double> B,
	double h) {

	//create vectors to keep track of input into next step
	vector<double> dyn1, dyn2, dyn3;

	//create matrices for values

	vector<double> temp1(size_A), temp2(size_A), temp3(size_A), temp4(size_A);

	editor(phonon, mg_alpha, mg_beta, matrixAdd, A, B, temp1, temp2, temp3, temp4);

	vector<double> k1ph = f_ph(phonon, mg_alpha, mg_beta, matrixAdd, temp1, temp2);
	vector<double> k1mg_a = f_mg_alpha(phonon, mg_alpha, mg_beta, matrixAdd, temp1, temp3);
	vector<double> k1mg_b = f_mg_beta(phonon, mg_alpha, mg_beta, matrixAdd, temp2, temp4);

	for (int i = 0; i < size_ph; i++) {
		dyn1.push_back(phonon[i] + (h / 2.0) * k1ph[i]);
	}
	for (int j = 0; j < qpoints; j++) {
		dyn2.push_back(mg_alpha[j] + (h / 2.0) * k1mg_a[j]);
		dyn3.push_back(mg_beta[j] + (h / 2.0) * k1mg_b[j]);
	}

	editor(dyn1, dyn2, dyn3, matrixAdd, A, B, temp1, temp2, temp3, temp4);

	vector<double> k2ph = f_ph(dyn1, dyn2, dyn3, matrixAdd, temp1, temp2);
	vector<double> k2mg_a = f_mg_alpha(dyn1, dyn2, dyn3, matrixAdd, temp1, temp3);
	vector<double> k2mg_b = f_mg_beta(dyn1, dyn2, dyn3, matrixAdd, temp2, temp4);

	for (int i = 0; i < size_ph; i++) {
		dyn1[i] = phonon[i] + h / 2.0 * k2ph[i];
	}
	for (int j = 0; j < qpoints; j++) {
		dyn2[j] = mg_alpha[j] + h / 2.0 * k2mg_a[j];
		dyn3[j] = mg_beta[j] + h / 2.0 * k2mg_b[j];
	}
	editor(dyn1, dyn2, dyn3, matrixAdd, A, B, temp1, temp2, temp3, temp4);

	vector<double> k3ph = f_ph(dyn1, dyn2, dyn3, matrixAdd, temp1, temp2);
	vector<double> k3mg_a = f_mg_alpha(dyn1, dyn2, dyn3, matrixAdd, temp1, temp3);
	vector<double> k3mg_b = f_mg_beta(dyn1, dyn2, dyn3, matrixAdd, temp2, temp4);

	for (int i = 0; i < size_ph; i++) {
		dyn1[i] = phonon[i] + h * k3ph[i];
	}
	for (int j = 0; j < qpoints; j++) {
		dyn2[j] = mg_alpha[j] + h * k3mg_a[j];
		dyn3[j] = mg_beta[j] + h * k3mg_b[j];
	}

	editor(dyn1, dyn2, dyn3, matrixAdd, A, B, temp1, temp2, temp3, temp4);

	vector<double> k4ph = f_ph(dyn1, dyn2, dyn3, matrixAdd, temp1, temp2);
	vector<double> k4mg_a = f_mg_alpha(dyn1, dyn2, dyn3, matrixAdd, temp1, temp3);
	vector<double> k4mg_b = f_mg_beta(dyn1, dyn2, dyn3, matrixAdd, temp2, temp4);

	for (int i = 0; i < size_ph; i++) {
		phonon[i] += h / 6.0 * (k1ph[i] + 2 * k2ph[i] + 2 * k3ph[i] + k4ph[i]);
	}
	for (int j = 0; j < qpoints; j++) {
		mg_alpha[j] += ((h / 6.0) * (k1mg_a[j] + 2 * k2mg_a[j] + 2 * k3mg_a[j] + k4mg_a[j]));
		mg_beta[j] += h / 6.0 * (k1mg_b[j] + 2 * k2mg_b[j] + 2 * k3mg_b[j] + k4mg_b[j]);
	}
}