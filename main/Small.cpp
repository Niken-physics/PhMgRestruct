#define _USE_MATH_DEFINES
#include<string>
#include <fstream> 
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <iomanip>
#include<algorithm>

using namespace std;
extern const double smear;

int add(int a, int b, int size_q) {
	vector<int> v, w, u;

	for (int i = 0; i < 3; i++) {
		v.push_back(a % size_q);
		w.push_back(b % size_q);
		u.push_back(v[i] - w[i]);
		if (u[i] < 0) {
			u[i] += size_q;
		}
		b /= size_q;
		a /= size_q;
	}
	return u[0] + u[1] * size_q + u[2] * size_q * size_q;
}

double dirac(double E) {
	return exp(-(E * E) / (smear * smear)) / sqrt(M_PI);
}

void readFiles(vector<double> & w, string title) {
	double element;
	ifstream in;
	in.open(title);
	if (in.is_open()) {
		while (in >> element) {
			w.push_back(element);
		}
	}
	in.close();
}
