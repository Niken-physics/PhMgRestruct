#include "Header.h"
#include "foo.h"

using namespace std;
using namespace Eigen;

int subtractIndex(int a, int b) {
	Vector3i v = get_vector(a, size_q);
	Vector3i w = get_vector(b, size_q);
	Vector3i u = v - w;
	move_inside_BZ(u);
	return u[0]+u[1]*size_q + u[2]*size_q*size_q;
}
int addIndex(int a, int b) {
	Vector3i v = get_vector(a, size_q);
	Vector3i w = get_vector(b, size_q);
	Vector3i u = v - w;
	move_inside_BZ(u);
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
Vector3i get_vector(int index, int size)
{
	Vector3i v;
	for (size_t i = 0; i < 3; i++)
	{
		v(i) = (index % size);
		index /= size;
	}
	return v;

}

void move_inside_BZ(Vector3i& v) {
	for (size_t i = 0; i < 3; i++)
	{
		//if (find(a.begin(), a.end(), v[i]) == a.end() && v[i]!=a.back())
		v[i] = v[i] % size_q;
		if (v[i] < 0) {
			v[i] += size_q;
		}
	}
}

int get_third(int a, int b, int q) {
	Vector3i vec_a = get_vector(a, q);
	Vector3i vec_b = get_vector(b, q);
	Vector3i vec_c;
	vec_c << 4, 4, 4;
	vec_c = vec_c - vec_a - vec_b;
	move_inside_BZ(vec_c);
	return (vec_c[0] + vec_c[1] * q + vec_c[2] * q * q);
}

vector<int> vec_flip(vector<int> v) {
	vector<int> tmp;
	tmp.push_back(v[0]);
	tmp.push_back(v[2]);
	tmp.push_back(v[1]);
	return tmp;
}

bool IsZero(int a, vector<vector<int>> irrep){
	if (find(irrep[0].begin(), irrep[0].end(), a) == irrep[0].end()) {
		return false;
	}
	return true;
}

MatrixE createSmallerM_Ftwo(vector<double> scatter, vector<int> index, vector<int> operation,vector<vector<int>> RED)
{
	MatrixE tmp;
	for (auto&& i : scatter)
	{
		auto count = &i - &scatter[0];
		if (abs(i) > compare) {
			int q = index[count] / size_ph;
			if (IsZero(q, RED)) {
				continue;
			}
			int k = index[count] % kpoints;
			if (IsZero(k, RED)) {
				continue;
			}
			int qP = operation[k + q * kpoints];
			if (IsZero(qP, RED)) {
				continue;
			}
			tmp.m.push_back(i);
			tmp.kb.push_back(index[count] % size_ph);
			tmp.qqP.push_back(q + qP * qpoints);
		}
	}
	return tmp;
}
MatrixE createSmallerM_Fone(vector<double> scatter, vector<int> index, vector<int> operation, vector<vector<int>> RED)
{
	MatrixE tmp;
	for (auto&& i : scatter)
	{
		auto count = &i - &scatter[0];
		if (abs(i) > compare) {
	
			int q = index[count] / size_ph;
			if (IsZero(q, RED)) {
				continue;
			}
			int k = index[count] % kpoints;
			if (IsZero(k, RED)) {
				continue;
			}
			int qP = operation[q + k * kpoints];
			if (IsZero(qP, RED)) {
				continue;
			}
			tmp.m.push_back(i);
			tmp.kb.push_back(index[count] % size_ph);
			tmp.qqP.push_back(q + qP * qpoints);
		}
	}
	return tmp;
}

MatrixE createSmallerM_Fthree(vector<double> scatter, vector<int> index, vector<int> operation, vector<vector<int>> RED)
{
	MatrixE tmp;
	for (auto&& i : scatter)
	{
		auto count = &i - &scatter[0];
		if (abs(i) > compare) {
			int q = index[count] / size_ph;
			if (IsZero(q, RED)) {
				continue;
			}
			int k = index[count] % kpoints;
			int b = index[count] / kpoints;
			if (IsZero(k, RED)) {
				continue;
			}
			Vector3i tmp_v = get_vector(k, size_q);
			tmp_v = -1 * tmp_v;
			move_inside_BZ(tmp_v);
			k = tmp_v[0] + tmp_v[1] * size_q + tmp_v[2] * size_q * size_q;
			int qP = operation[q + k * kpoints];
			if (IsZero(qP, RED)) {
				continue;
			}
			tmp.qqP.push_back(q + qP * qpoints);
			tmp.m.push_back(i);
			tmp.kb.push_back(k + kpoints * b);
		}
	}
	return tmp;
}


int mEnergyOne(int a, int b, int c, vector<double> w_ph) {
		int out = c;
		double temp = abs(w_ph[a] + w_ph[b] - w_ph[c]);
		for (int i = 1; i < branches; i++) {
			if (abs(w_ph[a] + w_ph[b] - w_ph[c + i * kpoints]) < temp) {
				int out = c + kpoints * i;
				double temp = abs(w_ph[a] + w_ph[b] - w_ph[c + i * kpoints]);
			}
		}
		return out;
	}

int mEnergyTwo(int a, int b, int c, vector<double> w_ph) {
	int out = c;
	double temp = abs(w_ph[a] - w_ph[b] + w_ph[c]);
	for (int i = 0; i < branches; i++) {
		if (abs(w_ph[a] - w_ph[b] + w_ph[c + i * kpoints]) < temp) {
			out = c + kpoints * i;
			temp = abs(w_ph[a] - w_ph[b] + w_ph[c + i * kpoints]);
		}
	}
	return out;
}

int mEnergyThree(int a, int b, int c, vector<double> w_ph) {
	int out = c;
	double temp = abs(w_ph[a] - w_ph[b] - w_ph[c]);
	for (int i = 0; i < branches; i++) {
		if (abs(w_ph[a] - w_ph[b] - w_ph[c + i * kpoints]) < temp) {
			out = c + kpoints * i;
			temp = abs(w_ph[a] - w_ph[b] - w_ph[c + i * kpoints]);
		}
	}
	return out;
}
