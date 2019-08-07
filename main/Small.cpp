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
#include<Eigen/Dense>
#include "foo.h"
#include "Header.h"

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

MatrixE createSmallerM(vector<double> scatter, vector<int> index, vector<int> operation,vector<vector<int>> RED)
{
	MatrixE tmp;
	for (auto&& i : scatter)
	{
		auto count = &i - &scatter[0];
		if (abs(i) > compare) {
			tmp.m.push_back(i);
			tmp.kb.push_back(index[count] % size_ph);
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
			tmp.qqP.push_back(q + qP * qpoints);
		}
	}
	return tmp;
}
