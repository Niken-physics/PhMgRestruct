#define _USE_MATH_DEFINES
#include <fstream> 
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>
#include <array>
#include <iomanip>
#include<algorithm>
#include <Eigen/Dense>
#include "foo.h"
#include "Header.h"

using namespace Eigen;
using namespace std;

MatrixE createSmallerM(vector<double> scatter,vector<int> index,vector<int> add)
{
	MatrixE tmp;
	for (auto&& i : scatter)
	{
		auto count = &i - &scatter[0];
		if (abs(i) > compare) {
			tmp.m.push_back(i);
			tmp.kb.push_back(index[count]%size_ph);
			int q = index[count] / size_ph;
			int k = index[count] % kpoints;
			int qP = add[k + q * kpoints];
			tmp.qqP.push_back(q+qP*qpoints);
		}
	}
	return tmp;
}