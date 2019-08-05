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
using namespace Eigen;
using namespace std;
using std::vector;
const double compare = pow(10, -6);
struct mIndex
{
	vector<double> m;
	vector<int> Index;
};
mIndex createSmallerM(mIndex arbit)
{
	mIndex tmp;
	for (auto&& i : arbit.m)
	{
		auto count = &i - &arbit.m[0];
		if (abs(i) > compare) {
			tmp.m.push_back(i);
			tmp.Index.push_back(arbit.Index[count]);
		}
	}
	return tmp;
}