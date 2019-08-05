// creatingIrredTriplets.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Find irreducible triplets (q1,q2,q3) such that q1+q2+q3=1 (vector sum)

//1. get irred points, matrices and degeneracy (weight)
//2. apply matrices over all points record degeneracy 
// 3. "compare triplets"

//useful functions: add reciprocal vectors to get back to 1:st BZ
//read matrices
//apply matrices
//get vector from index

#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
//declare variables
using std::vector;
const int size_q = 4; //number of points sampled in each direction
const int tot_size = size_q * size_q * size_q; // 3D
//const int branches //number of branches

Vector3i get_vector(int, int);
vector<int> vec_flip(vector<int>);
void move_inside_BZ(Vector3i&);
int get_third(int, int, int);

int main()
{
	vector<int> inside;
	for (size_t i = 0; i != tot_size; ++i)
	{
		inside.push_back(i);
	}
	// read matrices (how to store them?) Store them in a vector

	Matrix3i A0, A1, A2, A3, A4, A5;
	A0 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	A1 << 0, 1, 0, 0, 0, 1, 1, 0, 0;
	A2 << 0, 0, 1, 1, 0, 0, 0, 1, 0;
	A3 << 0, 0, -1, 0, -1, 0, -1, 0, 0;
	A4 << -1, 0, 0, 0, 0, -1, 0, -1, 0;
	A5 << 0, -1, 0, -1, 0, 0, 0, 0, -1;
	vector<Matrix3i> matrix; //store matrices
	matrix.push_back(A0);
	matrix.push_back(A1);
	matrix.push_back(A2);
	matrix.push_back(A3);
	matrix.push_back(A4);
	matrix.push_back(A5);
	matrix.push_back(-A0);
	matrix.push_back(-A1);
	matrix.push_back(-A2);
	matrix.push_back(-A3);
	matrix.push_back(-A4);
	matrix.push_back(-A5);



	vector<int> index_vector; //store the indices of points we have seen
	vector<int> irrep; //store the indices of irreducible points
	vector<vector<int>> matrix_index; //store indices of matrices
	for (auto&& i : inside)
	{
		if (find(index_vector.begin(), index_vector.end(), i) == index_vector.end()) {
			Vector3i tmp = get_vector(i, size_q);
			irrep.push_back(i);
			for (auto&& m : matrix)
			{
				auto j = &m - &matrix[0];

				if (j == 0) {
					index_vector.push_back(i);
					matrix_index.push_back({ j });
					continue;
				}
				Vector3i tmp1 = m * tmp;
				move_inside_BZ(tmp1);
				int tmp_index = tmp1[0] + tmp1[1] * size_q + tmp1[2] * size_q * size_q;

				if (tmp_index == i) {

					matrix_index.back().push_back(j);
					continue;
				}

				if (find(index_vector.begin(), index_vector.end(), tmp_index) == index_vector.end()) {
					index_vector.push_back(tmp_index);
				}

			}
		}
	}
	cout << "runs here" << endl;

	for (size_t i = 0; i < irrep.size(); i++)
	{
		//cout << irrep[i] << endl;
	}

	vector<Vector3i> store_trips;
	vector<int> degen;
	vector<int> numbOfTriplets;

	for (auto&& r : irrep) {
		numbOfTriplets.push_back(0);
		vector<int> indices;

		for (auto&& i : inside)
		{
			vector<vector<int>> store_tmp_trips;

			if (find(indices.begin(), indices.end(), i) != indices.end()) {
				continue;
			}

			numbOfTriplets.back() += 1;
			Vector3i storeTMPvector;
			storeTMPvector << r, i, get_third(r, i, size_q);
			store_trips.push_back(storeTMPvector);
			store_tmp_trips.push_back({ r,i, get_third(r, i,size_q) });
			degen.push_back(1);
			indices.push_back(i);

			//apply all matrices
			Vector3i tmp = get_vector(i, size_q);

			for (auto&& m : matrix_index[&r - &irrep[0]])
			{
				Vector3i tmp1 = matrix[m] * tmp;
				move_inside_BZ(tmp1);
				int tmp_index = tmp1[0] + tmp1[1] * size_q + tmp1[2] * size_q * size_q;

				if (find(indices.begin(), indices.end(), tmp_index) == indices.end()) {
					indices.push_back(tmp_index);
					store_tmp_trips.push_back({ r,i, get_third(r, tmp_index,size_q) });
					degen.back() += 1;
				}

			}
			for (auto&& v : store_tmp_trips)
			{
				vector<int> tmp_v = vec_flip(v);

				if (find(store_tmp_trips.begin(), store_tmp_trips.end(), tmp_v) == store_tmp_trips.end()) {
					degen.back() += 1;
					indices.push_back(tmp_v[1]);
				}

			}

		}
	}
	for (size_t i = 0; i < (irrep.size()); i++) {
		cout << "    " << endl;
		cout << irrep[i] << "   " << numbOfTriplets[i] << "   " << (matrix_index[i].size()) << endl;
		cout << "    " << endl;
		/*for (size_t j = 0; j < store_trips.size(); j++) {
			Vector3i a = store_trips[j];
			if (a(0) == irrep[i]) {
				cout << a[0] << "   " << a[1] << "   " << a[2] << "   " << degen[j] << endl;
			}
		}
		*/
	}

	return 0;
}

