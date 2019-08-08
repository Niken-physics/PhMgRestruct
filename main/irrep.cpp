#include "Header.h"
#include "foo.h"
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

using namespace std;
using namespace Eigen;


Total gen_total()
{
	vector<int> inside;
	for (size_t i = 0; i != qpoints; ++i)
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
	vector<vector<int>> RED; //store the indices of equivalent points
	for (auto&& i : inside)
	{
		if (find(index_vector.begin(), index_vector.end(), i) == index_vector.end()) {
			Vector3i tmp = get_vector(i, size_q);
			irrep.push_back(i);
			RED.push_back({});
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
					RED.back().push_back(tmp_index);
					
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
	/*for (size_t i = 0; i < (irrep.size()); i++) {
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
	vector<Vector4i> trip;
	for (auto&& d : degen) 
{
	auto i = &d - &degen[0];
	trip.push_back({ store_trips[i][0], store_trips[i][1],store_trips[i][2],d });
}
	int C = irrep.size();
	IRREP tmpIrrep{C,irrep,RED};
	Total total{ tmpIrrep,trip };

	return total;
}

