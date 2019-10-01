#include "Header.h"
#include "foo.h"

Total gen_total()
{
	vector<int> inside;
	for (size_t i = 0; i < qpoints; ++i)
	{
		inside.push_back(i);
	}
	//  Eigen matrices default to column-major storage order
	//  so this is used here

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
			Vector3i tmp = get_vector(i);
			irrep.push_back(i);
			RED.push_back({});
			for (auto&& m : matrix)
			{
				int j = &m - &matrix[0];

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
	//cout << "runs here" << endl;

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
		indices.push_back(-1);
		vector<vector<int>> store_tmp_trips;

		for (auto&& i : inside)
		{

			if (find(indices.begin(), indices.end(), i) != indices.end()) {
				continue;
			}

			numbOfTriplets.back() += 1;
			Vector3i storeTMPvector;
			storeTMPvector << r, i, get_third(r, i);
			store_trips.push_back(storeTMPvector);
			store_tmp_trips.push_back({ r,i, get_third(r, i) });
			degen.push_back(1);
			indices.push_back(i);

			//apply all matrices
			Vector3i tmp = get_vector(i);
			for (auto&& m : matrix_index[&r - &irrep[0]])
			{
				Vector3i tmp1 = matrix[m] * tmp;
				move_inside_BZ(tmp1);
				int tmp_index = tmp1[0] + tmp1[1] * size_q + tmp1[2] * size_q * size_q;

				if (find(indices.begin(), indices.end(), tmp_index) == indices.end()) {
					indices.push_back(tmp_index);
					store_tmp_trips.push_back({ r,tmp_index, get_third(r, tmp_index) });
					degen.back() += 1;
				}

			}
			int store_size = store_tmp_trips.size();
			for (int i=0;i<store_size;i++)
			{
				vector<int> v = store_tmp_trips[i];
				vector<int> tmp_v = vec_flip(v);
				if (find(store_tmp_trips.begin(), store_tmp_trips.end(), tmp_v) == store_tmp_trips.end()) {
					store_tmp_trips.push_back(tmp_v);
					degen.back() += 1;
					indices.push_back(tmp_v[1]);
				}

			}
		}
	}
	/*
	//This code will plot the triplets and stuff done
	int sum_1 = 0;
	for (size_t i = 0; i < numbOfTriplets.size(); i++) {
		sum_1 += numbOfTriplets[i];
	}
	std::cout << "Sum: " << sum_1 << std::endl;
	int Sum_all_triplets = 0;
	for (size_t i = 0; i < (irrep.size()); i++) {
		Sum_all_triplets += numbOfTriplets[i];
		cout << "    " << endl;
		cout << irrep[i] << "   " << numbOfTriplets[i] << "   " << (matrix_index[i].size()) << endl;
		cout << "    " << endl;
		for (size_t j = 0; j < store_trips.size(); j++) {
			Vector3i a = store_trips[j];
			if (a(0) == irrep[i]) {
				cout << a[0] << "   " << a[1] << "   " << a[2] << "   " << degen[j] << endl;
			}
		}
	}
	std::cout << "All triplets: " << Sum_all_triplets << std::endl;
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

