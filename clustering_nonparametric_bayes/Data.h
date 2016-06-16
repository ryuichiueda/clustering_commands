#ifndef __DATA_H___
#define __DATA_H___

#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <vector>
using namespace std;

class Data
{
public:
	vector<double> original_data;
	Eigen::VectorXd normalized_data;
	int cluster_id;
	bool target;

	Data()
	{
		cluster_id = 0;
		target = false;
	}

	void print(void)
	{
		for(auto e : original_data)
			cout << e << ' ';

		cout << cluster_id << endl;
	}
};

#endif
