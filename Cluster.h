#ifndef __CLUSTER_H_
#define __CLUSTER_H_

#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <random>
#include "ProbDistributions.h"
#include "DataSet.h"
using namespace std;

ProbDistributions pd;

class Cluster
{
public:
	vector<double> mean; //XXX: This should be Eigen::Vector2d
	Eigen::MatrixXd cov;
	vector<Data *> data;

	Cluster()
	{
		mean.push_back(pd.uniformRand(0.0,1.0));
		mean.push_back(pd.uniformRand(0.0,1.0));
		cov = Eigen::MatrixXd(2,2);
		cov << 0.5, 0.0, 0.0, 0.5;
	}

	void clear(void)
	{
		data.clear();
	}

	void regData(Data *d)
	{
		data.push_back(d);
	}

	void calcParams(void)
	{
		if(data.size() == 0)
			return;

		mean[0] = 0.0;
		mean[1] = 0.0;
		for(int j=0;j<2;j++){
			for(auto d : data){
				mean[j] += d->normalized_data[j];
			}
			mean[j] /= data.size();
		}

#if 0
		for(int j=0;j<2;j++){
			double sum = 0.0;
			for(int i=0;i<data.size();i++){
				double d = data[i]->normalized_data[j] - mean[j];
				sum += d*d;
			}
			if(j==0){
				cov(0,0) = sqrt(sum/data.size());
				if(cov(0,0) < 0.01)
					cov(0,0) = 0.01;
			}else{
				cov(1,1) = sqrt(sum/data.size());
				if(cov(1,1) < 0.01)
					cov(1,1) = 0.01;
			}
		}
#endif
		cov(0,0) = 0.02;
		cov(1,1) = 0.02;
	}

	void print(void)
	{
		if(data.size() == 0)
			return;

		cerr << setprecision(2);
		cerr << mean[0] << ',' << mean[1] << " cov: " << cov(0,0) << ',' << cov(1,1)
		<< " num: " << data.size() << endl;
	}
};

class Clusters
{
public:
	vector<Cluster> c;
	void calcParams(void)
	{
		for(auto &e : c){
			e.calcParams();
			e.print();
		}
	}
};

#endif
