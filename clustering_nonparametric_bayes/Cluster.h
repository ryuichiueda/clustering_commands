#ifndef __CLUSTER_H_
#define __CLUSTER_H_

#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <algorithm>
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
	Eigen::VectorXd mean;
	Eigen::MatrixXd cov;
	vector<Data *> data;

	Cluster(int dim)
	{
		mean.resize(dim);
		for(int i=0;i<dim;i++)
			mean(i) = pd.uniformRand(0.0,1.0);

		cov.resize(dim,dim);
		cov = Eigen::MatrixXd::Zero(dim,dim);
		for(int i=0;i<dim;i++)
			cov(i,i) = 0.5;

		dimension = dim;
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

		for(int i=0;i<dimension;i++)
			mean(i) = 0.0;

		for(int j=0;j<dimension;j++){
			for(auto d : data){
				mean[j] += d->normalized_data[j];
			}
			mean[j] /= data.size();
		}

		for(int i=0;i<dimension;i++)
			cov(i,i) = 0.02;
	}

	void print(void)
	{
		if(data.size() == 0)
			return;

		cerr << setprecision(2);
		for(int i=0;i<dimension;i++)
			cerr << mean(i) << ',';
		cerr << " cov: ";
		for(int i=0;i<dimension;i++)
			cerr << cov(i,i) << ',';

		cerr << " num: " << data.size() << endl;
	}

	int dimension;
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

	void remove(void)
	{
        	auto itrNewEnd = std::remove_if(c.begin(),c.end()
                                ,[](Cluster e)->bool
                                 { return e.data.size() == 0; });

        	c.erase(itrNewEnd,c.end());
	}
};

#endif
