#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include "ProbDistributions.h"
using namespace std;

ProbDistributions pd;

class Data
{
public:
	vector<double> axis;
	int cluster_id;
	bool target;

	Data()
	{
		cluster_id = 0;
		target = false;
	}

	void print(void)
	{
		for(auto e : axis)
			cout << e << ' ';
		cout << endl;
	}
};

class DataSet
{
public:
	vector<Data> x;

	void print(void)
	{
		for(auto e : x)
			e.print();
	}

	bool read(void)
	{
		string line;
		while(1){
			getline(cin,line);
			if(!cin)
				break;
			stringstream ss;
			ss << line;

			Data d;
			double tmp;
			while(1){
				ss >> tmp;
				if(!ss)
					break;
				d.axis.push_back(tmp);
			}
			x.push_back(d);
		}
		return true;
	}
};

class Cluster
{
public:
	vector<double> mean;
};

class Clusters
{
public:
	vector<Cluster> c;
};

int pickCandidateCluster(vector<int> &bincount,double d = 0.0)
{
	vector<double> prob_list;
	const double alpha = 0.8;
	int n = accumulate(bincount.begin(),bincount.end(),0.0);
	int n_k = 0;
	int i = 0;
	for(auto bin : bincount){
		n_k += bin;
		prob_list.push_back(((double)n_k - d*(i+1))/(n + alpha));
	}

	double r = pd.uniformRand(0.0,1.0);
	int candidate_cluster = 0;
	for(auto p : prob_list){
		if(r < p)
			break;
		candidate_cluster++;
	}

	return candidate_cluster;
}

void resampling(DataSet *ds, Clusters *cs, Data *d)
{
	d->target = true;

	vector<int> bincount(cs->c.size(),0);

	//どのクラスタに標本が幾つかる数える
	for(auto x : ds->x){
		if(!x.target)	
			bincount[x.cluster_id]++;
	}

	int candidate_cluster = pickCandidateCluster(bincount);
	//int old_cluster = d->cluster_id;
	
	double alpha;
	if(candidate_cluster == (int)cs->c.size()){//新しいクラスタ
		Cluster c;
		c.mean.push_back(pd.normalRand(0.0,g0_variance));//なんで原点周りなのかとても怪しい
		c.mean.push_back(pd.normalRand(0.0,g0_variance));//たぶん適当な標本周りでサンプリング
	}else{//既存のクラスタ
	}

	d->target = false;
}

void sweep(DataSet *ds,Clusters *cs)
{
	for(auto &target : ds->x){
		resampling(ds,cs,&target);
	}
}

int main(int argc, char const* argv[])
{

	Clusters cs;

	DataSet ds;
	ds.read();

	int repeat_times = 1;//50;
	int R = 1;//3;
	double g0_variance = 4.0;

	//最初のクラスタを作る。平均値は1軸ごとにガウス分布からサンプリング
	Cluster c;
	c.mean.push_back(pd.normalRand(0.0,g0_variance));
	c.mean.push_back(pd.normalRand(0.0,g0_variance));
	cs.c.push_back(c);

	for(int k=0;k<repeat_times;k++){
		for(int j=0;j<R;j++){
			sweep(&ds,&cs);
		}
	}
	
	exit(0);
}
