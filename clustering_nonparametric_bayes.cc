#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <vector>
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

Cluster getNewCluster(void)
{
	//double g0_variance = 4.0;
	Cluster c;
	c.mean.push_back(pd.uniformRand(0.0,1.0));
	c.mean.push_back(pd.uniformRand(0.0,1.0));
	return c;
}

double logDensityMultiNormal(Data &d, Cluster &c, Eigen::MatrixXd &cov)
{
	Eigen::MatrixXd info = cov.inverse();
	Eigen::Vector2d x(d.normalized_data[0],d.normalized_data[1]); 
	Eigen::Vector2d mu(c.mean[0],c.mean[1]); 
	Eigen::Vector2d diff = x - mu;

	double det = cov.determinant();
	double a = 1.0/ 2 * 3.151592 *  sqrt(det);
	double exp_part = -0.5 * diff.transpose() * info * diff;

	return log(a) + exp_part;
}

void resampling(DataSet *ds, Clusters *cs, Data *d)
{
	Eigen::MatrixXd cov(2,2);
	cov(0,0) = 0.2;
	cov(0,1) = 0.0;
	cov(1,0) = 0.0;
	cov(1,1) = 0.2;

	d->target = true;
	//どのクラスタに標本が幾つかる数える
	vector<int> bincount(cs->c.size(),0);
	for(auto x : ds->x){
		if(!x.target)	
			bincount[x.cluster_id]++;
	}
	d->target = false;

	int candidate_cluster = pickCandidateCluster(bincount);
	//int old_cluster = d->cluster_id;
	if(candidate_cluster == d->cluster_id)
		return;
	
	Cluster *old_cluster = &(cs->c[d->cluster_id]);
	double log_new = 0.0;
	double log_old = logDensityMultiNormal(*d,*old_cluster,cov);
	Cluster c;
	if(candidate_cluster == (int)cs->c.size()){//新しいクラスタ
		c = getNewCluster();
		log_new = logDensityMultiNormal(*d,c,cov);
	}else{//既存のクラスタ
		//c = cs->c[candidate_cluster];
		log_new = logDensityMultiNormal(*d,cs->c[candidate_cluster],cov);
	}

	double acceptance = exp(log_new - log_old);
	//double acceptance = 1.0 < alpha ? 1.0 : alpha;

	if(pd.uniformRand(0.0,1.0) < acceptance){
		d->cluster_id = candidate_cluster;
		if(candidate_cluster == (int)cs->c.size()){
			c.mean.clear();
			for(auto e : d->normalized_data)
				c.mean.push_back(e);
			cs->c.push_back(c);
		}
	}
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

	int sweep_num = 1;
	int chance = 3;

	//最初のクラスタを作る。平均値は1軸ごとにガウス分布からサンプリング
	cs.c.push_back(getNewCluster());

	for(int k=0;k<sweep_num;k++){
		for(int j=0;j<chance;j++){
			sweep(&ds,&cs);
		}

		//どのクラスタに標本が幾つかる数える
		vector<int> bincount(cs.c.size(),0);
		for(auto d : ds.x){
			if(!d.target)	
				bincount[d.cluster_id]++;
		}
		for(auto b : bincount){
			cout << b << endl;
		}
	}
	
	exit(0);
}
