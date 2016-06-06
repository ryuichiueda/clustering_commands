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
	vector<double> mean; //XXX: This should be Eigen::Vector2d
	Eigen::MatrixXd cov;
	vector<Data *> data;

	Cluster()
	{
		mean.push_back(pd.uniformRand(0.0,1.0));
		mean.push_back(pd.uniformRand(0.0,1.0));
		cov = Eigen::MatrixXd(2,2);
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

		mean.clear();
		mean.push_back(0.0);
		mean.push_back(0.0);
		for(int j=0;j<2;j++){
			for(int i=0;i<data.size();i++){
				mean[j] += data[i]->original_data[j];
			}
			mean[j] /= data.size();
		}
	}

	void print(void)
	{
		if(data.size() == 0)
			return;

		cout << mean[0] << ' ' << mean[1] << endl;
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
		log_new = logDensityMultiNormal(*d,c,cov);
	}else{//既存のクラスタ
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

	int sweep_num = 3;
	int chance = 3;

	//最初のクラスタを作る。平均値は1軸ごとにガウス分布からサンプリング
	cs.c.push_back(Cluster());

	for(int k=0;k<sweep_num;k++){
		cout << "sweep " << k << endl;
		for(int j=0;j<chance;j++){
			sweep(&ds,&cs);
		}

		//どのクラスタに標本が幾つかる数える
		for(auto &c : cs.c){
			c.clear();
		}
		for(auto &d : ds.x){
			cs.c[d.cluster_id].regData(&d);
		}
		cout << "%%" << endl;
		cs.calcParams();
		cout << "----" << endl;
		//ds.print();
	}
	
	exit(0);
}
