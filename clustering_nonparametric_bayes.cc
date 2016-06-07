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

		mean.clear();
		mean.push_back(0.0);
		mean.push_back(0.0);
		for(int j=0;j<2;j++){
			for(int i=0;i<data.size();i++){
				mean[j] += data[i]->normalized_data[j];
			}
			mean[j] /= data.size();
		}

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
			cov(0,0) = 0.02;
			cov(1,1) = 0.02;
		}
	}

	void print(void)
	{
		if(data.size() == 0)
			return;

		cout << setprecision(2);
		cout << mean[0] << ',' << mean[1] << " cov: " << cov(0,0) << ',' << cov(1,1)
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

double densityMultiNormal(Data &d, Cluster &c)
{
	Eigen::MatrixXd info = c.cov.inverse();
	Eigen::Vector2d x(d.normalized_data[0],d.normalized_data[1]); 
	Eigen::Vector2d mu(c.mean[0],c.mean[1]); 
	Eigen::Vector2d diff = x - mu;

	double det = c.cov.determinant();
	double a = 1.0/ (2 * 3.151592 *  sqrt(det));
	double exp_part = -0.5 * diff.transpose() * info * diff;

	return a*exp(exp_part);
}

bool resampling(DataSet *ds, Clusters *cs, Data *d)
{
	d->target = true;
	//どのクラスタに標本が幾つかる数える
	vector<int> bincount(cs->c.size(),0);
	for(auto x : ds->x){
		if(!x.target)	
			bincount[x.cluster_id]++;
	}
	d->target = false;

	int candidate_cluster = pickCandidateCluster(bincount);
	if(candidate_cluster == d->cluster_id)
		return false;
	
	Cluster *old_cluster = &(cs->c[d->cluster_id]);
	double eval_new = 0.0;
	double eval_old = densityMultiNormal(*d,*old_cluster);
	Cluster c;
	if(candidate_cluster == (int)cs->c.size()){//新しいクラスタ
		eval_new = densityMultiNormal(*d,c);
	}else{//既存のクラスタ
		eval_new = densityMultiNormal(*d,cs->c[candidate_cluster]);
/*
		if(cs->c[d->cluster_id].data.size() == 1){
			cout << "old";
			cs->c[d->cluster_id].print();
			cout << "new";
			cs->c[candidate_cluster].print();
		}
*/
	}

	double acceptance = eval_new/eval_old;
/*
	if(cs->c[d->cluster_id].data.size() == 1)
		cout << eval_old << " " << eval_new << " " << acceptance << endl;
*/
	if(pd.uniformRand(0.0,1.0) >= acceptance)
		return false;

	d->cluster_id = candidate_cluster;
	if(candidate_cluster == (int)cs->c.size()){
		c.mean.clear();
		for(auto e : d->normalized_data)
			c.mean.push_back(e);
		cs->c.push_back(c);
	}
	return false;
}

void sweep(DataSet *ds,Clusters *cs)
{
	int chance = 3;
	for(auto &target : ds->x){
		for(int i=0;i<chance;i++)
			if(resampling(ds,cs,&target))
				break;
	}
}

int main(int argc, char const* argv[])
{
	Clusters cs;

	DataSet ds;
	ds.read();

	int sweep_num = 50;

	//最初のクラスタを作る。平均値は1軸ごとにガウス分布からサンプリング
	cs.c.push_back(Cluster());
	for(auto &d : ds.x){
		cs.c[d.cluster_id].regData(&d);
	}
	cs.calcParams();
	cout << "----" << endl;

	for(int k=0;k<sweep_num;k++){
		cout << "sweep " << k << endl;
		sweep(&ds,&cs);

		//どのクラスタに標本が幾つかる数える
		for(auto &c : cs.c){
			c.clear();
		}
		for(auto &d : ds.x){
			cs.c[d.cluster_id].regData(&d);
		}
		cs.calcParams();
		cout << "----" << endl;
	}
	
	exit(0);
}
