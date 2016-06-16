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
#include "Cluster.h"
using namespace std;

int dim = 0;

int pickCandidateCluster(vector<int> *bincount,double d = 0.0)
{
	vector<double> prob_list;
	const double alpha = 0.8;
	int n = accumulate(bincount->begin(),bincount->end(),0.0);
	int n_k = 0;
	int i = 0;
	for(auto bin : *bincount){
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
	Eigen::VectorXd mu = c.mean;//[0],c.mean[1]); 
	Eigen::VectorXd diff = d.normalized_data - mu;

	double det = c.cov.determinant();
	double a = 1.0/ (pow(2 * 3.151592,c.dimension*0.5) *  sqrt(det));
	double exp_part = -0.5 * diff.transpose() * info * diff;

	return a*exp(exp_part);
}

bool resampling(DataSet *ds, Clusters *cs, Data *d, vector<int> &bin)
{
	int org_c_num = (int)cs->c.size();
	bin[d->cluster_id]--;
	int candidate_cluster = pickCandidateCluster(&bin);
	bin[d->cluster_id]++;

	if(candidate_cluster == d->cluster_id){
		return false;
	}
	
	Cluster *old_cluster = &(cs->c[d->cluster_id]);
	double eval_new = 0.0;
	double eval_old = densityMultiNormal(*d,*old_cluster);
	Cluster c(dim);
	if(candidate_cluster == org_c_num){//新しいクラスタ
		eval_new = densityMultiNormal(*d,c);
	}else{//既存のクラスタ
		eval_new = densityMultiNormal(*d,cs->c[candidate_cluster]);
	}

	double acceptance = eval_new/eval_old;
	if(pd.uniformRand(0.0,1.0) >= acceptance)
		return false;

	bin[d->cluster_id]--;
	d->cluster_id = candidate_cluster;
	if(candidate_cluster == org_c_num){
		cs->c.push_back(c);
		bin.push_back(1);
	}else
		bin[candidate_cluster]++;
	return true;
}

void sweep(DataSet *ds,Clusters *cs)
{
	vector<int> bincount(cs->c.size(),0);
	for(auto x : ds->x){
		bincount[x.cluster_id]++;
	}

	int chance = 3;
	for(auto &target : ds->x){
		for(int i=0;i<chance;i++){
			resampling(ds,cs,&target,bincount);
		}
	}
	//どのクラスタに標本が幾つかる数える
	for(auto &c : cs->c){
		c.clear();
	}
	for(auto &d : ds->x){
		cs->c[d.cluster_id].regData(&d);
	}
	cs->calcParams();
	cerr << "----" << endl;
}

int main(int argc, char const* argv[])
{
	Clusters cs;
	DataSet ds;
	dim = ds.read();
	if(dim <= 0)
		exit(1);

	int sweep_num = 50;

	//最初のクラスタを作る。平均値は1軸ごとにガウス分布からサンプリング
	cs.c.push_back(Cluster(dim));
	for(auto &d : ds.x){
		cs.c[d.cluster_id].regData(&d);
	}
	cs.calcParams();
	cerr << "----" << endl;

	for(int k=0;k<sweep_num;k++){
		cerr << "sweep " << k << endl;
		sweep(&ds,&cs);
	}
	ds.print();
	
	exit(0);
}
