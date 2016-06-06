#include "DataSet.h"
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;

void DataSet::print(void)
{
	for(auto e : x)
		e.print();
}

bool DataSet::read(void)
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
			d.original_data.push_back(tmp);
		}
		x.push_back(d);
	}

	int dim = (int)x[0].original_data.size();
	for(int i=0;i<dim;i++){
		double min = 1e100;
		double max = -1e100;

		for(Data d : x){
			if(min > d.original_data[i])
				min = d.original_data[i];
			else if(max < d.original_data[i])
				max = d.original_data[i];
		}
		if(max == min)
			max += 1.0; // all values of data becomes the value of min
		for(Data &d : x){
			d.normalized_data.push_back((d.original_data[i]-min)/(max - min));
		}
	}
	return true;
}
