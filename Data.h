#ifndef __DATA_H___
#define __DATA_H___

#include <iostream>
#include <vector>
using namespace std;

class Data
{
public:
	vector<double> original_data;
	vector<double> normalized_data;
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
		cout << endl;
	}
};

/*
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
				d.original_data.push_back(tmp);
			}
			x.push_back(d);
		}
		return true;
	}
};
*/

#endif
