#ifndef __DATA_SET_H___
#define __DATA_SET_H___

#include "Data.h"

class DataSet
{
public:
	vector<Data> x;
	void print(void);
	int read(void);
};

#endif
