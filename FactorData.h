#ifndef FACTORDATA_H
#define FACTORDATA_H

#include <string>

using namespace std;

class Factor {
public:

	// name of the factor
	string name;
	
	// level count
	int levels;
	
	// is numeric
	bool numeric;
	
	// names of the levels
	string *levelNames;
	
	// values of the levels (NULL if numeric is false)
	float *levelValues;
	
};

class FactorData {
private:
	// total factors (75ish for large dataset)
	int factorCount;
	
	// array of factor pointers
	Factor **factors;
	
public:
	
	FactorData(string file);
	
	Factor *getFactor(int factor_i);
	
	string getFactorName(int factor_i);
	
	string getFactorLevelName(int factor_i, int level_i);
	
};

#endif
