#ifndef LOCATINGARRAY_H
#define LOCATINGARRAY_H

#include <string>

using namespace std;

struct GroupingInfo {
	int levels;	// total levels for this factor
	
	bool grouped;
	int lessThanFactor;
	char *levelGroups;
};

class LocatingArray {
private:
	GroupingInfo **factorGrouping;
	
	char **levels;	// pointer to array of test levels
	
	int tests;		// count of tests in locating array
	int factors;	// count of factors in locating array
	
public:
	LocatingArray(string file);
	
	GroupingInfo **getGroupingInfo();
	
	char **getLevelMatrix();
	
	int getFactors();
	int getTests();
};

#endif