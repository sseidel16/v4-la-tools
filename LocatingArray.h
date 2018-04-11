#ifndef LOCATINGARRAY_H
#define LOCATINGARRAY_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

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
	
	vector <char*>levels;	// pointer to array of test levels (the main locating array)
	
	int tests;		// count of tests in locating array
	int factors;	// count of factors in locating array
	
	int t;			// covers t-way interactions
	
public:
	LocatingArray(string file);
	
	void addLevelRow(char *levelRow);
	char *remLevelRow();
	
	GroupingInfo **getGroupingInfo();
	
	char **getLevelMatrix();
	
	int getFactors();
	int getTests();
	
	int getT();
	
	void writeToFile(string file);
};

#endif