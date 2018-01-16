#include "LocatingArray.h"

#include <fstream>
#include <iostream>

LocatingArray::LocatingArray(string file) {
	int tempData;
	
	ifstream ifs(file.c_str(), ifstream::in);
	
	// read the first 2 lines of locating array
	ifs >> tests;
	ifs >> factors;
	
	cout << tests << ", " << factors << endl;
	factorGrouping = new GroupingInfo*[factors];
	levels = new char*[tests];
	
	// load the level counts for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		factorGrouping[factor_i] = new GroupingInfo();
		
		ifs >> tempData;
		factorGrouping[factor_i]->levels = tempData;
	}
	
	// load the grouping for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		// load grouped bool
		ifs >> tempData;
		factorGrouping[factor_i]->grouped = tempData;
		
		// load less than factor index
		ifs >> tempData;
		factorGrouping[factor_i]->lessThanFactor = tempData;
		
		// check grouping
		if (factorGrouping[factor_i]->grouped) {
			// allocate array for level grouping
			factorGrouping[factor_i]->levelGroups = new char[factorGrouping[factor_i]->levels];
		
			for (int level_i = 0; level_i < factorGrouping[factor_i]->levels; level_i++) {
				// load level group
				ifs >> tempData;
				factorGrouping[factor_i]->levelGroups[level_i] = tempData;
			}
		} else {
			factorGrouping[factor_i]->levelGroups = NULL;
		}
	}
	
	// load the tests now
	for (int test_i = 0; test_i < tests; test_i++) {
		levels[test_i] = new char[factors];
		
		// now load each factor level for this specific test
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			ifs >> tempData;
			levels[test_i][factor_i] = tempData;
		}
		
	}
	
	ifs.close();
	
}

GroupingInfo **LocatingArray::getGroupingInfo() {
	return factorGrouping;
}

char **LocatingArray::getLevelMatrix() {
	return levels;
}

int LocatingArray::getFactors() {
	return factors;
}

int LocatingArray::getTests() {
	return tests;
}
