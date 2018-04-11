#include "LocatingArray.h"

LocatingArray::LocatingArray(string file) {
	int tempData;
	
	ifstream ifs(file.c_str(), ifstream::in);
	
	// read the first 2 lines of locating array
	ifs >> tests;
	ifs >> factors;
	
	t = 2;
	
	cout << tests << ", " << factors << endl;
	factorGrouping = new GroupingInfo*[factors];
	
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
		char *levelRow = new char[factors];
		
		// now load each factor level for this specific test
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			ifs >> tempData;
			levelRow[factor_i] = tempData;
		}
		
		addLevelRow(levelRow);
		tests--;
	}
	
	ifs.close();
	
}

void LocatingArray::addLevelRow(char *levelRow) {
	levels.push_back(levelRow);
	tests++;
}

char *LocatingArray::remLevelRow() {
	char *levelRow = levels.back();
	levels.pop_back();
	tests--;
	return levelRow;
}

GroupingInfo **LocatingArray::getGroupingInfo() {
	return factorGrouping;
}

char **LocatingArray::getLevelMatrix() {
	return &levels[0];
}

int LocatingArray::getFactors() {
	return factors;
}

int LocatingArray::getTests() {
	return tests;
}

int LocatingArray::getT() {
	return t;
}

void LocatingArray::writeToFile(string file) {
	
	ofstream ofs(file.c_str());
	
	// initial tests and factors
	ofs << tests << "\t" << factors << endl;
	
	// write the level counts for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		ofs << (int)factorGrouping[factor_i]->levels << "\t";
	}
	ofs << endl;
	
	// write the grouping for the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		// load grouped bool
		ofs << factorGrouping[factor_i]->grouped << "\t";
		
		// load less than factor index
		ofs << factorGrouping[factor_i]->lessThanFactor << "\t";
		
		// check grouping
		if (factorGrouping[factor_i]->grouped) {
		
			for (int level_i = 0; level_i < factorGrouping[factor_i]->levels; level_i++) {
				// write level group
				ofs << factorGrouping[factor_i]->levelGroups[level_i] << "\t";
			}
		}
		ofs << endl;
	}
	
	// write the tests now
	char **levelMatrix = getLevelMatrix();
	for (int test_i = 0; test_i < tests; test_i++) {
		// now write each factor level for this specific test
		for (int factor_i = 0; factor_i < factors; factor_i++) {
			ofs << (int)levelMatrix[test_i][factor_i] << "\t";
		}
		ofs << endl;
	}
	
	ofs.close();
	
}
