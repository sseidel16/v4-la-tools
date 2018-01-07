
#include <fstream>

#include "FactorData.h"

FactorData::FactorData(string file) {
	ifstream ifs(file.c_str(), ifstream::in);
	
	// grab the total number of factors
	ifs >> factorCount;
	
	// allocate necessary memory
	factors = new Factor*[factorCount];
	
	// read in line by line
	for (int factor_i = 0; factor_i < factorCount; factor_i++) {
		// allocate the factor
		Factor *factor = new Factor;
		
		// read the factor name
		ifs >> factor->name;
		
		// read the number of levels
		ifs >> factor->levels;
		
		// read whether the factor has numeric levels
		ifs >> factor->numeric;
		
		// allocate memory for each level
		factor->levelNames = new string[factor->levels];
		
		// read individual level names
		for (int level_i = 0; level_i < factor->levels; level_i++) {
			ifs >> factor->levelNames[level_i];
		}
		
		// read the level values if they are numeric
		if (factor->numeric) {
			// allocate memory for each level value
			factor->levelValues = new float[factor->levels];
			
			// read individual level values
			for (int level_i = 0; level_i < factor->levels; level_i++) {
				ifs >> factor->levelValues[level_i];
			}
		} else {
			factor->levelValues = NULL;
		}
		
		factors[factor_i] = factor;
	}
}

Factor *FactorData::getFactor(int factor_i) {
	return factors[factor_i];
}

string FactorData::getFactorName(int factor_i) {
	return factors[factor_i]->name;
}

string FactorData::getFactorLevelName(int factor_i, int level_i) {
	return factors[factor_i]->levelNames[level_i];
}
