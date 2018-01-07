#ifndef CSMATRIX_H
#define CSMATRIX_H

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "FactorData.h"
#include "LocatingArray.h"

#define CS_COL_NONE		-1

using namespace std;

struct FactorSetting {
	
	bool grouped;
	char factor_i;
	char index; // level index (1st level if part of group)
	char levelsInGroup; // levels in group (1 if not a group)
	
};

class CSCol {
public:
	float *data;

	// column headers
	FactorSetting setting1;
	FactorSetting setting2;
	
};

class CSMatrix {
private:
	int rows;
	
	FactorData *factorData;
	
	vector <CSCol*>*data;	// m by n
	
	string getFactorLevelName(int factor_i, int level_i);
	string getFactorString(FactorSetting setting);
	
	void addOneWayInteraction(int factor_i, char level_i, Factor *factor, char **levelMatrix, vector <float>&sumOfSquares);
	
public:
	CSMatrix(LocatingArray *array, FactorData *factorData);
	
	int getRows();
	int getCols();
	
	float getDistanceToCol(int col_i, float *residuals);
	float getProductWithCol(int col_i, float *residuals);
	
	CSCol *getCol(int col_i);
	
	string getColName(int col_i);
	
	void print();
	
};

#endif
