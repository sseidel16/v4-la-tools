#ifndef CSMATRIX_H
#define CSMATRIX_H

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "FactorData.h"
#include "LocatingArray.h"
#include "VectorXf.h"

using namespace std;

struct FactorSetting {
	
	bool grouped;
	char factor_i;
	char index; // level index (1st level if part of group)
	char levelsInGroup; // levels in group (1 if not a group)
	
};
struct Mapping {
	int mappedTo;
	Mapping **mapping;
};
struct Path {
	int min;
	int max;
	
	Path *entryA;
	Path *entryB;
};

class CSCol {
public:
	// column data
	vector <float>dataVector;
	float *dataP; // pointer to the 1st vector element

	// number of contributing factors
	int factors;
	
	// column headers
	FactorSetting *setting;
	
};

class CSMatrix {
private:

	int rows;
	
	FactorData *factorData;
	LocatingArray *locatingArray;
	GroupingInfo **groupingInfo;
	
	vector <CSCol*>*data;	// m by n
	
	/* The Mapping structure helps map (factor + level) interactions to their
	appropriate column in the CS Matrix. First, a (factor + level) combination
	is mapped to an index using 'factorLevelMap'. The index for factor = 1,
	level = 2 can be found by: factorLevelMap[1][2]; This index represents a
	(factor + level) combination. A t-way interaction of these indeces is
	mapped to a column of the CS Matrix using 'mapping'. Indeces are placed in
	decreasing order. For example 10, 7, 2 represents a 3-way interaction and
	its column index would be found in the following way:
	mapping->mapping[10]->mapping[7]->mapping[2]->mappedTo;
	*/
	int **factorLevelMap;
	Mapping *mapping;
	
	void addRow(CSCol *csCol);
	
	string getFactorLevelName(int factor_i, int level_i);
	string getFactorString(FactorSetting setting);
	
	void addOneWayInteraction(int factor_i, char level_i, Factor *factor, char **levelMatrix, vector <float>&sumOfSquares);
	
	void addTWayInteractions(CSCol *csColA, int colBMax_i, int &col_i, int t,
		Mapping **mapping, vector <float>&sumOfSquares, GroupingInfo **groupingInfo, char **levelMatrix);
	int populateColumnData(CSCol *csCol, char **levelMatrix, int row_top, int row_len);
	void randomizePaths(CSCol **array, Path *path, int row_top, int k, long long int &score, list <Path*>*pathList);
	void randomizeRows(CSCol **backupArray, CSCol **array, long long int &csScore, int row_top, int row_len);
	void repopulateColumns(int setFactor_i, int setLevel_i, int row_top, int row_len);
	void repopulateColumns(int setFactor_i, int setLevel_i, int maxFactor_i, int t,
		Mapping *mapping, char **levelMatrix, int &lastCol_i, int row_top, int row_len);
	int getColIndex(CSCol *csCol);
	
	void swapColumns(CSCol **array, int col_i1, int col_i2);
	void swapRows(CSCol **array, int row_i1, int row_i2);
	void smartSort(CSCol **array, int sortedRows);
	void quickSort(CSCol **array, int min, int max, int row_top, int row_len);
	void rowSort(CSCol **array, int min, int max, int row_i, int row_len);
	void pathSort(CSCol **array, Path *path, int row_i, int &nPaths, list <Path*>*pathList);
	void deletePath(Path *path);
	void pathChecker(CSCol **array, Path *pathA, Path *pathB, int row_i, int k,
		long long int &score, FactorSetting *&settingToResample, long long int *rowContributions);
	int compare(CSCol *csCol1, CSCol *csCol2, int row_top, int row_len);
	long checkAdvanced(CSCol **array, int k, int min, int max, int row_top, int row_len, FactorSetting *&settingToResample);
	
	void addRow(CSCol **array, char *levelRow);
	void addRowFix(CSCol **array, long long int &csScore);
	long long int getArrayScore(CSCol **array);
	
public:
	CSMatrix(LocatingArray *locatingArray, FactorData *factorData);
	
	int getRows();
	int getCols();
	
	float getDistanceToCol(int col_i, float *residuals);
	float getProductWithCol(int col_i, float *residuals);
	
	CSCol *getCol(int col_i);
	
	string getColName(CSCol *csCol);
	
	void print();
	
	void reorderRows();
	void exactFix();
	void randomFix();
	
};

#endif
