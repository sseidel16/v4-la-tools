#include "CSMatrix.h"

#include <cmath>

CSMatrix::CSMatrix(LocatingArray *array, FactorData *factorData) {
	// assign factor data variable for use with grabbing column names
	this->factorData = factorData;
	
	// to be used when populating CSMatrix
	CSCol *csCol;
	int sum;
	
	// get number of factors in locating array
	int factors = array->getFactors();
	
	// get level counts (how many levels for each factor)
	GroupingInfo **groupingInfo = array->getGroupingInfo();
	
	// get the level matrix from locating array (this is the main data)
	char **levelMatrix = array->getLevelMatrix();
	
	// a row for every test
	rows = array->getTests();
	
	// initialize the data column vector
	data = new vector<CSCol*>;
	
	// sum of squares of each column for normalization
	vector <float>sumOfSquares;
	
	// populate compressive sensing matrix column by column
	int col_i = 1;
	
	// add the INTERCEPT
	csCol = new CSCol;
	csCol->data = new float[rows];
	sum = 0;
	
	csCol->setting1.grouped = false;
	csCol->setting1.factor_i = CS_COL_NONE;
	csCol->setting1.index = 0;
	csCol->setting1.levelsInGroup = 1;
	csCol->setting2.grouped = false;
	csCol->setting2.factor_i = CS_COL_NONE;
	csCol->setting2.index = 0;
	csCol->setting2.levelsInGroup = 1;
	
	// populate the intercept
	for (int row_i = 0; row_i < rows; row_i++) {
		csCol->data[row_i] = 1;
		sum += 1;
	}
	
	// push into vector
	data->push_back(csCol);
	sumOfSquares.push_back(sum);
	
	// go over all 1-way interactions
	for (char factor_i = 0; factor_i < factors; factor_i++) {
		
		Factor *factor = factorData->getFactor(factor_i);
		
		for (char level = 0; level < groupingInfo[factor_i]->levels; level++) {
			
			addOneWayInteraction(factor_i, level, factor, levelMatrix, sumOfSquares);
			// increment column index
			col_i++;
			
		}
		
	}
	
	int lastOneWay_i = col_i;
	// go over all 2-way interactions by permuting all 1-way columns (watch out for groups)
	CSCol *csCol1, *csCol2;
	for (int col_i1 = 1; col_i1 < lastOneWay_i; col_i1++) {
		
		// load pointer to column
		csCol1 = data->at(col_i1);
		int colsInGroup1 = 1;
		char groupIndex1 = -1;
		char levelIndex1 = csCol1->setting1.index;
		
		// check if this column is grouped
		if (groupingInfo[csCol1->setting1.factor_i]->grouped) {
			groupIndex1 = groupingInfo[csCol1->setting1.factor_i]->levelGroups[levelIndex1];
			
			for (int level_i = levelIndex1 + 1;
				level_i < groupingInfo[csCol1->setting1.factor_i]->levels &&
				groupingInfo[csCol1->setting1.factor_i]->levelGroups[level_i] == groupIndex1; level_i++) {
				
				colsInGroup1++;
				col_i1++;
			}
			csCol1 = data->at(col_i1);
		}
		
		for (int col_i2 = col_i1 + 1; col_i2 < lastOneWay_i; col_i2++) {
			
			// load pointer to column
			csCol2 = data->at(col_i2);
			int colsInGroup2 = 1;
			char groupIndex2 = -1;
			char levelIndex2 = csCol2->setting1.index;
			
			// check if this column is grouped
			if (groupingInfo[csCol2->setting1.factor_i]->grouped) {
				groupIndex2 = groupingInfo[csCol2->setting1.factor_i]->levelGroups[levelIndex2];
				
				for (int level_i = levelIndex2 + 1;
					level_i < groupingInfo[csCol2->setting1.factor_i]->levels &&
					groupingInfo[csCol2->setting1.factor_i]->levelGroups[level_i] == groupIndex2; level_i++) {
					
					colsInGroup2++;
					col_i2++;
				}
				csCol2 = data->at(col_i2);
			}
			
			// check for valid 2-way interaction
			if (csCol1->setting1.factor_i != csCol2->setting1.factor_i) {
				
				// create new column for CS Matrix
				csCol = new CSCol;
				csCol->data = new float[rows];
				sum = 0;
				
				// set the headers from the combining columns
				csCol->setting1.grouped = (groupIndex1 != -1);			// is 1st factor grouped
				csCol->setting1.factor_i = csCol1->setting1.factor_i;	// 1st factor
				csCol->setting1.index = levelIndex1;					// set level of 1st factor
				csCol->setting1.levelsInGroup = colsInGroup1;			// set last level if group
				csCol->setting2.grouped = (groupIndex2 != -1);			// is 2nd factor grouped
				csCol->setting2.factor_i = csCol2->setting1.factor_i;	// 2nd factor
				csCol->setting2.index = levelIndex2;					// set level of 2nd factor
				csCol->setting2.levelsInGroup = colsInGroup2;			// set last level if group
				
				// populate every row
				for (int row_i = 0; row_i < rows; row_i++) {
					
					// multiply the two combining columns
					//csCol->data[row_i] = csCol1->data[row_i] * csCol2->data[row_i];
					char rowData1 = csCol1->data[row_i];
					char rowData2 = csCol2->data[row_i];
					
					if (groupIndex1 != -1) {
						// col_i1 has been incremented to the last column in the group so subtract offset
						for (int colOffset = 1; colOffset < colsInGroup1 && rowData1 == -1; colOffset++) {
							rowData1 = data->at(col_i1 - colOffset)->data[row_i];
						}
					}
					
					if (groupIndex2 != -1) {
						// col_i2 has been incremented to the last column in the group so subtract offset
						for (int colOffset = 1; colOffset < colsInGroup2 && rowData2 == -1; colOffset++) {
							rowData2 = data->at(col_i2 - colOffset)->data[row_i];
						}
					}
					
					csCol->data[row_i] = (rowData1 == 1 && rowData2 == 1 ? 1 : -1);
					
					// add to sum of squares
					sum += csCol->data[row_i] * csCol->data[row_i];
				}
				
				// push into vector
				data->push_back(csCol);
				sumOfSquares.push_back(sum);
				
				// increment column index
				col_i++;
			}
			
		}
	}
	
	cout << "Went over " << col_i << " columns" << endl;
	
	// perform sqaure roots on sum of squares
	for (int col_i = 0; col_i < getCols(); col_i++) {
		sumOfSquares[col_i] = sqrt(sumOfSquares[col_i]);
	}
	
	// normalize all columns
	for (int col_i = 0; col_i < getCols(); col_i++) {
		csCol = data->at(col_i);
		for (int row_i = 0; row_i < getRows(); row_i++) {
		//	csCol->data[row_i] = csCol->data[row_i] / sumOfSquares[col_i];
		}
		data->at(col_i) = csCol;
	}
	
	cout << "Finished constructing CS Matrix" << endl;
	
}

int CSMatrix::getRows() {
	return rows;
}

int CSMatrix::getCols() {
	return data->size();
}

CSCol *CSMatrix::getCol(int col_i) {
	return data->at(col_i);
}

float CSMatrix::getDistanceToCol(int col_i, float *residuals) {
	float distanceSum = 0, subtractionResult = 0;
	CSCol *csCol = data->at(col_i);
	
	for (int row_i = 0; row_i < rows; row_i++) {
		subtractionResult = csCol->data[row_i] - residuals[row_i];
		distanceSum += subtractionResult * subtractionResult;
	}
	
	return distanceSum;
}

float CSMatrix::getProductWithCol(int col_i, float *residuals) {
	float dotSum = 0;
	CSCol *csCol = data->at(col_i);
	
	for (int row_i = 0; row_i < rows; row_i++) {
		dotSum += csCol->data[row_i] * residuals[row_i];
	}
	
	return abs(dotSum);
}

string CSMatrix::getColName(int col_i) {
	
	ostringstream colName;
	
	// get factor indices and levels from CS matrix
	CSCol *csCol = data->at(col_i);
	
	if (csCol->setting1.factor_i == CS_COL_NONE && csCol->setting2.factor_i == CS_COL_NONE) {
		
		colName << "INTERCEPT";
		
	} else {
		
		// print 1st factor if it exists
		if (csCol->setting1.factor_i != CS_COL_NONE) {
			colName << getFactorString(csCol->setting1);
		}
		
		// print 2nd factor if it exists
		if (csCol->setting2.factor_i != CS_COL_NONE) {
			colName << " & ";
			colName << getFactorString(csCol->setting2);
		}
		
	}
	
	return colName.str();
	
}

string CSMatrix::getFactorString(FactorSetting setting) {
	ostringstream factorString;
	
	factorString << factorData->getFactorName(setting.factor_i) << "=";
	
	if (setting.grouped) {
		factorString << "GROUP(";
	}
	
	for (int index = setting.index; index < setting.index + setting.levelsInGroup; index++) {
		if (index != setting.index) factorString << "|";
		factorString << getFactorLevelName(setting.factor_i, index);
	}
	
	if (setting.grouped) {
		factorString << ")";
	}
	
	return factorString.str();
}

string CSMatrix::getFactorLevelName(int factor_i, int level_i) {
	return factorData->getFactorLevelName(factor_i, level_i);
}

void CSMatrix::addOneWayInteraction(int factor_i, char level_i, Factor *factor, char **levelMatrix, vector <float>&sumOfSquares) {
	
	// create new column for CS Matrix
	CSCol *csCol = new CSCol;
	csCol->data = new float[rows];
	int sum = 0;
	
	// assign the headers
	csCol->setting1.grouped = false;
	csCol->setting1.factor_i = factor_i;	// single factor
	csCol->setting1.index = level_i;		// level of single factor
	csCol->setting1.levelsInGroup = 1;
	csCol->setting2.grouped = false;
	csCol->setting2.factor_i = CS_COL_NONE;
	csCol->setting2.index = 0;
	csCol->setting2.levelsInGroup = 1;
	
	// populate every row
	for (int row_i = 0; row_i < rows; row_i++) {
		
		// 1 or -1 depending on if the factor levels matched
		if (level_i == levelMatrix[row_i][factor_i])
			csCol->data[row_i] = 1;
		else
			csCol->data[row_i] = -1;
		
		// add to sum of squares
		sum += csCol->data[row_i] * csCol->data[row_i];
		
	}
	
	// push into vector
	data->push_back(csCol);
	sumOfSquares.push_back(sum);
	
}

void CSMatrix::print() {
	for (int col_i = 0; col_i < getCols(); col_i++) {
		cout << getColName(col_i) << "\t";
	}
	cout << endl;
	
	for (int row_i = 0; row_i < rows; row_i++) {
		
		for (int col_i = 0; col_i < getCols(); col_i++) {
			cout << data->at(col_i)->data[row_i] << "\t";
		}
		
		cout << endl;
	}
}
