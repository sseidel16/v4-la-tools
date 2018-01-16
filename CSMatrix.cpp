#include "CSMatrix.h"

// t-way interactions
CSMatrix::CSMatrix(LocatingArray *array, FactorData *factorData, int t) {
	
	srand(time(NULL));
	
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
	
	// initialize (factor + level)-to-index map
	factorLevelMap = new int*[factors];
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		factorLevelMap[factor_i] = new int[groupingInfo[factor_i]->levels];
	}
	
	// initialize index-to-column map (to be populated later)
	mapping = new Mapping;
	mapping->mappedTo = 0;
	mapping->mapping = NULL;
	
	// sum of squares of each column for normalization
	vector <float>sumOfSquares;
	
	// populate compressive sensing matrix column by column
	int col_i = 0;
	
	// Add the INTERCEPT
	csCol = new CSCol;
	csCol->data = new float[rows];
	sum = 0;
	
	csCol->factors = 0;
	csCol->setting = new FactorSetting[0];
	
	// populate the intercept
	for (int row_i = 0; row_i < rows; row_i++) {
		csCol->data[row_i] = 1;
		sum += 1;
	}
	
	// push into vector
	data->push_back(csCol);
	sumOfSquares.push_back(sum);
	col_i++;
	
	// go over all 1-way interactions
	for (char factor_i = 0; factor_i < factors; factor_i++) {
		
		Factor *factor = factorData->getFactor(factor_i);
		
		for (char level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
			
			addOneWayInteraction(factor_i, level_i, factor, levelMatrix, sumOfSquares);
			
			// set factor level map
			factorLevelMap[factor_i][level_i] = col_i - 1; // subtract 1 for INTERCEPT
			
			// increment column index
			col_i++;
			
		}
		
	}
	
	int lastOneWay_i = col_i;
	
	// initialize 1st level mapping
	mapping->mapping = new Mapping*[col_i - 1];
	
	//cout << "Adding t-way interactions" << endl;
	addTWayInteractions(csCol, col_i - 1, col_i, t, mapping->mapping, sumOfSquares, groupingInfo, levelMatrix);
	
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

void CSMatrix::addTWayInteractions(CSCol *csColA, int colBMax_i, int &col_i, int t,
		Mapping **mapping, vector <float>&sumOfSquares, GroupingInfo **groupingInfo, char **levelMatrix) {
	
	CSCol *csCol, *csColB, *csColC;
	int sum;
	
	// the offset is 1 for the INTERCEPT
	int colBOffset = 1;
	
	int colCMax_i = 0;
	
	// start after INTERCEPT, iterate over all main effects until colBMax_i
	for (int colB_i = 0; colB_i < colBMax_i; colB_i++) {
		
		// grab pointer to column
		csColB = data->at(colB_i + colBOffset);
		
		// the next 2 lines move colCMax_i to the 1st column with the same factor as csColB
		csColC = data->at(colCMax_i + colBOffset);
		if (csColB->setting[0].factor_i > csColC->setting[0].factor_i) colCMax_i = colB_i;
		
		// set mapping for this column
		mapping[colB_i] = new Mapping;
		mapping[colB_i]->mapping = (t > 1 ? new Mapping*[colCMax_i] : NULL);
		mapping[colB_i]->mappedTo = (csColA->factors > 0 ? col_i : colB_i + colBOffset);
		Mapping *groupMapping = mapping[colB_i];
		
		int colsInGroupB = 1;
		char groupIndexB = -1;
		char levelIndexB = csColB->setting[0].index;
		
		// check if this column is grouped
		if (groupingInfo[csColB->setting[0].factor_i]->grouped) {
			groupIndexB = groupingInfo[csColB->setting[0].factor_i]->levelGroups[levelIndexB];
			
			for (int level_i = levelIndexB + 1;
				level_i < groupingInfo[csColB->setting[0].factor_i]->levels &&
				groupingInfo[csColB->setting[0].factor_i]->levelGroups[level_i] == groupIndexB; level_i++) {
				
				colsInGroupB++;
				colB_i++;
				
				// set to group mapping unless main effect (no grouping on main effects)
				if (csColA->factors > 0) {
					mapping[colB_i] = groupMapping;
				} else {
					mapping[colB_i] = new Mapping;
					mapping[colB_i]->mapping = groupMapping->mapping;
					mapping[colB_i]->mappedTo = colB_i + colBOffset;
				}
			}
			
			csColB = data->at(colB_i + colBOffset);
		}
		
		// create new column for CS Matrix
		csCol = new CSCol;
		csCol->data = new float[rows];
		
		// set the headers from the combining columns
		csCol->factors = csColA->factors + 1;
		csCol->setting = new FactorSetting[csCol->factors];
		
		// copy previous factors
		for (int setting_i = 0; setting_i < csColA->factors; setting_i++) {
			csCol->setting[setting_i] = csColA->setting[setting_i];
		}
		
		// assign new setting
		csCol->setting[csColA->factors].grouped = (groupIndexB != -1);			// is 1st factor grouped
		csCol->setting[csColA->factors].factor_i = csColB->setting[0].factor_i;	// 1st factor
		csCol->setting[csColA->factors].index = levelIndexB;					// set level of 1st factor
		csCol->setting[csColA->factors].levelsInGroup = colsInGroupB;			// set last level if group
		
		// populate the actual column data of CS matrix
		sum = populateColumnData(csCol, levelMatrix);
		
		if (csCol->factors > 1) {
			// push into vector
			data->push_back(csCol);
			sumOfSquares.push_back(sum);
			
			col_i++;
			
		}
		
		if (t > 1) {
			// recursive call
			addTWayInteractions(csCol, colCMax_i, col_i, t - 1,
				mapping[colB_i]->mapping, sumOfSquares, groupingInfo, levelMatrix);
		}
	}
	
}

int CSMatrix::populateColumnData(CSCol *csCol, char **levelMatrix) {
	int sum = 0;
	
	// populate every row
	bool rowData; // start with true, perform AND operation
	for (int row_i = 0; row_i < rows; row_i++) {
		
		rowData = true; // start with true, perform AND operation
		
		for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
			
			rowData &= levelMatrix[row_i][csCol->setting[setting_i].factor_i] >=
				csCol->setting[setting_i].index;
			
			rowData &= levelMatrix[row_i][csCol->setting[setting_i].factor_i] <
				csCol->setting[setting_i].index + csCol->setting[setting_i].levelsInGroup;
			
		}
		
		// AND operation
		csCol->data[row_i] = (rowData ? 1 : -1);
		
		// add to sum of squares
		sum += csCol->data[row_i] * csCol->data[row_i];
	}
	
	return sum;
}

void CSMatrix::randomFix(LocatingArray *array, int t) {
	
	GroupingInfo **groupingInfo = array->getGroupingInfo();
	char **levelMatrix = array->getLevelMatrix();
	
	// check advanced
	CSCol **matrix = new CSCol*[data->size()];
	for (int col_i = 0; col_i < data->size(); col_i++) {
		matrix[col_i] = data->at(col_i);
	}
	
	FactorSetting *settingToResample = NULL;
	
	int k = 3;
	cout << "Checking advanced for " << k << " differences" << endl;
	
	long score = checkAdvanced(matrix, k, 0, data->size() - 1, 0, rows, settingToResample);
	
	int *oldLevels = new int[rows];
	
	while (true) {
		
		int row_i = rand() % rows;
		cout << "Resampling setting: " << getFactorString(*settingToResample) << " row: " << row_i << endl;
		
		oldLevels[row_i] = levelMatrix[row_i][settingToResample->factor_i];
		
		levelMatrix[row_i][settingToResample->factor_i] =
			rand() % groupingInfo[settingToResample->factor_i]->levels;
		
		int lastCol_i = -1;
		repopulateColumns(settingToResample->factor_i, array->getFactors() - 1, 2,
			mapping, groupingInfo, levelMatrix, lastCol_i);
		
		long newScore = checkAdvanced(matrix, k, 0, data->size() - 1, 0, rows, settingToResample);
		
		if (newScore > score) {
			cout << "Rollback" << endl;
			levelMatrix[row_i][settingToResample->factor_i] = oldLevels[row_i];
			
			int lastCol_i = -1;
			repopulateColumns(settingToResample->factor_i, array->getFactors() - 1, 2,
				mapping, groupingInfo, levelMatrix, lastCol_i);
			
			score = checkAdvanced(matrix, k, 0, data->size() - 1, 0, rows, settingToResample);
		} else {
			score = newScore;
			
			cout << "Score updated: " << score << endl;
		}
		
	}
	
}

void CSMatrix::repopulateColumns(int setFactor_i, int maxFactor_i, int t,
		Mapping *mapping, GroupingInfo **groupingInfo, char **levelMatrix, int &lastCol_i) {
	
	if (setFactor_i > maxFactor_i && mapping->mappedTo != lastCol_i) {
//		cout << "Repopulating " << getColName(getCol(mapping->mappedTo)) << ": " << mapping->mappedTo << endl;
		populateColumnData(data->at(mapping->mappedTo), levelMatrix);
		lastCol_i = mapping->mappedTo;
	}
	
	if (t == 0) return;
	
	int minFactor_i;
	int minLevel_i, maxLevel_i;
	
	if (setFactor_i > maxFactor_i) {
		minFactor_i = t - 1;
	} else {
		minFactor_i = setFactor_i;
		
		// if this is the last interaction, make sure setFactor shows up
		if (t == 1) maxFactor_i = minFactor_i;
	}
	
	for (int factor_i = minFactor_i; factor_i <= maxFactor_i; factor_i++) {
		
		for (int level_i = 0; level_i < groupingInfo[factor_i]->levels; level_i++) {
			repopulateColumns(setFactor_i, factor_i - 1, t - 1,
				mapping->mapping[factorLevelMap[factor_i][level_i]], groupingInfo, levelMatrix, lastCol_i);
		}
		
	}
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

string CSMatrix::getColName(CSCol *csCol) {
	
	ostringstream colName;
	
	if (csCol->factors == 0) colName << "INTERCEPT";
	
	// add all factor strings
	for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
		colName << (setting_i == 0 ? "" : " & ") << getFactorString(csCol->setting[setting_i]);
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
	csCol->factors = 1; // 1 contributing factor (1-way interaction)
	csCol->setting = new FactorSetting[1];
	csCol->setting[0].grouped = false;
	csCol->setting[0].factor_i = factor_i;	// single factor
	csCol->setting[0].index = level_i;		// level of single factor
	csCol->setting[0].levelsInGroup = 1;
	
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
	cout << "CSMatrix:" << endl;
	for (int col_i = 0; col_i < getCols(); col_i++) {
		cout << getColName(data->at(col_i)) << "\t";
	}
	cout << endl;
	
	for (int row_i = 0; row_i < rows; row_i++) {
		
		for (int col_i = 0; col_i < getCols(); col_i++) {
			cout << data->at(col_i)->data[row_i] << "\t";
		}
		
		cout << endl;
	}
	
	// verify that the column mappings are correct
	for (int col_i = 0; col_i < getCols(); col_i++) {
		if (col_i != getColIndex(data->at(col_i))) {
			cout << "Invalid Mapping!!" << endl;
		}
	}
	
}

int CSMatrix::getColIndex(CSCol *csCol) {
	
	Mapping *m = mapping;
	for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
		if (m == NULL) return -1;
		m = m->mapping[factorLevelMap[csCol->setting[setting_i].factor_i][csCol->setting[setting_i].index]];
	}
	return m->mappedTo;
	
}

void CSMatrix::swapColumns(CSCol **array, int col_i1, int col_i2) {
	CSCol *tempCol = array[col_i1];
	array[col_i1] = array[col_i2];
	array[col_i2] = tempCol;
}

void CSMatrix::quickSort(CSCol **array, int min, int max, int row_top, int row_len) {	
	if (min == max) return;
	
	int tempMin = min;
	int tempMax = max;
	CSCol *midCol = array[(min + max) / 2];
	
	while (true) {
		
		while (tempMin <= tempMax && tempMin < max && compare(array[tempMin], midCol, row_top, row_len) < 0) {
			tempMin++;
		}
		
		while (tempMin <= tempMax && tempMax > min && compare(array[tempMax], midCol, row_top, row_len) > 0) {
			tempMax--;
		}
		
		if (tempMin == tempMax) {
			if (tempMin == min) tempMin++;
			else tempMax--;
			break;
		} else if (tempMin <= tempMax) {
			swapColumns(array, tempMin, tempMax);
			tempMin++;
			tempMax--;
		} else {
			break;
		}
		
	}
	
	quickSort(array, min, tempMax, row_top, row_len);
	quickSort(array, tempMin, max, row_top, row_len);
}

int CSMatrix::compare(CSCol *csCol1, CSCol *csCol2, int row_top, int row_len) {
	int length = row_len * sizeof(float);
	return memcmp(&csCol1->data[row_top], &csCol2->data[row_top], length);
}

long CSMatrix::checkAdvanced(CSCol **array, int k, int min, int max, int row_top, int row_len, FactorSetting *&settingToResample) {
	
	// no more differences to find or only 1 column
	if (k <= 0 || min >= max) {
		return 0;
	} else if (row_len == 0) {
		
//		cout << "Looking for " << k << " differences, but there are " << row_len << " rows left" << endl;
		
		CSCol *csCol1 = array[min];
		CSCol *csCol2 = array[max];
		
		if (settingToResample == NULL) {
			// choose random factor+level for resampling
			CSCol *csColToResample = NULL;
			
			if (csCol1->factors != 0 && csCol2->factors != 0) {
				// neither has 0 factors, choose randomly
				if (rand() % 2) {
					csColToResample = csCol1;
				} else {
					csColToResample = csCol2;
				}
			} else if (csCol1->factors != 0) {
				// csCol1 has factors so choose it (csCol2 has no factors)
				csColToResample = csCol1;
			} else if (csCol2->factors != 0) {
				// csCol2 has factors so choose it (csCol1 has no factors)
				csColToResample = csCol2;
			}
			
			// chose a random factor setting to resample
			if (csColToResample != NULL) {
				settingToResample = &csColToResample->setting[rand() % csColToResample->factors];
			}
		}
/*
		cout << "Factor1: " << getColName(csCol1) << endl;
		cout << "Factor2: " << getColName(csCol2) << endl;
		
		for (int row_i = 0; row_i < rows; row_i++) {
			if (csCol1->data[row_i] != csCol2->data[row_i]) {
				cout << "Difference at row " << row_i << endl;
			}
		}
*/
		int cols = max - min + 1;
		return (long) pow((double) (cols * cols), (double) k);
	} else {
		
		long score = 0;
		
//		cout << "Sorting top row" << endl;
//		cout << "Data size is " << data->size() << endl;
		quickSort(array, min, max, row_top, 1);
		
		int divide;
		for (divide = min; divide < max; divide++) {
			if (compare(array[divide], array[divide + 1], row_top, 1) != 0) {
				break;
			}
		}
		for (int col_i = divide + 1; col_i < max; col_i++) {
			if (compare(array[col_i], array[col_i + 1], row_top, 1) != 0) {
				cout << "Integrity issue" << endl;
			}
		}
//		cout << "Not the same at " << divide << endl;
		
		// check the left side
		score += checkAdvanced(array, k, min, divide, row_top + 1, row_len - 1, settingToResample);
		
		// check the right side
		score += checkAdvanced(array, k, divide + 1, max, row_top + 1, row_len - 1, settingToResample);
		
		// check both
		score += checkAdvanced(array, k - 1, min, max, row_top + 1, row_len - 1, settingToResample);
		
		return score;
	}
}

