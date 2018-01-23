#include "CSMatrix.h"

CSMatrix::CSMatrix(LocatingArray *locatingArray, FactorData *factorData) {
	
	srand(time(NULL));
	
	this->locatingArray = locatingArray;
	
	// assign factor data variable for use with grabbing column names
	this->factorData = factorData;
	
	// to be used when populating CSMatrix
	CSCol *csCol;
	int sum;
	
	// get number of factors in locating array
	int factors = locatingArray->getFactors();
	
	// get level counts (how many levels for each factor)
	groupingInfo = locatingArray->getGroupingInfo();
	
	// get the level matrix from locating array (this is the main data)
	char **levelMatrix = locatingArray->getLevelMatrix();
	
	// a row for every test
	rows = locatingArray->getTests();
	
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
	sum = 0;
	
	csCol->factors = 0;
	csCol->setting = new FactorSetting[0];
	
	// populate the intercept
	for (int row_i = 0; row_i < rows; row_i++) {
		addRow(csCol);
		csCol->dataP[row_i] = 1;
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
	
	cout << "Adding t-way interactions" << endl;
	addTWayInteractions(csCol, col_i - 1, col_i, locatingArray->getT(),
		mapping->mapping, sumOfSquares, groupingInfo, levelMatrix);
	
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

void CSMatrix::addRow(CSCol *csCol) {
	csCol->dataVector.push_back(0);
	csCol->dataP = &csCol->dataVector[0];
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
		for (int row_i = 0; row_i < rows; row_i++) {
			addRow(csCol);
		}
		
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
		sum = populateColumnData(csCol, levelMatrix, 0, rows);
		
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

int CSMatrix::populateColumnData(CSCol *csCol, char **levelMatrix, int row_top, int row_len) {
	int sum = 0;
	
	// populate every row
	bool rowData; // start with true, perform AND operation
	for (int row_i = row_top; row_i < row_top + row_len; row_i++) {
		
		rowData = true; // start with true, perform AND operation
		
		for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
			
			rowData &= levelMatrix[row_i][csCol->setting[setting_i].factor_i] >=
				csCol->setting[setting_i].index;
			
			rowData &= levelMatrix[row_i][csCol->setting[setting_i].factor_i] <
				csCol->setting[setting_i].index + csCol->setting[setting_i].levelsInGroup;
			
		}
		
		// AND operation
		csCol->dataP[row_i] = (rowData ? 1 : -1);
		
		// add to sum of squares
		sum += csCol->dataP[row_i] * csCol->dataP[row_i];
	}
	
	return sum;
}

void CSMatrix::exactFix() {
	
	// create a work array
	CSCol **array = new CSCol*[data->size()];
	for (int col_i = 0; col_i < data->size(); col_i++) {
		array[col_i] = data->at(col_i);
	}
	
	quickSort(array, 0, data->size() - 1, 0, rows);
	
	long long int csScore = getArrayScore(array);
	
	while (csScore > 0) {
		addRow(array, csScore);
	}
	
	delete [] array;
	
}

void CSMatrix::randomFix() {
	
	char **levelMatrix = locatingArray->getLevelMatrix();
	
	// check advanced
	CSCol **matrix = new CSCol*[data->size()];
	for (int col_i = 0; col_i < data->size(); col_i++) {
		matrix[col_i] = data->at(col_i);
	}
	
	FactorSetting *settingToResample = NULL;
	
	int k = 2;
	cout << "Checking advanced for " << k << " differences" << endl;
	
	long score = checkAdvanced(matrix, k, 0, data->size() - 1, 0, rows, settingToResample);
	cout << "Score: " << score << endl;
	int *oldLevels = new int[rows];
	/*
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
	*/
	
	delete [] oldLevels;
}

void CSMatrix::repopulateColumns(int setFactor_i, int row_top, int row_len) {
	
	int lastCol_i = -1;
	repopulateColumns(setFactor_i, locatingArray->getFactors() - 1, locatingArray->getT(),
		mapping, locatingArray->getLevelMatrix(), lastCol_i, row_top, row_len);
}

void CSMatrix::repopulateColumns(int setFactor_i, int maxFactor_i, int t,
		Mapping *mapping, char **levelMatrix, int &lastCol_i, int row_top, int row_len) {
	
	if (setFactor_i > maxFactor_i && mapping->mappedTo != lastCol_i) {
		populateColumnData(data->at(mapping->mappedTo), levelMatrix, row_top, row_len);
		
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
				mapping->mapping[factorLevelMap[factor_i][level_i]], levelMatrix, lastCol_i, row_top, row_len);
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
		subtractionResult = csCol->dataP[row_i] - residuals[row_i];
		distanceSum += subtractionResult * subtractionResult;
	}
	
	return distanceSum;
}

float CSMatrix::getProductWithCol(int col_i, float *residuals) {
	float dotSum = 0;
	CSCol *csCol = data->at(col_i);
	
	for (int row_i = 0; row_i < rows; row_i++) {
		dotSum += csCol->dataP[row_i] * residuals[row_i];
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
		addRow(csCol);
		
		// 1 or -1 depending on if the factor levels matched
		if (level_i == levelMatrix[row_i][factor_i])
			csCol->dataP[row_i] = 1;
		else
			csCol->dataP[row_i] = -1;
		
		// add to sum of squares
		sum += csCol->dataP[row_i] * csCol->dataP[row_i];
		
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
			cout << data->at(col_i)->dataP[row_i] << "\t";
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

// sort the array, given that some rows are already sorted
void CSMatrix::smartSort(CSCol **array, int sortedRows) {
	int min, max;
	min = 0;
	
	for (int col_i = 1; col_i < data->size(); col_i++) {
		// check if the streak ended
		if (compare(array[col_i - 1], array[col_i], 0, sortedRows) < 0) {
			max = col_i - 1;
			if (max - min > 0) quickSort(array, min, max, sortedRows, rows - sortedRows);
			min = col_i;
		}
	}
	
	// add the final streak
	max = data->size() - 1;
	if (max - min > 0) quickSort(array, min, max, sortedRows, rows - sortedRows);
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
	return memcmp(&csCol1->dataP[row_top], &csCol2->dataP[row_top], length);
}

long CSMatrix::checkAdvanced(CSCol **array, int k, int min, int max, int row_top, int row_len, FactorSetting *&settingToResample) {
	
	// no more differences to find or only 1 column
	if (k <= 0 || min >= max) {
		return 0;
	} else if (k > row_len) {
		
		cout << "Looking for " << k << " differences, but there are " << row_len << " rows left" << endl;
		
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

// FUNCTIONS FOR 'fixla'

void CSMatrix::addRow(CSCol **array, long long int &csScore) {
	// the maximum seconds without finding a better score
	int secsMax = 2;
	
	// total factors in locating array
	int factors = locatingArray->getFactors();
	
	long long int bestScore;
	long long int newScore;
	CSCol *csCol, *bestCol;
	
	int cols = getCols();
	
	// for backing up the order of array
	CSCol **backupArray = new CSCol*[cols];
	
	// increment rows
	rows++;
	
	cout << "The matrix now has " << rows << " rows" << endl;
	
	// add a row to each column of the CS matrix
	for (int col_i = 0; col_i < cols; col_i++) {
		addRow(array[col_i]);
	}
	
	// allocate memory for new row of locating array
	char *levelRow = new char[factors];
	char *oldLevelRow = new char[factors];
	char *newLevelRow = new char[factors];
	char *bestLevelRow = new char[factors];
	
	// track if factors of new row are finalized
	bool *finalized = new bool[factors];
	
	// add a random row non-finalized to locating array
	locatingArray->addLevelRow(levelRow);
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		finalized[factor_i] = false;
		
		levelRow[factor_i] = rand() % groupingInfo[factor_i]->levels;
	}
	
	// change the columns of the CS matrix affected by the factors
	for (int factor_i = 0; factor_i < factors; factor_i++) {
		repopulateColumns(factor_i, rows - 1, 1);
	}
	
	// smartly sort the array and score it
	smartSort(array, rows - 1);
	csScore = getArrayScore(array);
	
	cout << "Score after random non-finalized row: " << csScore << endl;
	
	// used to track duplicate columns in CS matrix
	bool duplicate = false;
	bool lastPairMatched = false;
	
	// continue adding finalizing factors while CS score decreases (improves)
	while (true) {
		
		// find the best column index (with the best score)
		bestScore = csScore;	// the original best is the original score
		bestCol = NULL;			// no best column yet
		
		// backup the sorted order of the array
		memcpy(backupArray, array, sizeof(CSCol*) * cols); // backup array
		
		// grab initial time
		struct timespec start;
		struct timespec finish;
		float elapsedTime;
		clock_gettime(CLOCK_REALTIME, &start);
		
		// go through all columns of CS matrix
		for (int col_i = 0; col_i < data->size() - 1; col_i++) {
			
			// check current time
			clock_gettime(CLOCK_REALTIME, &finish);
			// get elapsed seconds
			elapsedTime = (finish.tv_sec - start.tv_sec);
			// add elapsed nanoseconds
			elapsedTime += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
			if (elapsedTime > secsMax) break;
			
			// check for duplicate column
			if (compare(array[col_i], array[col_i + 1], 0, rows) >= 0) {
				duplicate = true;
				lastPairMatched = true;
			} else if (lastPairMatched) {
				duplicate = true;
				lastPairMatched = false;
			} else {
				duplicate = false;
			}
			
			// retrieve CSCol (and duplicate) from sorted array
			csCol = array[col_i];
			
			// if duplicate and the column has factors to change (column is not the INTERCEPT)
			if (duplicate && csCol->dataP[rows - 1] != 1 && csCol->factors > 0) {
				
				// grab duplicate column
				CSCol *csDup = array[col_i + 1];
				
				// the goal is to make the last row of csCol be a 1
				// all factors must be changed, check if this change conflicts with finalized factors
				bool changeAllowed = true;
				int factor_i;
				
				// go through all relevant factors for this column
				for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
					factor_i = csCol->setting[setting_i].factor_i;
					
					// backup the factor level
					oldLevelRow[factor_i] = levelRow[factor_i];
					
					// check if finalized yet
					if (!finalized[factor_i]) {
						// update factor to match the current column
						newLevelRow[factor_i] = csCol->setting[setting_i].index +
							(rand() % csCol->setting[setting_i].levelsInGroup);
					} else {
						// ensure this will actually make this column a 1
						changeAllowed &= (newLevelRow[factor_i] >= csCol->setting[setting_i].index &&
							newLevelRow[factor_i] < csCol->setting[setting_i].index + csCol->setting[setting_i].levelsInGroup);
					}
					
				}
				
				// check if this column change is allowed
				if (changeAllowed) {
					// CAREFUL!!!
					// col_i cannot be used after smartSort and before the rollback
					// because the array is reordered. Use csCol
					
					// try making the changes and see if the CS score gets smaller
					newScore = csScore;
					for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
						factor_i = csCol->setting[setting_i].factor_i;
						
						if (levelRow[factor_i] != newLevelRow[factor_i]) {
							// update columns
							levelRow[factor_i] = newLevelRow[factor_i];
							repopulateColumns(factor_i, rows - 1, 1);
						}
					}
					
					// smart sort and get new score
					smartSort(array, rows - 1);
					newScore = getArrayScore(array);
					
					// check if we have a better score
					if (newScore < bestScore) {
						// store the better score as the best, and save the column
						bestScore = newScore;
						bestCol = csCol;

						for (int factor_i = 0; factor_i < factors; factor_i++) {
							bestLevelRow[factor_i] = newLevelRow[factor_i];
						}
					}
					
					// rollback the change
					for (int setting_i = 0; setting_i < csCol->factors; setting_i++) {
						factor_i = csCol->setting[setting_i].factor_i;
						
						if (oldLevelRow[factor_i] != levelRow[factor_i]) {
							// rollback columns and get new score
							levelRow[factor_i] = oldLevelRow[factor_i];
							repopulateColumns(factor_i, rows - 1, 1);
						}
					}
					
					// restore old order
					memcpy(array, backupArray, sizeof(CSCol*) * data->size());
					
				}
			}
		}
		
		// check if we found a better score
		if (bestScore < csScore) {
			cout << "Found better score! " << bestScore << endl;
			
			// make the change to improve the score
			int factor_i;
			for (int setting_i = 0; setting_i < bestCol->factors; setting_i++) {
				factor_i = bestCol->setting[setting_i].factor_i;
				
				if (bestLevelRow[factor_i] != levelRow[factor_i]) {
					// update columns and get new score
					levelRow[factor_i] = bestLevelRow[factor_i];
					repopulateColumns(factor_i, rows - 1, 1);
				}
				
				// finalize the changes
				finalized[factor_i] = true;
			}
			
			// do a final smart sort and finalize the changed factors
			smartSort(array, rows - 1);
			csScore = getArrayScore(array);
			
		} else {
			break;
		}
		
	}
	
	cout << "Score after finalized row: " << csScore << endl;
	
	delete[] oldLevelRow;
	delete[] newLevelRow;
	delete[] bestLevelRow;
	delete[] backupArray;
	delete[] finalized;
	
}

long long int CSMatrix::getArrayScore(CSCol **array) {
	long long int streak = 0;
	long long int squaredSum = 0;
	
	for (int col_i = 0; col_i < data->size() - 1; col_i++) {
		
		// increment streak
		streak++;
		
		// check if the streak ended
		if (compare(array[col_i], array[col_i + 1], 0, rows) < 0) {
			squaredSum += streak * streak;
			streak = 0;
		} else if (compare(array[col_i], array[col_i + 1], 0, rows) > 0) {
			cout << "Mistake in array" << endl;
		}
	}
	
	// add the final streak
	streak++;
	squaredSum += streak * streak;
	
	return (squaredSum - data->size());
}
