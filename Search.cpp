#include <cmath>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <sys/types.h>

#include "CSMatrix.h"
#include "FactorData.h"
#include "Model.h"
#include "LocatingArray.h"
#include "VectorXf.h"

using namespace std;

// static member of the class
WorkSpace *Model::workSpace = NULL;

// loader section
VectorXf *loadResponseVector(string directory, string column, bool performLog) {
	struct dirent *dp;
	
	// open the directory
	DIR *dir = opendir(directory.c_str());
	if (dir == NULL) {
		cout << "Cannot open directory \"" << directory << "\"" << endl;
		exit(0);
	}
	
	int rows, cols, col_i;
	string header;
	
	string tempString;
	float tempFloat;
	
	float *data, *responseData;
	
	// the response vector and response count vector
	VectorXf *response = NULL;
	VectorXf *responseCount = NULL;
	
	while ((dp = readdir(dir)) != NULL) {
		
		// go through not hidden files
		if (dp->d_name[0] != '.') {
			
			string file = directory + "/" + dp->d_name;
			cout << "File: \"" << file << "\"" << endl;
			
			// initialize the file input stream
			ifstream ifs(file.c_str(), ifstream::in);
			
			// read the headers
			ifs >> rows;
			
			// allocate memory for response vector if necessary
			if (response == NULL) {
				// allocate memory
				response = new VectorXf(rows);
				responseCount = new VectorXf(rows);
				
				// grab data vectors
				data = response->getData();
				responseData = responseCount->getData();
				
				// initialize vectors
				for (int row_i = 0; row_i < rows; row_i++) {
					data[row_i] = 0;
					responseData[row_i] = 0;
				}
			}
			
			// make sure the rows match
			if (response->getLength() != rows) {
				cout << "Row mismatch: " << file << endl;
				cout << "Expected " << response->getLength() << " but received " << rows << endl;
				exit(0);
			}
			
			// find the relevant column
			col_i = -1;		// this is the column index we are lookin for
			string line = "";
			string value = "";
			getline(ifs, line); // terminating new line of rows
			getline(ifs, line); // headers
			stringstream colStream(line);
			for (int tcol_i = 0; colStream >> value; tcol_i++) {
				if (value == column) {
					col_i = tcol_i;
					break;
				}
			}

			// ensure we found the relevant column
			if (col_i == -1) {
				cout << "Could not find column \"" << column << "\" in \"" << file << "\"" << endl;
			} else {
				// read into the response vector struct
				data = response->getData();
				for (int row_i = 0; row_i < rows; row_i++) {
					getline(ifs, line);
					stringstream myStream(line);
					
					for (int tcol_i = 0; getline(myStream, value, '\t'); tcol_i++) {
						if (tcol_i == col_i) {
							if (value == "") {
								cout << "No value found" << endl;
							} else {
								tempFloat = atof(value.c_str());
								
								if (tempFloat == 0) {
									cout << "Warning: Float value response was 0 in " << file << " in row " << row_i << endl;
								}
								
								// process the data
								responseData[row_i] += 1;
								if (performLog) {
									data[row_i] += log(tempFloat);
								} else {
									data[row_i] += tempFloat;
								}
								
							}
								
						}
					}
					
				}
				
			}
			
			// close the file
			ifs.close();
			
		}
		
	}
	
	// average the responses
	data = response->getData();
	for (int row_i = 0; row_i < rows; row_i++) {
		data[row_i] /= responseData[row_i];
		//cout << data[row_i] << " from " << responseData[row_i] << " responses " << endl;
	}
	
	// delete response count vector
	delete responseCount;
	
	// close the directory
	closedir(dir);
	
	// calculate SStot for r-squared calculations
	response->calculateSStot();
	
	cout << "Loaded responses" << endl;
	
	return response;
}

void createModels(VectorXf *response, CSMatrix *csMatrix, int maxTerms, int models_n) {
	cout << "Creating Models..." << endl;
	Model::setupWorkSpace(response->getLength(), maxTerms);
	
	// work variables
	Model *model;		// current model we are working on
	
	Model **topModels = new Model*[models_n];
	Model **nextTopModels = new Model*[models_n];
	
	// populate initial top models
	topModels[0] = new Model(response, maxTerms, csMatrix);
	for (int model_i = 1; model_i < models_n; model_i++)
		topModels[model_i] = NULL;
	
	// struct for linked list of top columns of CS matrix
	struct ColDetails {
		int termIndex;
		float dotProduct;
		bool used;
		
	} *colDetails;
	
	// allocate memory for top columns
	colDetails = new ColDetails[csMatrix->getCols()];
	
	while (topModels[0] != NULL && topModels[0]->getTerms() < maxTerms) {
		// LOOP HERE
		
		// make sure all next top models are NULL
		for (int model_i = 0; model_i < models_n; model_i++) {
			nextTopModels[model_i] = NULL;
		}
		
		// grab the models from the topModels priority queue
		for (int model_i = 0; model_i < models_n; model_i++) {
			
			// we are done finding the next top models if we hit a NULL model
			if (topModels[model_i] == NULL) break;
			
			// grab the model from top models queue
			model = topModels[model_i];
			
			// grab the distances to columns in cs matrix
			for (int col_i = 0; col_i < csMatrix->getCols(); col_i++) {
				colDetails[col_i].dotProduct =
					csMatrix->getProductWithCol(col_i, model->getResiVec());
				colDetails[col_i].termIndex = col_i;
				colDetails[col_i].used = model->termExists(col_i);
			}
			
			// find the columns with the largest dot products (at most as many as we have models)
			for (int colsUsed = 0; colsUsed < models_n; colsUsed++) {
				
				float largestDotProduct = 0;
				int bestCol_i = -1;
				
				// find the unused column with smallest distance
				for (int col_i = 0; col_i < csMatrix->getCols(); col_i++) {
										
					// check if this column has larger dot product and should be marked as best
					if (!colDetails[col_i].used &&
						(bestCol_i == -1 || colDetails[col_i].dotProduct > largestDotProduct)) {
						bestCol_i = col_i;
						largestDotProduct = colDetails[col_i].dotProduct;
					}
					
				}
				
				// check if all columns have been used already
				if (bestCol_i == -1) {
					break;
				}
				
				// mark the term as used
				colDetails[bestCol_i].used = true;
				
				// add the term to the model temporarily
				model->addTerm(bestCol_i);
				model->leastSquares();
				
				// find a possible next top model to replace
				float worstRSquared = 0;
				int modelToReplace = -1;
				for (int model_i = 0; model_i < models_n; model_i++) {
					
					if (nextTopModels[model_i] == NULL) {
						// there is extra space (we found a NULL entry), insert immediately
						modelToReplace = model_i;
						break;
					} else if (nextTopModels[model_i]->getRSquared() == model->getRSquared()) {
						// matching r-squared means same model already added
						cout << "Duplicate Model!!!" << endl;
						modelToReplace = -1;
						break;
					} else if (nextTopModels[model_i]->getRSquared() < model->getRSquared() &&
							(modelToReplace == -1 || nextTopModels[model_i]->getRSquared() < worstRSquared)) {
						// better r-squared and the model to replace is the worst we have seen
						modelToReplace = model_i;
						worstRSquared = nextTopModels[model_i]->getRSquared();
					}
				}
				
				// replace the model if a better one was found
				if (modelToReplace != -1) {
					// make sure we deallocate any older next top model
					if (nextTopModels[modelToReplace] != NULL) {
						delete nextTopModels[modelToReplace];
						nextTopModels[modelToReplace] = NULL;
					}
					
					//cout << "Now we have a model with r-squared: " << model->getRSquared() << endl;
					
					// insert the new next top model
					nextTopModels[modelToReplace] = new Model(model);
				}
				
				// remove the temporarily added term
				model->removeTerm(bestCol_i);
				
			}
			
			// delete the model from the top models queue since it has been processed
			delete model;
			
		}
		
		// copy next top models to top models
		for (int model_i = 0; model_i < models_n; model_i++) {
			topModels[model_i] = nextTopModels[model_i];
		}
		
		// find the top model
		int topModel_i = -1;
		float bestRSquared = 0;
		for (int model_i = 0; model_i < models_n; model_i++) {
			if (topModels[model_i] != NULL &&
				(topModel_i == -1 || topModels[model_i]->getRSquared() > bestRSquared)) {
				topModel_i = model_i;
				bestRSquared = topModels[model_i]->getRSquared();
			}
		}
		
		if (topModel_i != -1) {
			cout << "Top Model (" << bestRSquared << "):" << endl;
			topModels[topModel_i]->printModelFactors();
		} else {
			cout << "No Model" << endl;
		}
		
	}
	
	cout << endl;
	cout << "Final Models Ranking: " << endl;
	for (int ranking = 1;; ranking++) {
		// find the top model
		int topModel_i = -1;
		float bestRSquared = 0;
		for (int model_i = 0; model_i < models_n; model_i++) {
			if (topModels[model_i] != NULL &&
				(topModel_i == -1 || topModels[model_i]->getRSquared() > bestRSquared)) {
				topModel_i = model_i;
				bestRSquared = topModels[model_i]->getRSquared();
			}
		}
		
		if (topModel_i != -1) {
			cout << "Model " << ranking << " (" << bestRSquared << "):" << endl;
			topModels[topModel_i]->printModelFactors();
			delete topModels[topModel_i];
			topModels[topModel_i] = NULL;
		} else {
			break;
		}
	}
	
}

int main(int argc, char **argv) {
	
	LocatingArray *array = new LocatingArray("LA.tsv");
	VectorXf *response = loadResponseVector("responses", "Response", false);
	FactorData *factorData = new FactorData("Factors.tsv");
	CSMatrix *matrix = new CSMatrix(array, factorData);
	createModels(response, matrix, 5, 5);
	
	cout << endl;
	cout << "Other Stuff" << endl;
	cout << "Response range: " << response->getData()[0] << " to " << response->getData()[response->getLength() - 1] << endl;
	cout << "Columns in CSMatrix: " << matrix->getCols() << endl;
	
//	matrix->print();
	
	if (argc > 1 && strcmp(argv[1], "memchk") == 0) {
		int exit;
		cin >> exit;
	}
	
	return 0;
	
}
