#ifndef VECTORXF_H
#define VECTORXF_H

#include <iostream>

class VectorXf {
	
	int length;
	float *data;
	
	// indicates whether SStot has been calculated
	bool usesSStot;
	
	// SStot for calculating r-squared
	float SStot;
	
public:
	
	VectorXf(int length);
	
	float *getData();
	int getLength();
	
	void calculateSStot();
	
	float getSStot();
	
	~VectorXf();
};

#endif
