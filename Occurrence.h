#ifndef OCCURRENCE_H
#define OCCURRENCE_H

struct Occurrence {
	
	int *factorList;
	int factorList_n;
	
	int count;
	
	float magnitude;
	
	Occurrence *list;
	
};

#endif