FLAG = 
CPLR = g++
SRCS = FactorData.cpp FactorData.h LocatingArray.cpp LocatingArray.h CSMatrix.cpp CSMatrix.h Model.cpp Model.h Occurrence.h Search.cpp VectorXf.cpp VectorXf.h
EXEC = Search

$(EXEC): $(SRCS)
	$(CPLR) -o $@ $^ $(FLAG)

clean:
	rm -f ./*.o $(EXEC)