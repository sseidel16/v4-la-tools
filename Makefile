FLAG = 
CPLR = g++
SRCS = CSMatrix.cpp CSMatrix.h FactorData.cpp FactorData.h LocatingArray.cpp LocatingArray.h Model.cpp Model.h Noise.cpp Noise.h Occurrence.h Search.cpp VectorXf.cpp VectorXf.h
EXEC = Search

$(EXEC): $(SRCS)
	$(CPLR) -o $@ $^ $(FLAG)

clean:
	rm -f ./*.o $(EXEC)