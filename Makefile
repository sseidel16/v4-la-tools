CPLR = g++
EXEC = Search
DEPS = CSMatrix.h FactorData.h LocatingArray.h Model.h Noise.h Occurrence.h VectorXf.h
FLAGS =

%.o: %.cpp $(DEPS)
	$(CPLR) -c $< $(FLAGS)

$(EXEC): ConstraintGroup.o CSMatrix.o FactorData.o LocatingArray.o Model.o Noise.o Search.o VectorXf.o VectorXf.h
	$(CPLR) -o $@ $^

clean:
	rm -f ./*.o $(EXEC)
