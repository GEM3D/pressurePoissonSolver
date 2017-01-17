C      = g++
CFLAGS  = -Wall -pedantic -std=c++11 -O2
LDFLAGS = -lstdc++ -lm -llapack

all: heat clean

heat: heat.cpp TriDiagSolver.o Domain.o
		$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

TriDiagSolver.o: TriDiagSolver.cpp
		$(CC) -c $(CFLAGS) $<

Domain.o: Domain.cpp
		$(CC) -c $(CFLAGS) $<
clean:
		rm *.o

