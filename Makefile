C      = g++
CFLAGS  = -Wall -pedantic -std=c++11 -O2
LDFLAGS = -lstdc++ -lm

all: heat clean

heat: heat.cpp TriDiagSolver.o Grid.o
		$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

TriDiagSolver.o: TriDiagSolver.cpp
		$(CC) -c $(CFLAGS) $<

Grid.o: Grid.cpp
		$(CC) -c $(CFLAGS) $<
clean:
		rm *.o

