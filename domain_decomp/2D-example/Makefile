C      = g++
include EigenDir.mk
CFLAGS  = -Wall -pedantic -std=c++11 -O2 -I $(EIGEN_DIR)
LDFLAGS = -lstdc++ -lm -llapack

all: heat clean

heat: heat.cpp Domain.o
		$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

Domain.o: Domain.cpp
		$(CC) -c $(CFLAGS) $<

clean:
		rm *.o
