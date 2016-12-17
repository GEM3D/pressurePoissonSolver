C      = g++
include EigenDir.mk
CFLAGS  = -Wall -pedantic -std=c++11 -O2 -I $(EIGEN_DIR)
LDFLAGS = -lstdc++ -lm -llapack

all: heat clean

heat: heat.cpp 
		$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
#		rm *.o
