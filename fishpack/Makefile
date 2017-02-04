TARGET = test_hwscrt
FC = gfortran
FFLAGS = -O3
LFLAGS =
LINK = gfortran

OBJECTS = \
	./$(TARGET).o \
	./hwscrt_dble.o

all : $(OBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) -o $(TARGET)

clean :
	rm *.o $(TARGET)

### DO NOT remove this line - make depends on it ###
