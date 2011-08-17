XX = g++
OPTIONS = -std=c++0x -I. -Iinclude -I/home/feng/include -Wall  -O2 -march=native
#OPTIONS = -g -std=c++0x -I. -I/home/feng/include -Wall -O2 -march=native  -D_GLIBCXX_PARALLEL -fopenmp
SOURCES = test/test.cc

OBJECTS=$(SOURCES:.cc=.o)

EXECUTABLE=bin/test

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(XX) $(OBJECTS) -o $@  -lpthread

.cc.o:
	$(XX) -c $(OPTIONS) $< -o $@ 
	
clean:
	rm -f *.o
	rm -r test

