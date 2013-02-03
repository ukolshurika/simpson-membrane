CC=g++

CFLAGS=-Wall -I. -I$(GTEST_DIR)/include -std=c++0x -ffast-math -msse2 -O2
LDFLAGS=$(GTEST_DIR)/libgtest.a -lpthread -lrt -pg
GLFLAGS=-lglut -lGL -lGLU

EXECUTABLE=membrane
UNITTESTS=unittests
ANIMATE=animation

$(EXECUTABLE): main.o membrane.o matrix_surface.o ideal_sliding.o free_deformation.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(UNITTESTS): unittests.o simpson_unittest.o membrane_unittest.o ideal_sliding_unittest.o matrix_surface_unittest.o membrane.o matrix_surface.o ideal_sliding.o free_deformation.o 
	$(CC) $^ -o $@ $(LDFLAGS)

$(ANIMATE): drawer2.cpp
	$(CC) $^ -o $@ $(GLFLAGS)

membrane.o: membrane.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

matrix_surface.o: matrix_surface.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

free_deformation.o: free_deformation.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

ideal_sliding.o: ideal_sliding.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

simpson_unittest.o: simpson_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

membrane_unittest.o: membrane_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

ideal_sliding_unittest.o: ideal_sliding_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

matrix_surface_unittest.o: matrix_surface_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

unittests.o: unittests.cc
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean

clean:
	rm -rf  *.o $(EXECUTABLE) $(UNITTESTS) $(ANIMATE)
