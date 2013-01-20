CC=g++-4.6

CFLAGS=-g -Wall -I. -I$(GTEST_DIR)/include -std=c++0x
LDFLAGS=$(GTEST_DIR)/libgtest.a -lpthread -lrt
GLFLAGS=-lglut -lGL -lGLU

EXECUTABLE=membrane
UNITTESTS=unittests
ANIMATE=animation

$(EXECUTABLE): main.o membrane.o matrix_surface.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(UNITTESTS): unittests.o simpson_unittest.o membrane_unittest.o membrane.o matrix_surface.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(ANIMATE): drawer2.cpp
	$(CC) $^ -o $@ $(GLFLAGS)

membrane.o: membrane.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

matrix_surface.o: matrix_surface.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

simpson_unittest.o: simpson_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

membrane_unittest.o: membrane_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

unittests.o: unittests.cc
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean

clean:
	rm -rf  *.o $(EXECUTABLE) $(UNITTESTS) $(ANIMATE)
