CC=g++

CFLAGS=-g -Wall -I. -I$(GTEST_DIR)/include
LDFLAGS=$(GTEST_DIR)/libgtest.a -lpthread
GLFLAGS=-lglut -lGL -lGLU

EXECUTABLE=membrane
UNITTESTS=unittests
ANIMATE=animation

$(EXECUTABLE): main.o membrane.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(UNITTESTS): unittests.o simpson_unittest.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(ANIMATE): drawer2.cpp
	$(CC) $^ -o $@ $(GLFLAGS)

membrane.o: membrane.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

simpson_unittest.o: simpson_unittest.cc
	$(CC) $(CFLAGS) $^ -c -o $@

unittests.o: unittests.cc
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean

clean:
	rm -rf  *.o $(EXECUTABLE) $(UNITTESTS) $(ANIMATE)
