CC=g++

CFLAGS=-g -Wall -I. -I$(GTEST_DIR)/include
LDFLAGS=$(GTEST_DIR)/libgtest.a -lpthread

EXECUTABLE=membrane
UNITTESTS=unittests

$(EXECUTABLE): main.o membrane.o
	$(CC) $^ -o $@ $(LDFLAGS)

$(UNITTESTS): unittests.o simpson_unittest.o
	$(CC) $^ -o $@ $(LDFLAGS)

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
	rm -rf  *.o $(EXECUTABLE) $(UNITTESTS)
