CC=g++
CFLAGS=-g -Wall -pedantic -I.
LDFLAGS=
EXECUTABLE=membrane

$(EXECUTABLE): main.o membrane.o
	$(CC) $(LDFLAGS) $^ -o $@

membrane.o: membrane.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) $^ -c -o $@

.PHONY: clean

clean:
	rm -rf  *.o $(EXECUTABLE)
