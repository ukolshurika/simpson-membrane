CC=g++

SDIR=src
BDIR=bin

SOURCES=$(addprefix $(SDIR)/,main.cpp bound.cpp membrane.cpp matrix.cpp utils.cpp)
HEADERS=$(wildcard $(SDIR)/*.h)
OBJECTS=$(addprefix $(BDIR)/,$(notdir $(SOURCES:.cpp=.o)))

CFLAGS=-Wall -I. -I$(GTEST_DIR)/include -std=c++0x -msse2 -DDEBUG -g
LDFLAGS=$(GTEST_DIR)/libgtest.a -lpthread -lrt
GLFLAGS=-lglut -lGL -lGLU

EXECUTABLE=membrane
UNITTESTS=unittests
ANIMATE=animation
PLOT_DATA=plot

$(BDIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(BDIR)/$(PLOT_DATA): $(filter-out $(BDIR)/main.o,$(OBJECTS)) $(addprefix $(BDIR)/,plot.o)
	$(CC) $^ -o $@ $(LDFLAGS)

$(BDIR)/$(ANIMATE): $(SDIR)/drawer2.cpp
	$(CC) $^ -o $@ $(GLFLAGS)

$(addprefix $(BDIR)/,%.o): $(addprefix $(SDIR)/,%.cpp) $(HEADERS)
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: clean

clean:
	-rm $(OBJECTS) $(BDIR)/$(EXECUTABLE) $(BDIR)/$(ANIMATE) $(BDIR)/$(PLOT_DATA)
