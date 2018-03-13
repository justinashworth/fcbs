CXX = g++

# compiler warning flags
WFLAGS = -Wall -W -Wextra -pedantic -Wno-attributes

# debug build flags
#DBFLAGS = -ggdb -pg
DBFLAGS = -ggdb -static
CXXFLAGS = $(WFLAGS) $(DBFLAGS)

# performance build flags
OFLAGS = -O3
CXXFLAGS = $(WFLAGS) $(OFLAGS)

EXE = fcbs
OBJECTFILES = main.o

# external libraries
LDLIBS = -lstdc++
#LDFLAGS = -pg

# build targets
all: $(EXE)

$(EXE): $(OBJECTFILES)
#	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE) $(LDFLAGS)
	$(CXX) $(OBJECTFILES) $(LDLIBS) -o $(EXE)

clean:
	-rm *.o $(EXE)

.PHONY: tags
tags:
	ctags *.cpp *.h
