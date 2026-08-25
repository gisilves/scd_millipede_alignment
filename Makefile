CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

# ROOT configuration
ROOTCFLAGS := $(shell root-config --cflags 2>/dev/null)
ROOTLIBS   := $(shell root-config --libs 2>/dev/null)
ROOTLIBS   += -lMinuit

CXXFLAGS   += $(ROOTCFLAGS)

MILLESRC = ./src/Mille.cc

TARGETS = track_generator display_hits_root align

all: $(TARGETS)

track_generator: ./src/track_generator.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

display_hits_root: ./src/display_hits_root.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS)

align: ./src/alignment.cpp $(MILLESRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: all clean
