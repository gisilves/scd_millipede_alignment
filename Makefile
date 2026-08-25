CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -Itarget

# ROOT configuration
ROOTCFLAGS := $(shell root-config --cflags 2>/dev/null)
ROOTLIBS   := $(shell root-config --libs 2>/dev/null)
ROOTLIBS   += -lMinuit

CXXFLAGS   += $(ROOTCFLAGS)

MILLE_OBJ = ./target/Mille.o

TARGETS = track_generator display_hits_root align

all: $(TARGETS)

./target/Mille.o: ./target/Mille.cc ./target/Mille.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

track_generator: ./src/track_generator.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

display_hits_root: ./src/display_hits_root.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS)

align: ./src/alignment.cpp $(MILLE_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ROOTLIBS)

clean:
	rm -f $(TARGETS) ./src/*.o ./target/*.o

.PHONY: all clean