CXX = g++

BOOST_VERSION=1.38.0

#PATHS
CURRENT_DIR	= $(shell pwd)/
PROPOSAL = $(CURRENT_DIR)/../build/
PROPOSALSRC = $(CURRENT_DIR)/../src/
#LOG4CPLUS = /usr/carguelles/programs/log4cplus/build/

SOURCES 	= $(wildcard *.cpp)

OBJECTS 	= $(SOURCES:.cpp=.o)

INCLUDE_PATH += -I/usr/local/include
INCLUDE_PATH += -I$(PROPOSALSRC)/public

#Libraries

#Compiler
CC 		= gcc
CXX 		= g++

#Dynamic Library
#Flags
CXX_FLAGS = $(INCLUDE_PATH) -I.  -O3 -g -fPIC -std=c++11

LD 		= g++
LD_FLAGS 	= -L/usr/local/lib/ -L/usr/lib
LD_FLAGS  +=-L$(PROPOSAL)/lib
#LD_FLAGS  +=-L$(LOG4CPLUS)/lib -L$(BOOST)/lib
LD_FLAGS  += -lPROPOSAL -lboost_program_options -llog4cplus -lboost_serialization
#LD_FLAGS += -lboost_histogram
#LD_FLAGS  += -lPROPOSAL -llog4cplus -lboost
LD_FLAGS  += `root-config --libs`

.PHONY: all clean

all: generate_taus.exe

generate_taus.exe: $(OBJECTS)
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

%.o:%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

clean:
	rm *.o *.exe
