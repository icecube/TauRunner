CXX = g++

BOOST_VERSION=1.38.0

#PATHS
CURRENT_DIR	= $(shell pwd)/
PROPOSAL = /home/austin/PROPOSAL/build/
BOOST_HIST = /home/austin/histogram/
#LOG4CPLUS = /home/carguelles/programs/log4cplus/build/
#BOOST = /data/user/carguelles/programs/boost/build/

SOURCES 	= $(wildcard *.cpp)

OBJECTS 	= $(SOURCES:.cpp=.o)

INCLUDE_PATH	+= -I$(PROPOSAL)/include #-I$(BOOST)/include
INCLUDE_PATH += -I$(BOOST_HIST)/include

#Libraries

#Compiler
CC 		= gcc
CXX 		= g++

#Dynamic Library
#Flags
CXX_FLAGS = $(INCLUDE_PATH) -I.  -O3 -fPIC -std=c++11

LD 		= g++
LD_FLAGS 	= -L/usr/local/lib/ -L/usr/lib
LD_FLAGS  +=-L$(PROPOSAL)/lib
LD_FLAGS  +=-L$(BOOST_HIST)/lib
#LD_FLAGS  +=-L$(LOG4CPLUS)/lib -L$(BOOST)/lib
LD_FLAGS  += -lPROPOSAL -lboost_program_options -llog4cplus -lboost_serialization
#LD_FLAGS += -lboost_histogram
#LD_FLAGS  += -lPROPOSAL -llog4cplus -lboost
LD_FLAGS  += `root-config --libs`

.PHONY: all clean

all: generate_muons.exe

generate_muons.exe: $(OBJECTS)
	$(LD) $^ $(LIBS) $(LD_FLAGS) -o $@

%.o:%.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@

clean:
	rm *.o *.exe
