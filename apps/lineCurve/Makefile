#######################################
# MAKEFILE for this directory         #
#######################################

# determine modus of compilation: YES,[unset implies release]
ifneq ($(origin DEBUG),environment)
  DEBUG = 
endif


# compile/link executables and flags
SHELL = /bin/bash
CXX   = g++
LD    = g++

# includes
INCLUDES = -I ../../tools/

# libraries
LIBS = 

# common flags
CPPFLAGS = $(INCLUDES)

# set flags depending on mode
ifeq ($(DEBUG),YES) 
  CPPFLAGS += -g3 -Wall -Wextra -ansi -pedantic -Wno-unused-parameter -Wno-variadic-macros -Wconversion -std=c++0x 
  LDFLAGS  += -g3 
else
  CPPFLAGS += -g -O3 
  LDFLAGS  += -g -O3 
endif

# the final target
TARGET = lineCurve

# find the sources and object files
CPPFILES = $(wildcard *.cpp)
CPPOBJS = $(CPPFILES:.cpp=.o)

#define sources and objects
SRCS = $(CPPFILES)
OBJS = $(CPPOBJS)

# for cleaning up
RM = -rm -f

#######################################
# the rules
#

# if nothing is given, make the target
all: $(TARGET)

# link the target
$(TARGET): $(OBJS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) 

# the source dependencies
depend: 
	$(DEPEND) $(DEPFLAGS) $(SRCS)

# implicit rule for makin object files
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

# clean
clean:
	$(RM) $(OBJS) $(TARGET)

#############################################################################
# DO NOT DELETE

