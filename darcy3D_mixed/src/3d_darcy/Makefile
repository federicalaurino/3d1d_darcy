include ../../../config.mk

CPPFLAGS=-I../../include -I$(mkGetfemInc) -I$(mkBoostInc) 
CXXFLAGS+=-std=c++14
ifeq ($(WITH_SAMG),1)
  CXXFLAGS += -I ${SAMG}/
  CXXFLAGS+= -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS
endif
# -D=M3D1D_VERBOSE_

ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif
LDFLAGS+=-L../../lib -L$(mkGetfemLib) -L$(mkBoostLib) -Wl,-rpath=../../lib
LIBRARIES+=-lproblem3d1d -lgetfem -lutil -lboost_iostreams -lboost_system -lboost_filesystem
ifeq ($(WITH_SAMG),1)
  LDFLAGS += -L/opt/lib/samg/
  LIBRARIES += -lamg -liomp5
endif
LINKLIBS := ../../lib/libproblem3d1d.a



SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)
EXEC=M3D1D

OUTDIR=vtk

.PHONY: all clean distclean

all: library $(EXEC) 
	@echo
	@echo Compilation completed!

%.o: %.cpp $(LINKLIBS) 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<

$(EXEC): $(OBJS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBRARIES)

library:
	$(MAKE) -C ../../include

clean:
	$(RM) $(OBJS) $(EXEC) *~ *.log

distclean: clean
	$(RM) *.txt $(OUTDIR)/*
