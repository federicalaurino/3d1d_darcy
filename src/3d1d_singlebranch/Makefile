CPPFLAGS=-I../../3d1d_geomat/include -I$(mkGetfemInc)
CXXFLAGS=-std=c++14 -D=M3D1D_VERBOSE_

#CXXFLAGS += -I ${SAMG}/
#CXXFLAGS+= -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS

ifeq ($(DEBUG),yes)
  OPTFLAGS=-g -Wall
else
  OPTFLAGS=-O3 -march=native
  CPPFLAGS+=-DNDEBUG
endif
LDFLAGS=-L../../3d1d_geomat/lib -L$(mkGetfemLib)

#LDFLAGS += -L/opt/lib/samg/

LIBRARIES=-lproblem3d1d -lgetfem -lblas


#LIBRARIES += -lamg -liomp5

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
	$(MAKE) -C ../../3d1d_geomat/include

clean:
	$(RM) $(OBJS) $(EXEC) *~ *.log

distclean: clean
	$(RM) *.txt $(OUTDIR)/*
