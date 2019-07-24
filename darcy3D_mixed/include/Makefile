include ../../config.mk

CPPFLAGS=-I. -I$(GETFEM_PREFIX)/include
CXXFLAGS+=-std=c++14
OPTFLAGS=-O3 -DNDEBUG -march=native

LIBDIR=../lib
LIBNAME=problem3d1d
LIBFILE=lib$(LIBNAME).a

LIBSRC=$(wildcard *.cpp)
LIBOBJS=$(LIBSRC:.cpp=.o)
LIBHEADERS=$(wildcard *.hpp)

.PHONY: all clean distclean library

all: library
	@echo
	@echo Library installed!

library: $(LIBOBJS)
	install -d $(LIBDIR)
	ar -r $(LIBDIR)/$(LIBFILE) $(LIBOBJS) 

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OPTFLAGS) -o $@ -c $<

clean:
	$(RM) $(LIBOBJS) *~

distclean: clean
	$(RM) $(LIBDIR)/$(LIBFILE)

