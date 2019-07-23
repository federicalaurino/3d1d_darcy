# A GetFEM++ interface for SAMG library with applications
**Politecnico di Milano (ITALY)**

**Author** : Federica Laurino

**Mailto** : federica.laurino@polimi.it

**Date** : July 2019

### Previous projects

**Authors** : Stefano Brambilla, Luca Possenti, Simone Di Gregorio, Giorgio Raimondi, Fannie Gerosa, Daniele Cerroni

**Mailto** : s.brambilla93@gmail.com

**Date** : January 2018

**Github Page** : https://github.com/stefano-brambilla-853558

### How to install and run the program

**THE PACKAGE**

    '3d1d_geomat' : contains the code which solves 3D-1D Darcy in the pressure variable solely using the SAMG solver as described in https://doi.org/10.1007/s13137-019-0115-9 (with D.Cerroni)
    'darcy3D_mixed' : contains the code for solving a simple 3D mixed Darcy with preconditioned GMRES 
    'darcy3d1d_HT': contains the code which solves velocity and pressure problem in a vascular network taking into account the effects of the presence of the RBCs (Possenti, Di Gregorio, Gerosa, Cerroni, Brambilla)
    
    
    In each folder, 
    
- `include/` : General include files

- `lib/`     : Main library (to be generated)

- `src/`     : Example sources

- `config.mk`: Specify the variable GETFEM_PREFIX for GetFEM++ library

- `Makefile` : Instruction to install the project (see INSTALLATION)

- `utilities`: other files for problem set up

**INSTALLATION**

Prerequisites

* You need the open source finite element library "GetFEM++" See http://download.gna.org/getfem/html/homepage Version >= 5.1 is necessary. You have to modify the path to the GetFEM library in `config.mk`:

    GETFEM_PREFIX=/home/.../path/to/.../getfem

  Alternatively, at MOX cluster use the `module.sh` file:
    ```
    $ source configure.sh
    ```

* SAMG LICENCE is not required but recommended for more efficient simulations. To abilitate the usage of SAMG library in the code modify in `configure..sh` the flag

    WITH_SAMG = 1

  and use
  
    ```
    $ source configure.sh
    ```

* Gnuplot: Gnuplot is NOT required, but it can be used to visualize residuals. To use it, uncomment lines within the code and see the GNUPLOT_Istruzioni_installazione to install. https://sourceforge.net/projects/gnuplot/files/gnuplot/

BEWARE: Recall to add the GetFEM library path to LD_LIBRARY_PATH. Example:

    $ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib

======================

### Before installation

You must modify, in `config.mk`, the path to the GetFEM library
``` 
GETFEM_PREFIX=/home/.../path/to/.../getfem
``` 
In `config.mk`, you can also set the flags for optimized or debug installation, and for activating the comments at runtime:
``` 
DEBUG= no
VERBOSE= yes
``` 

Finally, you need to export to LD_LIBRARY_PATH the paths to GetFEM, Boost and Qhull libraries;
this can be done using the modules system (from a MOX computer) or setting manually the paths in the file configure.sh.
Before compiling, call:
``` 
$ source configure.sh
``` 

BEWARE: configure.sh contains the details of SAMG LICENCE: if you do not have a SAMG LICENCE, you should use:
``` 
$ WITH_SAMG = 0
```
======================

### Installation Build in each folder:

In the `include/` folder to build the  (static) library "libproblem3d1d"
```
$ make
``` 
In the `src/` folder to build the examples
```
$ make 
```

BEWARE: 
If you want non-optimized program type:
``` 
$ make DEBUG=yes 
``` 
By default DEBUG=no.

The following macro are defined and exported
``` 
CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include

CXXFLAGS=-std=c++11 

OPTFLAGS=-O3 -DNDEBUG -march=native

LDFLAGS=-L$(GETFEM_PREFIX)/lib

LIBRARIES=-lgetfem
```


DEV ENVIRONMENT

OS : CentOS Linux 7 64-bit

Compiler : g++-5.2.1

GetFEM lib : 5.2

Compiler : g++-5.2.1

GetFEM lib : 5.2

``### Before installation

You must modify, in `config.mk`, the path to the GetFEM library
``` 
GETFEM_PREFIX=/home/.../path/to/.../getfem
``` 
In `config.mk`, you can also set the flags for optimized or debug installation, and for activating the comments at runtime:
``` 
DEBUG= no
VERBOSE= yes
``` 

Finally, you need to export to LD_LIBRARY_PATH the paths to GetFEM, Boost and Qhull libraries;
this can be done using the modules system (from a MOX computer) or setting manually the paths in the file configure.sh.
Before compiling, call:
``` 
$ source configure.sh
``` 

BEWARE: configure.sh contains the details of SAMG LICENCE: if you do not have a SAMG LICENCE, you should use:
``` 
$ source configure_nosamg.sh
``` 

======================

### Installation
Build the whole project with:
``` 
$ make
``` 
It first build the (static) library "libproblem3d1d" by calling
the Makefile in `include/`:
``` 
$ make -C include/
``` 
Then, it calls the inner makefiles provided for all examples.

It is also possible to build a single example, e.g. "2_singlebranch", with:
``` 
$ make library

$ make -C src/2_singlebranch
``` 

BEWARE: 
If you want non-optimized program type:
``` 
$ make DEBUG=yes 
``` 
By defaul DEBUG=no.

The following macro are defined and exported
``` 
CPPFLAGS=-I../../include -I$(GETFEM_PREFIX)/include

CXXFLAGS=-std=c++11 -D=M3D1D_VERBOSE_

OPTFLAGS=-O3 -DNDEBUG -march=native

LDFLAGS=-L$(GETFEM_PREFIX)/lib

LIBRARIES=-lgetfem
``` 
Recall that any macro may be overrlued by specifying it when calling 
make. Example: 
``` 
$ make CXXFLAGS+=-DSOMETHING OPTFLAGS=-g
``` 

======================

### Documentation (only implemented for Notaro's code)
The documentation is produced by doxygen. The file Doxyfile contains 
the common doxygen configuration for all examples.
Build the code documentation with:
``` 
$ make pdf
``` 
It first fills doc/ with code documentation ($ make doc) and then compile
the .tex files to produce a portable file ($ pdflatex doc/latex/refman.tex).
You can visualize the documentation with
``` 
$ evince doc/latex/refman.pdf
``` 

## MAKE OPTIONS
All examples are provided with a Makefile which accepts the following
options:
-  all       : makes the example
-  clean     : as it says
-  distclean : clean and also deletes temporary file and local doc directory
Being "all" the first target of the makefile, to compile the examples is
sufficient to type make. 
In addition the external Makefile (./Makefile) has the following options:
-  doc       : produces the documentation (html, tex)
-  pdf       : produces a guide in portable format
- library    : build the library from files in include/

## RUN EXAMPLES
To run a specific example, go to the related subdirectory
``` 
$ cd src/2_singlebranch
``` 
Build the program
``` 
$ make






