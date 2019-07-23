# A Getfem interface for SAMG library with applications
**Politecnico di Milano (ITALY)**

**Author** : Federica Laurino

**Mailto** : federica.laurino@polimi.it

**Date** : July 2019

### Previous projects

**Authors** : Stefano Brambulla, Luca Possenti, Simone Di Gregorio, Giorgio Raimondi, Fannie Gerosa, Daniele Cerroni

**Mailto** : s.brambilla93@gmail.com

**Date** : January 2018

**Github Page** : https://github.com/stefano-brambilla-853558

### How to install and run the program
**THE PACKAGE**

    'fluid_ht_curvature' : The folder containing the code which solves velocity and pressure problem (Notaro, Possenti, Di Gregorio code)
    'transport': The folder containing the code which solves transport problem (Brambilla code)

The details of the two libraries are in their respective README.md

**INSTALLATION**

Prerequisites

* You need the open source finite element library "GetFEM++" See http://download.gna.org/getfem/html/homepage Version >= 5.1 is necessary

* SAMG LICENCE is not required but recommended for more efficient simulations.

* Gnuplot: Gnuplot is NOT required, but it can be used to visualize residuals. To use it, uncomment lines within the code and see the GNUPLOT_Istruzioni_installazione to install. https://sourceforge.net/projects/gnuplot/files/gnuplot/

BEWARE: Recall to add the GetFEM library path to LD_LIBRARY_PATH. Example:

$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib

Finally, modify the file config.mk as described in the README.md.

======================

### Installation Build the whole project with:
```
$ make
``` 
It first build the shared library "libproblem3d1d.so" and his examples by calling the Makefile in fluid_ht_curvature/:
```
$ make -C fluid_ht_curvature/
```
Then, it build the shared library "libtransport3d1d.so" and his examples by calling the Makefile in transport/:
```
$ make -C transport/
```
It is also possible to build a single example, by building the two libraries with:
```
$ make library
```
and then calling in the folder example:
```
$ make
```
Alternatively, you can build the two libraries separately, calling, in the rispective folders,
```
$ make
```
Remember to modify the config.mk file, following the instructions of the README.md files.
This can be helpful if a new release of "libproblem3d1d" is available. See: https://github.com/lpossenti/MANworks_ht_curvature for the latest release.

DEV ENVIRONMENT

OS : CentOS Linux 7 64-bit

Compiler : g++-5.2.1

GetFEM lib : 5.2

Compiler : g++-5.2.1

GetFEM lib : 5.2
