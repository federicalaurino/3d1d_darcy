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

**INSTALLATION**

Prerequisites

* You need the open source finite element library "GetFEM++" See http://download.gna.org/getfem/html/homepage Version >= 5.1 is necessary

* SAMG LICENCE is not required but recommended for more efficient simulations.

* Gnuplot: Gnuplot is NOT required, but it can be used to visualize residuals. To use it, uncomment lines within the code and see the GNUPLOT_Istruzioni_installazione to install. https://sourceforge.net/projects/gnuplot/files/gnuplot/

BEWARE: Recall to add the GetFEM library path to LD_LIBRARY_PATH. Example:

$ export LD_LIBRARY_PATH=/home/...path/to.../getfem/lib

Finally, modify the file config.mk as described in the README.md.

======================

### Installation Build in each folder:

In the **include** folder to build the library
```
$ make
``` 
In the **src** folder to build the examples
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
