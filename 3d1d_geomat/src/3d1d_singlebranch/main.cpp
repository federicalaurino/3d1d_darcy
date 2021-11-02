/*! 
  @file   main.cpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   June 2019.
  @brief  Main program for test simulations.
 */

#include <iostream>
#include <getfem/bgeot_config.h> // for FE_ENABLE_EXCEPT
#include <darcy3d1d.hpp>
 
using namespace getfem;

//! main program
int main(int argc, char *argv[]) 
{

	GMM_SET_EXCEPTION_DEBUG; // Exceptions make a memory fault, to debug.
	FE_ENABLE_EXCEPT;        // Enable floating point exception for Nan.

	try {   
		// Declare a new problem 
		getfem::darcy3d1d p; 
		// Initialize the problem
		p.init(argc, argv);
		// assemble        
		p.assembly();    
		// solve      
		#if (WITH_SAMG == 1)
	 		if (!p.solve_samg()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step 
		#else     
     			if (!p.solve()) GMM_ASSERT1(false, "solve procedure has failed");  // the export is in the solve at each time step
   		#endif
				      
  		p.export_vtk();
    
   		// p.test();
	}

	GMM_STANDARD_CATCH_ERROR;

	return 0; 
	
} /* end of main program */

