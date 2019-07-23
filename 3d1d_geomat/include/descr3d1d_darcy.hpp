/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d1d Darcy problem - IJGE
======================================================================*/
/*! 
  @file   Adescr3d1d_darcy.cpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Definition of the aux class for algorithm description strings for transport problem.
 */
 

#ifndef M3D1D_DESCR3D1D_HPP_
#define M3D1D_DESCR3D1D_HPP_

#include <string>

namespace getfem {

//! Class to import the descriptors of the coupled 3D/1D solver

struct descr3d1d_darcy {  
    
    //TODO: IMPORT MESH_FILET
   
	//! Absolute path to the vessel mesh file
	std::string MESH_FILEV;
    //! Identifier of vessel mesh tipe
	std::string MESH_TYPET;
	//! Identifier of vessel mesh tipe
	std::string MESH_TYPEV;
	//! Identifier of tissue pressure's FEM type
	std::string FEM_TYPET_P;
	//! Identifier of vessel pressure's FEM type
	std::string FEM_TYPEV_P;
	//! Identifier of tissue coefficients' FEM type
	std::string FEM_TYPET_DATA;
    //! Identifier of vessel coefficients' FEM type
	std::string FEM_TYPEV_DATA;
    //! Identifier of vessel integration method type
	std::string IM_TYPET;
	//! Identifier of vessel integration method type
	std::string IM_TYPEV;
	//! Output directory for transport problem
	std::string OUTPUT;

	// Solver information
	//! Identifief of the monolithic solver for transport problem
	std::string SOLVE_METHOD;
	//! Maximum number of iterations (iterative solvers)
	size_type   MAXITER;
	//! Mamimum residual (iterative solvers)
	scalar_type RES; 
	//! Number of target points for the tissue-to-vessel average
	size_type   NInt;
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Import algorithm specifications from file .param
	void import(ftool::md_param & fname) 
	{
		FILE_ = fname;
		MESH_FILEV  = FILE_.string_value("MESH_FILEV","1D points file for transport problem");
        MESH_TYPET  = FILE_.string_value("MESH_TYPET","3D mesh type");
		MESH_TYPEV  = FILE_.string_value("MESH_TYPEV","1D mesh type");
		FEM_TYPET_P   = FILE_.string_value("FEM_TYPET_P","FEM 3D tissue - pressure");
		FEM_TYPEV_P   = FILE_.string_value("FEM_TYPEV_P","FEM 1D vessel - pressure");
        FEM_TYPET_DATA   = FILE_.string_value("FEM_TYPET_DATA","FEM 3D tissue - data");
		FEM_TYPEV_DATA   = FILE_.string_value("FEM_TYPEV_DATA","FEM 1D vessel - data");
        IM_TYPET 	= FILE_.string_value("IM_TYPET","Name of integration method");
		IM_TYPEV 	= FILE_.string_value("IM_TYPEV","Name of integration method");
		SOLVE_METHOD = FILE_.string_value("SOLVE_METHOD", "Monolithic Solver"); 
		if (SOLVE_METHOD != "SuperLU") { // iterative solver
			MAXITER  = FILE_.int_value("MAXITER", "Max number of sub-iterations");
			RES = FILE_.real_value("RES"); if (RES == 0.) RES = 2.0e-10;
		}
		NInt = size_type(FILE_.int_value("NInt", "Node numbers on the circle for the nonlocal term"));  
		OUTPUT = FILE_.string_value("OUTPUT","Output Directory");
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const descr3d1d_darcy & descr
		)
	{ 
		cout << "---- PROBLEM DESCRIPTORS--------------------------" << endl;
		cout << " MESH FILE 1D problem      : " << descr.MESH_FILEV  << endl;
		cout << " IM  TYPE  3D problem      : " << descr.IM_TYPET    << endl;
		cout << " IM  TYPE  1D problem      : " << descr.IM_TYPEV	<< endl;
		cout << " FEM TYPE  3D pressure     : " << descr.FEM_TYPET_P << endl;
		cout << " FEM TYPE  1D pressure     : " << descr.FEM_TYPEV_P << endl;
		cout << " FEM TYPE  3D data         : " << descr.FEM_TYPET_DATA << endl;
		cout << " FEM TYPE  1D data         : " << descr.FEM_TYPEV_DATA << endl;
        cout << " INTEEGRATION METHOD 3d    : " << descr.IM_TYPET << endl;
        cout << " INTEEGRATION METHOD 1d    : " << descr.IM_TYPEV << endl;
        cout << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
