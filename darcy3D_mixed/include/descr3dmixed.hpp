/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d mixed Darcy 
======================================================================*/
/*! 
  @file   descr3dmixed.hpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Definition of the aux class for algorithm description strings for transport problem.
 */
 
/** @defgroup input User-defined parameters  */

#ifndef M3D1D_DESCR3D1D_TRANSP_HPP_
#define M3D1D_DESCR3D1D_TRANSP_HPP_

#include <string>

namespace getfem {

//! Class to import the descriptors of the 3D problem
/*!
	\ingroup input
 */
struct descr3dmixed {  
    
    //TODO: IMPORT MESH_FILET
   
    //! Identifier of vessel mesh tipe
	std::string MESH_TYPET;
    //! Identifier of tissue velocity's FEM type
	std::string FEM_TYPET_U;
	//! Identifier of tissue pressure's FEM type
	std::string FEM_TYPET_P;
	//! Identifier of tissue coefficients' FEM type
	std::string FEM_TYPET_DATA;
    //! Identifier of vessel integration method type
	std::string IM_TYPET;
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
	    MESH_TYPET  = FILE_.string_value("MESH_TYPET","3D mesh type");
        FEM_TYPET_U   = FILE_.string_value("FEM_TYPET_U","FEM 3D tissue - velocity");
		FEM_TYPET_P   = FILE_.string_value("FEM_TYPET_P","FEM 3D tissue - pressure");
        FEM_TYPET_DATA   = FILE_.string_value("FEM_TYPET_DATA","FEM 3D tissue - data");
        IM_TYPET 	= FILE_.string_value("IM_TYPET","Name of integration method");
		SOLVE_METHOD = FILE_.string_value("SOLVE_METHOD", "Monolithic Solver"); 
		if (SOLVE_METHOD != "SuperLU") { // iterative solver
			MAXITER  = FILE_.int_value("MAXITER", "Max number of sub-iterations");
			RES = FILE_.real_value("RES"); if (RES == 0.) RES = 2.0e-10;
		}
		NInt = size_type(FILE_.int_value("NInt", "Node numbers on the circle for the nonlocal term"));  
		OUTPUT = FILE_.string_value("OUTPUT","Output Directory");
		//STATIONARY = size_type(FILE_.int_value("STATIONARY", "Flag to make the problem stationary")); 
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const descr3dmixed & descr
		)
	{ 
		cout << "---- PROBLEM DESCRIPTORS--------------------------" << endl;
		cout << " IM  TYPE  3D problem      : " << descr.IM_TYPET    << endl;
        cout << " FEM TYPE  3D velocity     : " << descr.FEM_TYPET_U << endl;
		cout << " FEM TYPE  3D pressure     : " << descr.FEM_TYPET_P << endl;
		cout << " FEM TYPE  3D data         : " << descr.FEM_TYPET_DATA << endl;
        cout << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
