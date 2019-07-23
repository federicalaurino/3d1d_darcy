/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
3d1d Darcy problem - IJGE
======================================================================*/
/*! 
  @file   utilities.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_HPP_
#define M3D1D_UTILITIES_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh_fem.h>
#include <gmm/gmm.h>

namespace getfem {

	//! Aux function to compute the diameter of an element
	scalar_type 
	estimate_h(const mesh & mesh, const size_type i); 

	//! Aux function to read an array of string split by a delim. 
	//! Store the results in a pre-constructed vector
	std::vector <std::string> &
	split(const std::string & s, 
		  char delim, 
		  std::vector<std::string> & elems
		  ) ;

	//! Aux function to read an array of string split by a delim. 
	//! Return a new vector
	std::vector <std::string>
	split(const std::string & s, 
		  char delim
		  ) ;

	//! Build the integral of FE base functions
	//! @f$ \Phi = \int_{\Omega} \phi(x)~dx @f$
	/*!
	 @param V    The array of computed values
	 @param mim  The integration method to be used
	 @param mf   The finite element method
	 @param rg   The region where to integrate
	 */ 
	template
	<typename VEC>
	void 
	asm_basis_function(VEC & V,
		 			   const mesh_im & mim,
		 			   const mesh_fem & mf,
		 			   const mesh_region & rg = mesh_region::all_convexes())		
	{
		GMM_ASSERT1(mf.get_qdim() == 1, 
			"invalid data mesh fem for pressure (Qdim=1 required)");
		generic_assembly 
		assem("V$1(#1)+=comp(Base(#1));");
		assem.push_mi(mim);
		assem.push_mf(mf);
		assem.push_vec(V);
		assem.assembly(rg);
	}

	//! Aux function to extract the radius of the ith branch, R[i] 
	template
	<typename VEC>
	scalar_type
	compute_radius(const mesh_im & mim,
		 		   const mesh_fem & mf_coef,
		 		   const VEC & R,
		 		   const size_type & rg) 
	{
		VEC rbasis(mf_coef.nb_dof());
		asm_basis_function(rbasis, mim, mf_coef, rg);
		size_type firstcv = mf_coef.linked_mesh().region(rg).index().first_true();
		size_type rg_size = mf_coef.linked_mesh().region(rg).size();
		gmm::scale(rbasis, 1.0/estimate_h(mf_coef.linked_mesh(), firstcv));
		return gmm::vect_sp(rbasis, R)/rg_size;
	}

} /* end of namespace */

#endif
