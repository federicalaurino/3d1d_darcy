/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d mixed Darcy 
======================================================================*/
/*! 
  @file   darcy3d1d.hpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Definition of the aux class for physical parameters.
 */
 
#ifndef M3D1D_PARAM3D1D_TRANSP_HPP_
#define M3D1D_PARAM3D1D_TRANSP_HPP_

#include <utilities_darcy.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
struct param3dmixed {

	// TODO: Dimensional physical parameters (microcirc applications)
    // TODO Dimensionless physical parameters (test-cases)
    

    // Tissue permeability TODO Generalize to tensor
    vector_type Kp_;
    
    getfem::mesh_fem mf_coeft_;

	// Utils
	//! File .param
	ftool::md_param FILE_;
    
	// Methods
	void build ( ftool::md_param & fname, getfem::mesh_fem mf_coeft ) 
	{
		FILE_ = fname;
        mf_coeft_ = mf_coeft;
        size_type dof_coeft = mf_coeft_.nb_dof();
        
        scalar_type Kp_const = FILE_.real_value("Kp", "Tissue permeability");
        Kp_.assign(dof_coeft, Kp_const);
    }
	
	//! Get Kp
	vector_type & Kp (void) { return Kp_; }
	

}; /* end struct */

} /* end namespace */

#endif
