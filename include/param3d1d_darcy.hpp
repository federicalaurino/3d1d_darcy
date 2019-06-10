/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   param3d1d_transp.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for physical parameters.
 */
 
#ifndef M3D1D_PARAM3D1D_TRANSP_HPP_
#define M3D1D_PARAM3D1D_TRANSP_HPP_

#include <mesh1d_darcy.hpp>    // import_network_radius
#include <utilities_darcy.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
struct param3d1d_darcy {

	// TODO: Dimensional physical parameters (microcirc applications)
    // TODO Dimensionless physical parameters (test-cases)
    
    // Radius for each node of the 1d mesh1d_darcy
    vector_type R_;
    // Tissue permeability TODO Generalize to tensor
    vector_type Kp_;
    // Vessel permeability TODO Generalize to tensor
    vector_type Kv_;
    // interface permeability 
    vector_type kappa_;
    
    getfem::mesh_fem mf_coeft_;
    getfem::mesh_fem mf_coefv_;

	// Utils
	//! File .param
	ftool::md_param FILE_;
    
	// Methods
	void build ( ftool::md_param & fname, getfem::mesh_fem mf_coeft, getfem::mesh_fem mf_coefv ) 
	{
		FILE_ = fname;
        mf_coeft_ = mf_coeft;
		mf_coefv_ = mf_coefv;
        size_type dof_coeft = mf_coeft_.nb_dof();
        size_type dof_coefv = mf_coefv_.nb_dof();
        
        // TODO adimensionalization
        // TODO add the possibility to import R from file 

        scalar_type Rav_ = FILE_.real_value("RADIUS", "Vessel average radius");
        R_.assign(dof_coefv, Rav_);
        
        scalar_type Kp_const = FILE_.real_value("Kp", "Tissue permeability");
        Kp_.assign(dof_coeft, Kp_const);
        scalar_type Kv_const = FILE_.real_value("Kv", "Vessel permeability");
        Kv_.assign(dof_coefv, Kv_const);
        scalar_type kappa_const = FILE_.real_value("kappa", "Interface permeability");
		kappa_.assign(dof_coefv, kappa_const);
	}
	
	//! Get the radius at a given dof
	inline scalar_type R  (size_type i) { return R_[i];  } const
	
	//! Get the radius at a given mesh_region
	scalar_type R  (const getfem::mesh_im & mim, const size_type rg) 
    { 
		return compute_radius(mim, mf_coefv_, R_, rg);  
	}
	//! Get the radius
	vector_type & R (void) { return R_; }
	//! Get Kp
	vector_type & Kp (void) { return Kp_; }
	//! Get Kv
	vector_type & Kv (void) { return Kv_; }
	//! Get kappa
	vector_type & kappa (void) { return kappa_; }
	

}; /* end struct */

} /* end namespace */

#endif
