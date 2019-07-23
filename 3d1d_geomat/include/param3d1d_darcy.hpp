/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d1d Darcy problem - IJGE
======================================================================*/
/*! 
  @file   AMG_Interface.cpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Definition of the aux class for physical parameters.
 */
 
#ifndef M3D1D_PARAM3D1D_DARCY_HPP_
#define M3D1D_PARAM3D1D_DARCY_HPP_

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
	        
	        // Flag to import dimensionless param
	        bool IMPORT_DIMLESS = FILE_.int_value("IMPORT_DIMLESS");
	        bool IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS");

	        if (IMPORT_RADIUS) {   /* case R' = const */
	            std::string RFILE = FILE_.string_value("RFILE"); 
	            cout << "  Importing radius values from file " << RFILE << " ..." << endl;
	            std::ifstream ist(RFILE);
	            if (!ist) cerr << "impossible to read from file " << RFILE << endl;
	            import_network_radius(R_, ist, mf_coefv_);
	        }
	        else{
	            scalar_type Rav_ = FILE_.real_value("RADIUS", "Vessel average radius");
	            if (IMPORT_DIMLESS)
	                R_.assign(dof_coefv, Rav_);
	            else{
	                scalar_type L = FILE_.real_value("L", "Characteristic length");
	                R_.assign(dof_coefv, Rav_/L);  
	                } 
	          }

	        scalar_type Kp_const = FILE_.real_value("Kp", "Tissue permeability");
	        scalar_type Kv_const = FILE_.real_value("Kv", "Vessel permeability");
	        scalar_type kappa_const = FILE_.real_value("kappa", "Interface permeability");

	        if (IMPORT_DIMLESS){
	            Kp_.assign(dof_coeft, Kp_const);
	            Kv_.assign(dof_coefv, Kv_const);          
	    		kappa_.assign(dof_coefv, kappa_const);
	        }
	        else{
	            scalar_type L = FILE_.real_value("L", "Characteristic length");
	            scalar_type K = FILE_.real_value("K", "Reference permeability");
	            Kp_.assign(dof_coeft, Kp_const/K);
	            Kv_.assign(dof_coefv, Kv_const/K);          
	            kappa_.assign(dof_coefv, kappa_const * L/K);
	        }
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
