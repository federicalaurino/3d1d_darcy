/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d1d Darcy problem - IJGE
======================================================================*/
/*! 
  @file   darcy3d1d.hpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Declaration of the main class for the 3D/1D coupled Darcy problem.
 */
 
#ifndef M3D1D_DARCY3D1D_HPP_
#define M3D1D_DARCY3D1D_HPP_
#define M3D1D_VERBOSE_

// GetFem++ libraries
#include <getfem/getfem_assembling.h> 
//#include <getfem/getfem_import.h>
#include <getfem/getfem_export.h>   
#include <getfem/getfem_regular_meshes.h>
//#include <getfem/getfem_mesh.h>
//#include <getfem/getfem_derivatives.h>
//#include <getfem/getfem_superlu.h>
//#include <getfem/getfem_partial_mesh_fem.h>
#include <getfem/getfem_interpolated_fem.h>
#include <gmm/gmm.h>
//#include <gmm/gmm_matrix.h>
#include <gmm/gmm_inoutput.h>
#include <gmm/gmm_iter_solvers.h>
// Standard libraries
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
// Project headers
#include <defines.hpp>      
#include <mesh1d_darcy.hpp>
#include <utilities_darcy.hpp>
#include <assembling3d1d_darcy.hpp>
#include <node.hpp>
#include <dof3d1d_darcy.hpp>
#include <descr3d1d_darcy.hpp>
#include <param3d1d_darcy.hpp>
#include <AMG_Interface.hpp>

#include <cmath>

//TODO not tissue and vessel but reservoir and well

#define M3D1D_VERBOSE_
 
 namespace getfem {

//!	Main class defining the coupled 3D/1D Darcy problem.
class darcy3d1d{ 

public:
	darcy3d1d(void) : 
		mimt(mesht),  mimv(meshv),
		mf_Pt(mesht), mf_Pv(meshv),
		mf_coeft(mesht),mf_coefv(meshv)
	{} 
	//! Initialize the problem
	void init (int argc, char *argv[]);
	//! Assemble the problem
	void assembly (void);
	//! Solve the problem
	bool solve (void);
	bool solve_samg (void);
	//! Export the solution
 	void export_vtk (void); 
	    
	/*//! Compute mean concentration in the tissue
	inline scalar_type mean_ct (void){ 
		return asm_mean(mf_Ct, mimt, 
			gmm::sub_vector(UM_transp, gmm::sub_interval(0,  dof_transp.Ct()))); 
	}
	void test(void);
	void test2(void);*/

protected:
	 
	//! Mesh for the interstitial tissue @f$\Omega@f$ (3D)
	mesh mesht;
	//! Mesh for the vessel network @f$\Lambda@f$ (1D)
	mesh meshv;
	//! Intergration Method for the interstitial tissue @f$\Omega@f$
	mesh_im mimt;
	//! Intergration Method for the vessel network @f$\Lambda@f$
	mesh_im mimv;
	//! Finite Element Method for the interstitial pressure @f$p_t@f$
	mesh_fem mf_Pt;
    //! Finite Element Method for the vessel pressure @f$p_v@f$
	mesh_fem mf_Pv; 
	//! Finite Element Method for PDE coefficients defined on the interstitial volume
	mesh_fem mf_coeft;  
	//! Finite Element Method for PDE coefficients defined on the network
	mesh_fem mf_coefv;
	
    //! Dimension of the tissue domain (3)
	size_type DIMT;
	//! Number of vertices per branch in the vessel network
	vector_size_type nb_vertices;
	//! Number of branches in the vessel network
	size_type nb_branches;
	//! Number of extrema of the vessel network
	size_type nb_extrema;
	//! Number of junctions of the vessel network
	size_type nb_junctions;
    
    //! Input file
	ftool::md_param PARAM;
	//! Algorithm description strings (mesh files, FEM types, solver info, ...) 
	descr3d1d_darcy descr_darcy;
	//! Physical parameters (dimensionless)
	param3d1d_darcy param_darcy;
	//! Number of degrees of freedom
	dof3d1d_darcy dof_darcy;
		
	
	//! List of BC nodes of the network
	vector< node > BCv_darcy;	
	//! List of BC nodes of the tissue
	vector< node > BCt_darcy;
    //! List of junction nodes of the network
	vector< node > Jv_darcy;
		
	//! Monolithic matrix for the coupled problem
	sparse_matrix_type AM_darcy;
	//! Monolithic array of unknowns for the coupled problem
	vector_type UM_darcy;
	//! Monolithic right hand side for the coupled problem
	vector_type FM_darcy;
    
	// Aux methods for init
	//! Import algorithm specifications
	void import_data(void);
	//! Import mesh for tissue (3D) and vessel (1D)  
	void build_mesh(void); 
	//! Set finite elements methods and integration methods 
	void set_im_and_fem(void);
	//! Build problem parameters
	void build_param(void);
	//! Build the list of tissue boundary data 
	/*!	Face numbering:
		  0 : {x = 0 }  "back"
		  1 : {x = Lx}  "front"
		  2 : {y = 0 }  "left"
		  3 : {y = Ly}  "right"
		  4 : {z = 0 }  "bottom"
		  5 : {z = Lz}  "top"
	 */
	void build_tissue_boundary(void);
	//! Build the list of vessel boundary (and junctions) data 
	void build_vessel_boundary(void);
    //! Build the monolithic matrix AM_transp by blocks
	
    void assembly3d_darcy (void);
    void assembly1d_darcy (void);
    void assembly3d1d_darcy (void);
    void assembly_bc (void);
    
    template<typename MAT, typename VEC>
    void asm_tissue_bc
        (VEC & F,
        MAT & M,
        const mesh_im & mim,
        const mesh_fem & mf_p,
        const mesh_fem & mf_data,
        const std::vector<getfem::node> & BC,
        const scalar_type beta
        );
    
    template<typename MAT, typename VEC>
    void asm_network_bc
        (VEC & F, MAT & M, 
        const mesh_im & mim,
        const mesh_fem & mf_c,
        const mesh_fem & mf_data,
        const std::vector<getfem::node> & BC,
        const scalar_type beta,
        const VEC & R);
	//! Build the monolithic rhs FM_transp by blocks
	//void assembly_rhs(void);
	//void update(vector_type);
    
	
	
}; //end of class darcy3d1d

}  //end of namespace

#endif
