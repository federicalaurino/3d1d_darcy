/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d mixed Darcy 
======================================================================*/
/*! 
  @file   darcy3d1d.hpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Declaration of the main class for the 3D/1D coupled Darcy problem.
 */
 
#ifndef M3D1D_DARCY3D1D_HPP_
#define M3D1D_DARCY3D1D_HPP_

#define WITH_AMG

// GetFem++ libraries
#include <getfem/getfem_mesh_fem.h> 
#include <gmm/gmm.h>

// Project headers
#include <defines.hpp>      
#include <node.hpp>
#include <dof3dmixed.hpp>
#include <descr3dmixed.hpp>
#include <param3dmixed.hpp>
#include <darcy_precond.hpp>

#define M3D1D_VERBOSE_
 
namespace getfem {

	//!	Main class defining the 3D Darcy mixed problem
	class darcy3dmixed{ 

		public:
		darcy3dmixed(void) : 
			mimt(mesht), 
			mf_Pt(mesht), mf_Ut(mesht),
			mf_coeft(mesht)
		{} 
		//! Initialize the problem
		void init (int argc, char *argv[]);
		//! Assemble the problem
		void assembly (void);
		//! Solve the problem
		bool solve (void);
		//! Export the solution
		void export_vtk (void); 

		    
		protected:
		 
		//! Mesh for the interstitial tissue @f$\Omega@f$ (3D)
		mesh mesht;
		//! Intergration Method for the interstitial tissue @f$\Omega@f$
		mesh_im mimt;
		//! Finite Element Method for the interstitial pressure @f$p_t@f$
		mesh_fem mf_Pt;
		//! Finite Element Method for the velocities
		mesh_fem mf_Ut;
		//! Finite Element Method for PDE coefficients defined on the interstitial volume
		mesh_fem mf_coeft;  

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
		descr3dmixed descr_darcy;
		//! Physical parameters (dimensionless)
		param3dmixed param_darcy;
		//! Number of degrees of freedom
		dof3dmixed dof_darcy;
			
		//! List of BC nodes of the tissue
		vector< node > BC_darcy;
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
		void build_tissue_boundary(void);
		//! Build the list of vessel boundary (and junctions) data 
		void build_vessel_boundary(void);

		// Aux methods for assembly
		//! Assembling the 3D Darcy block
		void assembly3d_darcy (void);
		//! Assembling the BC
		void assembly_bc (void);

		//! Aux funtions for assembly_bc
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
		void asm_tissue_bc
		    (MAT & M, VEC & F,
		    const mesh_im & mim,
		    const mesh_fem & mf_u,
		    const mesh_fem & mf_data,
		    const std::vector<getfem::node> & BC
		    );
		    
		//! Build the monolithic rhs FM_transp by blocks
		void assembly_rhs(void);
	}; //end of class darcy3dmixed

}  //end of namespace

#endif
