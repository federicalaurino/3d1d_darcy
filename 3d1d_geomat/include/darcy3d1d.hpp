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

// Getfem++ libraries   
#include <getfem/getfem_mesh_fem.h>
#include <gmm/gmm.h>

// Project headers
#include <defines.hpp>      
#include <mesh1d_darcy.hpp>
#include <utilities_darcy.hpp>
#include <assembling3d1d_darcy.hpp>
#include <node.hpp>
#include <dof3d1d_darcy.hpp>
#include <descr3d1d_darcy.hpp>
#include <param3d1d_darcy.hpp>

 
namespace getfem {

	//!	Main class defining the coupled 3D/1D Darcy problem.
	class darcy3d1d{ 

		public:
			darcy3d1d () : 
				mimt(mesht),  mimv(meshv),
				mf_Pt(mesht), mf_Pv(meshv),
				mf_coeft(mesht),mf_coefv(meshv)
			{} 
			//! Initialize the problem
			void init (int argc, char *argv[]);
			//! Assemble the problem
			void assembly ();
			//! Solve the problem
			bool solve ();

		#if (WITH_SAMG == 1)
			bool solve_samg ();
		#endif
			//! Export the solution
		 	void export_vtk (); 
		    //! Test with a given source on the 1D manifold
		    void test ();
		
		protected:
		 
			//! Mesh for the 3D outer domain @f$\Omega@f$ (3D)
			mesh mesht;
			//! Mesh for the 1D inner domain @f$\Lambda@f$ (1D)
			mesh meshv;
			//! Intergration Method for the outer domain @f$\Omega@f$
			mesh_im mimt;
			//! Intergration Method for inner domain @f$\Lambda@f$
			mesh_im mimv;
			//! Finite Element Method for the external pressure @f$p_t@f$
			mesh_fem mf_Pt;
		    //! Finite Element Method for the internal pressure @f$p_v@f$
			mesh_fem mf_Pv; 
			//! Finite Element Method for PDE coefficients defined on the outer domain
			mesh_fem mf_coeft;  
			//! Finite Element Method for PDE coefficients defined on the inner domain
			mesh_fem mf_coefv;
			
		    //! Dimension of the tissue domain (3)
			size_type DIMT;
			//! Number of vertices per branch in the inner domain
			vector_size_type nb_vertices;
			//! Number of branches in the inner domain
			size_type nb_branches;
			//! Number of extrema of the inner domain
			size_type nb_extrema;
			//! Number of junctions of the inner domain
			size_type nb_junctions;
		    
		    //! Input file
			ftool::md_param PARAM;
			//! Algorithm description strings (mesh files, FEM types, solver info, ...) 
			descr3d1d_darcy descr_darcy;
			//! Physical parameters (dimensionless)
			param3d1d_darcy param_darcy;
			//! Number of degrees of freedom
			dof3d1d_darcy dof_darcy;
				
			
			//! List of BC nodes of the inner domain
			std::vector <node> BCv_darcy;	
			//! List of BC nodes of the outer domain
			std::vector <node> BCt_darcy;
		    //! List of junction nodes of inner domain
			std::vector <node> Jv_darcy;
				
			//! Monolithic matrix for the coupled problem
			sparse_matrix_type AM_darcy;
			//! Monolithic array of unknowns for the coupled problem
			vector_type UM_darcy;
			//! Monolithic right hand side for the coupled problem
			vector_type FM_darcy;
		    
			// Aux methods for init
			//! Import algorithm specifications
			void import_data ();
			//! Import mesh for outer (3D) and inner (1D) domains  
			void build_mesh (); 
			//! Set finite elements methods and integration methods 
			void set_im_and_fem ();
			//! Build problem parameters
			void build_param ();
			//! Build the list of 3D boundary data 
			void build_3d_boundary ();
			//! Build the list of 1D boundary (and junctions) data 
			void build_1d_boundary ();
		    
		    //! Aux methods for assembly 
			//! Assembling the 3D Darcy block
		    void assembly3d_darcy ();
		    //! Assembling the 1D Darcy block
		    void assembly1d_darcy ();
		    //! Assembling the coupling terms
		    void assembly3d1d_darcy ();
		    //! Assembling the BC
		    void assembly_bc ();
		    
		    //! Aux funtions for assembly_bc
		    template<typename MAT, typename VEC>
		    void asm_bc_3d
		        (VEC & F,
		        MAT & M,
		        const mesh_im & mim,
		        const mesh_fem & mf_p,
		        const mesh_fem & mf_data,
		        const std::vector<getfem::node> & BC,
		        const scalar_type beta
		        );
		    template<typename MAT, typename VEC>
		    void asm_bc_1d
		        (VEC & F, MAT & M, 
		        const mesh_im & mim,
		        const mesh_fem & mf_c,
		        const mesh_fem & mf_data,
		        const std::vector<getfem::node> & BC,
		        const scalar_type beta,
		        const VEC & R);

	}; //end of class darcy3d1d

}  //end of namespace

#endif
