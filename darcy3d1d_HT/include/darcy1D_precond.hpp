/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d1d - HT
======================================================================*/
/*! 
  @file   darcy1D_precond.hpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Class for darcy preconditioner for the 3D darcy problem: 
 */

/*  Class for darcy preconditioner for the 3D darcy problem:
    P = [diag^-1 0 ; 0 schur^-1]
    where diag is the diagonal of the velocity mass matrix and 
    schur is an approximation of the schur complement given by the jumps of the pressure
    across the internal faces + a boundary term
*/

#ifndef M3D1D_DARCY1D_PRECOND_HPP_
#define M3D1D_DARCY1D_PRECOND_HPP_

#include <gmm/gmm.h>
#include <gmm/gmm_precond.h>
//#include <getfem/getfem_mesh.h>
#include <getfem/getfem_generic_assembly.h>
#include <defines.hpp>

    
template <typename Matrix> struct darcy1D_precond{
  
    gmm::diagonal_precond <Matrix> D;
    Matrix schur; 
    
    getfem::size_type dof_u;
    getfem::size_type dof_p;
    
    getfem::size_type nrows = dof_u + dof_p;
    getfem::size_type ncols = dof_u + dof_p;
    
   // Case of Dirichlet boundary condition
    void build_with ( const Matrix &A, const getfem::vector<getfem::mesh_fem> & mf_ui, const getfem::mesh_fem & mf_p, const getfem::mesh_im & mim){
//        std::cout << "---------Entering constructor of Precv" << std::endl;
        dof_u = 0;
        for( gmm::size_type i=0; i< mf_ui.size(); i++) 
            dof_u += mf_ui[i].nb_dof();
        dof_p = mf_p.nb_dof();
        
        //build the preconditioner  blocks from the matrix A
        
        //diagonal block

        //Extract the mass matrix
        Matrix M;
        M.resize(dof_u, dof_u);
        gmm::copy (gmm::sub_matrix(A,
            gmm::sub_interval (0, dof_u),
                                   gmm::sub_interval (0, dof_u)
            ), M);

//        std::cout << "---------Build diagonal" << std::endl;
        D.build_with(M);
//        std::cout << "---------End building diagonal" << std::endl;

        
        //approximated schur block
 
        gmm::resize(schur, dof_p, dof_p); gmm::clear(schur);
            

        const getfem::mesh &mesh = mf_p.linked_mesh();
        getfem::mesh_region outer_faces;
        getfem::outer_faces_of_mesh(mesh, outer_faces);

        
        getfem::ga_workspace wp;
        std::vector<double> p(mf_p.nb_dof());
        wp.add_fem_variable("p", mf_p, gmm::sub_interval(0,mf_p.nb_dof()), p);  
        wp.add_expression( "Grad_p.Grad_Test_p",mim);
        // to remove in case of mix/neumann condition
        wp.add_expression("1*p * Test_p  ",
                            mim, outer_faces);
        wp.assembly(2);
//        std::cout << "---------end workspace" << std::endl;
        gmm::copy(wp.assembled_matrix(), schur);           
//        std::cout << "---------end copy schur" << std::endl;
        
    } // end of build with
    
    // default constructor
    darcy1D_precond(void) {}  
    
    //constructor
    darcy1D_precond(const Matrix &A, const getfem::vector<getfem::mesh_fem> & mf_ui, const getfem::mesh_fem & mf_p, const getfem::mesh_im & mim) {
        build_with(A, mf_ui, mf_p, mim); 
    } 
    
    }; // end of struct
    
    
// Multiplication of the prec by a vector: P*vec -> res.
// In our case we build P but what we have to do is P^-1*vec, therefore we solve the linear system P res =  vec 

template < typename Matrix, typename V1, typename V2>
void mult(const darcy1D_precond<Matrix> &P, const V1 &vec, V2 &res) {
    
    // multiplication of the diag term -> use the mult method defined for diagonal_precond
    gmm::mult(P.D,
              gmm::sub_vector(vec, gmm::sub_interval(0, P.dof_u)), 
              gmm::sub_vector(res, gmm::sub_interval(0, P.dof_u)));
    
    // multiplication of the Schur term
    double cond;
    gmm::SuperLU_solve(P.schur, 
                       gmm::sub_vector(res, gmm::sub_interval(P.dof_u, P.dof_p)), 
                       gmm::sub_vector(vec, gmm::sub_interval(P.dof_u, P.dof_p)), 
                       cond);
} // end mult

//Redefinition of some gmm functions for darcy_precond

namespace gmm{
    
template <typename MATRIX>
 struct linalg_traits<::darcy1D_precond<MATRIX>> {
        using this_type = ::darcy1D_precond<MATRIX>;
        using sub_orientation = owned_implementation;

        static gmm::size_type nrows(const this_type &m) { return m.nrows; }
        static gmm::size_type ncols(const this_type &m) { return m.ncols; }
    };   
    
}// end of namespace gmm

#endif
