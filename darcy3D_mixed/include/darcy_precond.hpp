/*  Class for darcy preconditioner for the 3D darcy problem:
    P = [diag^-1 0 ; 0 schur^-1]
    where diag is the diagonal of the velocity mass matrix and 
    schur is an approximation of the schur complement given by the jumps of the pressure
    across the internal faces + a boundary term
*/

#ifndef M3D1D_DARCY_PRECOND_HPP_
#define M3D1D_DARCY_PRECOND_HPP_

#include <defines.hpp>
#include <gmm/gmm.h>
#include <gmm/gmm_precond.h>
#include <getfem/getfem_generic_assembly.h>
#include <chrono>
#if WITH_SAMG == 1
    #include <AMG_Interface.hpp>
#endif   

  
template <typename Matrix> struct darcy_precond{
  
    gmm::diagonal_precond <Matrix> D;
    Matrix schur; 

    getfem::size_type dof_u;
    getfem::size_type dof_p;
    
    getfem::size_type nrows = dof_u + dof_p;
    getfem::size_type ncols = dof_u + dof_p;
    
    // Case of Dirichlet boundary condition
    void build_with (const Matrix &A, const getfem::mesh_fem & mf_u, const getfem::mesh_fem & mf_p, const getfem::mesh_im & mim){

        dof_u = mf_u.nb_dof();
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
        
        D.build_with(M);
        
        //approximated schur block
 
        gmm::resize(schur, dof_p, dof_p); gmm::clear(schur);
        
        const getfem::mesh &mesh = mf_p.linked_mesh();
        getfem::mesh_region inner_faces = getfem::inner_faces_of_mesh(mesh);
        getfem::mesh_region outer_faces;
        getfem::outer_faces_of_mesh(mesh, outer_faces);
        
        getfem::ga_workspace wp;
        std::vector<double> p(mf_p.nb_dof());
        wp.add_fem_variable("p", mf_p, gmm::sub_interval(0,mf_p.nb_dof()), p);  

        wp.add_expression( "1 / element_size * (p - Interpolate(p, neighbour_elt))"
                            " * (Test_p - Interpolate(Test_p, neighbour_elt))",
                            mim, inner_faces);
        // to remove in case of mix/neumann condition
        wp.add_expression("1/element_size*p*Test_p",
                            mim, outer_faces);
        //#endif
        wp.assembly(2);

        gmm::copy(wp.assembled_matrix(), schur);
        
    } // end of build with
    
    
    // default constructor
    darcy_precond(void) {} 
    //constructor
    darcy_precond(const Matrix &A, 
                  const getfem::mesh_fem & mf_u, 
                  const getfem::mesh_fem & mf_p, 
                  const getfem::mesh_im & mim) {
        dof_u = mf_u.nb_dof();
        dof_p = mf_p.nb_dof();
        build_with(A, mf_u, mf_p, mim); 
    } 
    
}; // end of struct
    
    
// Multiplication of the prec by a vector: P*vec -> res.
// In our case we build P but what we have to do is P^-1*vec -> we solve the linear system P res =  vec 

template < typename Matrix, typename V1, typename V2>
void mult(const darcy_precond<Matrix> &P, const V1 &vec, V2 &res) {
    
    // multiplication of the diag term -> use the mult method defined for diagonal_precond
    gmm::mult(P.D,
              gmm::sub_vector(vec, gmm::sub_interval(0, P.dof_u)), 
              gmm::sub_vector(res, gmm::sub_interval(0, P.dof_u)));
    
    // multiplication of the Schur term
    #if (WITH_SAMG==1)
        // AMG solver
        gmm::csr_matrix <scalar_type> schur_csr;
        gmm::copy(P.schur, schur_csr);
        std::vector <double> rhs; gmm::resize(rhs, P.dof_p);
        gmm::copy(gmm::sub_vector(vec, gmm::sub_interval(P.dof_u, P.dof_p)), rhs);
        std::cout << "------------------Solving with SAMG" << std::endl;
        AMG sys(schur_csr, rhs,1);
        sys.csr2samg();
        sys.solve();
        //solution
        for(int i=0; i< P.dof_p; i++)
            res[P.dof_u + i] = sys.getsol()[i];
    #else 
        double cond;
        // Direct solver
        gmm::SuperLU_solve(P.schur, 
                       gmm::sub_vector(res, gmm::sub_interval(P.dof_u, P.dof_p)), 
                       gmm::sub_vector(vec, gmm::sub_interval(P.dof_u, P.dof_p)), 
                       cond);
    #endif
        
} // end mult

//Redefinition of some gmm functions for darcy_precond

namespace gmm{
    
struct owned_implementation {};

template <>
    struct principal_orientation_type<owned_implementation> {
        using potype = owned_implementation;
    }; 

template <typename MATRIX>
 struct linalg_traits<::darcy_precond<MATRIX>> {
        using this_type = ::darcy_precond<MATRIX>;
        using sub_orientation = owned_implementation;

        static gmm::size_type nrows(const this_type &m) { return m.nrows; }
        static gmm::size_type ncols(const this_type &m) { return m.ncols; }
    };   
    
template <typename L1, typename L2, typename L3>
    inline
    void mult_spec(const L1 &m, const L2 &src, L3 &dst, owned_implementation tag)
    {
       mult(m, src, dst);
    }

}    
 // end of namespace gmm

#endif