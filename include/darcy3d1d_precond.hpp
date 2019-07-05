/*  Class for darcy preconditioner for the 3D1D darcy problem:
    P = [diag_t^-1 0  0 0; 0 schur_t^-1; 0 0 diag_v^-1 schur_v^-1]
    where diag is the diagonal of the velocity mass matrix and 
    schur is an approximation of the schur complement given by the jumps of the pressure
    across the internal faces + a boundary term
*/

#include <darcy_precond.hpp>
    
template <typename Matrix> struct darcy3d1d_precond{
    
    //TODO add mixed bc for tissue !!!!!
    
    //gmm::diagonal_precond <Matrix> D;
    //Matrix schur; 
    
    darcy_precond <Matrix> Prec_t;
    darcy_precond <Matrix> Prec_v;
    
    getfem::size_type dof_ut;
    getfem::size_type dof_pt;
    getfem::size_type dof_uv;
    getfem::size_type dof_pv;
    
    //getfem::size_type nrows = dof_u + dof_p;
    //getfem::size_type ncols = dof_u + dof_p;
    
    
    // default constructor
    darcy3d1d_precond(void) {} 
    //constructor
    darcy3d1d_precond(const Matrix &A, const getfem::mesh_fem & mf_ut, const getfem::mesh_fem & mf_pt, const getfem::mesh_im & mimt,
        const getfem::mesh_fem & mf_uv, const getfem::mesh_fem & mf_pv, const getfem::mesh_im & mimv)
        {
        dof_ut = mf_ut.nb_dof();
        dof_pt = mf_pt.nb_dof();
        dof_uv = mf_uv.nb_dof();
        dof_pv = mf_pv.nb_dof();
        
        // build the tissue precond 
        Prec_t.build_with(gmm::sub_matrix(A,
                                         gmm::sub_interval(0, dof_ut + dof_pt), gmm::sub_interval(0, dof_ut + dof_pt)),
        
            mf_ut, mf_pt, mimt);
        // build the vessel precond
        Prec_v.build_with(gmm::sub_matrix(A,
                                         gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv), gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv)),
        
            mf_uv, mf_pv, mimv);        
    } 
    
    }; // end of struct
    
    
// Multiplication of the prec by a vector: P*vec -> res.
// In our case we build P but what we have to do is P^-1*vec, therefore we solve the linear system P res =  vec 

template < typename Matrix, typename V1, typename V2>
void mult(const darcy3d1d_precond<Matrix> &P, const V1 &vec, V2 &res) {
    
    // multiplication of the tissue block
    mult(P.Prec_t, gmm::sub_vector(vec, (0, P.dof_ut + P.dof_pt)), gmm::sub_vector(res, (0, P.dof_ut + P.dof_pt)));

    // multiplication of the vessel block
    mult(P.Prec_v, gmm::sub_vector(vec, (0, P.dof_uv + P.dof_pv)), gmm::sub_vector(res, (0, P.dof_uv + P.dof_pv)));
    
} // end mult

//Redefinition of some gmm functions for darcy_precond

namespace gmm{
    
/*struct owned_implementation {};

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
*/
}    
 // end of namespace gmm
