/*  Class for darcy preconditioner for the 3D1D darcy problem:
    P = [diag_t^-1 0  0 0; 0 schur_t^-1; 0 0 diag_v^-1 schur_v^-1]
    where diag is the diagonal of the velocity mass matrix and 
    schur is an approximation of the schur complement given by the jumps of the pressure
    across the internal faces + a boundary term
*/

#include <darcy_precond.hpp>
    
template <typename Matrix> struct darcy3d1d_precond{
    
    //TODO add the possibility to have mixed bc !!!!!
    
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
        const getfem::vector<getfem::mesh_fem> & mf_uvi, const getfem::mesh_fem & mf_pv, const getfem::mesh_im & mimv)
        {
        dof_ut = mf_ut.nb_dof();
        dof_pt = mf_pt.nb_dof();
        dof_uv = 0;
        for (gmm::size_type i = 0; i<mf_uvi.size(); i++)
            dof_uv += mf_uvi[i].nb_dof();
        dof_pv = mf_pv.nb_dof();
        
        Matrix M;
        gmm::copy(gmm::sub_matrix(A,
                        gmm::sub_interval(0, dof_ut + dof_pt), gmm::sub_interval(0, dof_ut + dof_pt)), M);
        // build the tissue precond 
        Prec_t.build_with(M, mf_ut, mf_pt, mimt);
        // build the vessel precond
        gmm::clear(M);
        gmm::copy(gmm::sub_matrix(A,
                                         gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv), gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv)), M);
        Prec_v.build_with(M, mf_uvi, mf_pv, mimv);        
    } 
    
    // constructor for preconditioner with conditions for tissue
        darcy3d1d_precond(const Matrix &A, const getfem::mesh_fem & mf_ut, const getfem::mesh_fem & mf_pt, const getfem::mesh_im & mimt,
                          const getfem::vector_type k,
                        const getfem::vector<getfem::mesh_fem> & mf_uvi, const getfem::mesh_fem & mf_pv, const getfem::mesh_im & mimv)
        {
        dof_ut = mf_ut.nb_dof();
        dof_pt = mf_pt.nb_dof();
        dof_uv = 0;
        for (gmm::size_type i = 0; i<mf_uvi.size(); i++)
            dof_uv += mf_uvi[i].nb_dof();
        dof_pv = mf_pv.nb_dof();
        
        Matrix M;
        gmm::copy(gmm::sub_matrix(A,
                        gmm::sub_interval(0, dof_ut + dof_pt), gmm::sub_interval(0, dof_ut + dof_pt)), M);
        // build the tissue precond with mixed condition 
        Prec_t.build_with(M, mf_ut, mf_pt, mimt, k);
        // build the vessel precond
        gmm::clear(M);
        gmm::copy(gmm::sub_matrix(A,
                                         gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv), gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv)), M);
        Prec_v.build_with(M, mf_uvi, mf_pv, mimv);        
    }
    
    // constructor for preconditioner with mixed conditions for vessel
        darcy3d1d_precond(const Matrix &A, const getfem::mesh_fem & mf_ut, const getfem::mesh_fem & mf_pt, const getfem::mesh_im & mimt,
                        const getfem::vector<getfem::mesh_fem> & mf_uvi, const getfem::mesh_fem & mf_pv, const getfem::mesh_im & mimv,
                        const getfem::vector_type k
                         )
        {
        dof_ut = mf_ut.nb_dof();
        dof_pt = mf_pt.nb_dof();
        dof_uv = 0;
        for (gmm::size_type i = 0; i<mf_uvi.size(); i++)
            dof_uv += mf_uvi[i].nb_dof();
        dof_pv = mf_pv.nb_dof();
        
        Matrix M;
        gmm::copy(gmm::sub_matrix(A,
                        gmm::sub_interval(0, dof_ut + dof_pt), gmm::sub_interval(0, dof_ut + dof_pt)), M);
        // build the tissue precond 
        Prec_t.build_with(M, mf_ut, mf_pt, mimt);
        // build the vessel precond with mixed condition
        gmm::clear(M);
        gmm::copy(gmm::sub_matrix(A,
                                         gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv), gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv)), M);
        Prec_v.build_with(M, mf_uvi, mf_pv, mimv, k);        
    }
    
// constructor for preconditioner with mixed conditions for vessel and tissue
        darcy3d1d_precond(const Matrix &A, const getfem::mesh_fem & mf_ut, const getfem::mesh_fem & mf_pt, const getfem::mesh_im & mimt,
                           const getfem::vector_type kt,
                        const getfem::vector<getfem::mesh_fem> & mf_uvi, const getfem::mesh_fem & mf_pv, const getfem::mesh_im & mimv,
                        const getfem::vector_type kv
                         )
        {
        dof_ut = mf_ut.nb_dof();
        dof_pt = mf_pt.nb_dof();
        dof_uv = 0;
        for (gmm::size_type i = 0; i<mf_uvi.size(); i++)
            dof_uv += mf_uvi[i].nb_dof();
        dof_pv = mf_pv.nb_dof();
        
        Matrix M;
        gmm::copy(gmm::sub_matrix(A,
                        gmm::sub_interval(0, dof_ut + dof_pt), gmm::sub_interval(0, dof_ut + dof_pt)), M);
        // build the tissue precond with mixed condition
        Prec_t.build_with(M, mf_ut, mf_pt, mimt, kt);
        // build the vessel precond with mixed condition
        gmm::clear(M);
        gmm::copy(gmm::sub_matrix(A,
                                         gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv), gmm::sub_interval(dof_ut + dof_pt, dof_uv + dof_pv)), M);
        Prec_v.build_with(M, mf_uvi, mf_pv, mimv, kv);        
    }
    
    
    }; // end of struct
    
    
// Multiplication of the prec by a vector: P*vec -> res.
// In our case we build P but what we have to do is P^-1*vec, therefore we solve the linear system P res =  vec 

template < typename Matrix, typename V1, typename V2>
void mult(const darcy3d1d_precond<Matrix> &P, const V1 &vec, V2 &res) {
    
    // multiplication of the tissue block
    std::vector<gmm::size_type> ipvt(P.dof_ut + P.dof_pt);
    for(gmm::size_type i=0; i< P.dof_ut + P.dof_pt; i++ ) ipvt[i] = i;
    gmm::unsorted_sub_index interval;
    interval = gmm::unsorted_sub_index(ipvt);
    mult(P.Prec_t, gmm::sub_vector(vec, interval), gmm::sub_vector(res, interval));

    // multiplication of the vessel block
    std::vector<gmm::size_type> ipvt1(P.dof_uv + P.dof_pv);
    for(gmm::size_type i=P.dof_ut + P.dof_pt; i< P.dof_ut + P.dof_pt + P.dof_uv + P.dof_pv; i++ ) ipvt1[i] = i;
    gmm::unsorted_sub_index interval1;
    interval1 = gmm::unsorted_sub_index(ipvt1);
    mult(P.Prec_v, gmm::sub_vector(vec, interval1), gmm::sub_vector(res, interval1));
    
} // end mult

