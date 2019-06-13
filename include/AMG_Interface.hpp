// AMG_Interrface.h
//#ifdef WITH_SAMG
//#ifndef aamg
//#define aamg


#include "gmm/gmm.h"
#include "getfem/getfem_mesh.h"
//#include "samg.h"
/* default 4 Byte integer types */
//#ifndef APPL_INT
//#define APPL_INT int
//#endif

using bgeot::scalar_type;
using bgeot::size_type;

class AMG {
    
  //private:
/*std::vector<scalar_type> sol_vec; /// solution vector
std::vector<int> _pt2uk; /// point to oknow vector
int _q_dof;              /// quadratic dof
int _l_dof;				/// linear dof
int _Pt, _Ut, _Pv, _Uv;
int _npts;
bool first_=true;*/



private:
std::vector<scalar_type> a_samg_;   // stores the entries of the matrix in SAMG format
std::vector<unsigned int> ja_samg_;    // stores the column indices of the entries of a_samg_  
std::vector<unsigned int> ia_samg_;    // row repartition on a_samg_ and ja_samg_
//APPL_INT nnu_,nna_;

public:
  
gmm::csr_matrix <scalar_type> A_csr_;
std::vector <scalar_type> U_;
std::vector <scalar_type> F_;

  // ======== costructor the class ========================
  AMG(std::string name);
  AMG(std::string name, gmm::csr_matrix<scalar_type> A_csr, 
      std::vector <scalar_type> U_, std::vector <scalar_type> F_);
// void set_dof(int , int, int, int);
    // ======== destructor the class ========================
  ~AMG();  // This is the destructor: declaration
  // ======== generation af matrix
  void csr2samg(void);
    // ======== solver of the class ========================
  void solve(void);
  // =========== set to point to uknown vector ========================
//  void set_pt2uk(int * dofpt , int q_dof, int l_dof, int npts);
  
    // =========== return the solution ========================
//  std::vector<scalar_type> getsol(){return sol_vec;}
};

//#endif
//  #endif // WITH_SAMG
