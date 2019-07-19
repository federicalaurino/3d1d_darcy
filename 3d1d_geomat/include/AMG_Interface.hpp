// AMG_Interrface.h
//#ifdef WITH_SAMG
//#ifndef aamg
//#define aamg


#include "gmm/gmm.h"
#include "getfem/getfem_mesh.h"

/* default 4 Byte integer types */
#ifndef APPL_INT
#define APPL_INT int
#endif

using bgeot::scalar_type;
using bgeot::size_type;

class AMG {
    
private:
    
std::vector<scalar_type> a_samg_;   // stores the entries of the matrix in SAMG format
std::vector<unsigned int> ja_samg_;    // stores the column indices of the entries of a_samg_  
std::vector<unsigned int> ia_samg_;    // row repartition on a_samg_ and ja_samg_
std::vector<scalar_type> sol; //samg solution   


public:
//NOTE I thinnk it should not be necessery to pass also U_
gmm::csr_matrix <scalar_type> A_csr_;
std::vector <scalar_type> U_;
std::vector <scalar_type> F_;

  // ======== costructor the class ========================
  AMG(std::string name);
  AMG(std::string name, gmm::csr_matrix<scalar_type> A_csr, 
      std::vector <scalar_type> U_, std::vector <scalar_type> F_);
    // ======== destructor the class ========================
  ~AMG();  // This is the destructor: declaration
  // ======== generation af matrix
  void csr2samg(void);
    // ======== solver of the class ========================
  void solve(void);  
    // =========== return the solution ========================
  std::vector<scalar_type> getsol(){return sol;}
};

//#endif
//  #endif // WITH_SAMG
