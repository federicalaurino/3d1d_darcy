#include "gmm/gmm.h"
#include "getfem/getfem_mesh.h"

using bgeot::scalar_type;
using bgeot::size_type;

class AMG {
    
private:

  gmm::csr_matrix <scalar_type> A_csr_; // csr matrix defining the linear system
  std::vector <scalar_type> F_; // rhs defining the linear system
  std::vector<scalar_type> a_samg_;   // stores the entries of the matrix in SAMG format
  std::vector<unsigned int> ja_samg_;    // stores the column indices of the entries of a_samg_  
  std::vector<unsigned int> ia_samg_;    // row repartition on a_samg_ and ja_samg_
  std::vector<scalar_type> sol; //samg solution   

public:

    //! Costructor the class
    AMG();
    AMG(gmm::csr_matrix<scalar_type> A_csr, std::vector <scalar_type> F_);

    //! Destructor the class 
    ~AMG();  // This is the destructor: declaration
    //! Converting the matrix from csr to samg format
    void csr2samg(void);
    //! SAMG solver
    void solve(void);  
    //! Getting the solution
    std::vector<scalar_type> getsol(){return sol;}
};


