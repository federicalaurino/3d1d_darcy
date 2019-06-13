//#ifdef WITH_SAMG
#include "AMG_Interface.hpp" 
//#include <memory>

AMG::AMG(std::string name)
{
	std::cout<<"Build class AMG for "<< name << std::endl;
		
}

AMG::AMG(std::string name, gmm::csr_matrix<scalar_type> A_csr, 
      std::vector <scalar_type> F, std::vector <scalar_type> U);)
{
    std::cout<<"Build class AMG for "<< name << std::endl;
    gmm::copy(A_csr, A_csr_);
    gmm::resize(F_, F.size());
    gmm::copy(F, F_);
    gmm::resize(U_, U.size());
    gmm::copy(U, U_);
}

AMG::~AMG(void) {
	// delete [] 
	std::cout << "AMG deleted" << std::endl;
		
}

// Convert a csr matrix into SAMG format

void AMG::csr2samg(void)
{
    std::cout << "Conversion to SAMG format" << std::endl;
    
    gmm::resize(a_samg_, A_csr_.pr.size()); gmm::clear(a_samg_);
    gmm::copy(A_csr_.pr, a_samg_);
    gmm::resize(ja_samg_, A_csr_.ir.size()); gmm::clear(ja_samg_);
    gmm::copy(A_csr_.ir, ja_samg_);
    gmm::resize(ia_samg_, A_csr_.jc.size()); gmm::clear(ia_samg_);
    gmm::copy(A_csr_.jc, ia_samg_);
        
    
    size_type nnu = ia_samg_.size()-1;
    for (size_type row=0; row < nnu; row ++)
    {   
        size_type first = ia_samg_[row];     // index of the first element of row
        size_type last = ia_samg_[row+1]-1;    // index of the last element of row
        
        std::cout << "check diag" << std::endl; 
        //check if row has a nonnull diag element
        bool has_diag_el = false;  
        //TODO can I improve it with iterators? not to have so many for cycles
        size_type k=first;
        while (k<=last && ja_samg_[k] != row)
            k++;
        if (ja_samg_[k] == row)
            has_diag_el=true;
        
        std::cout << "if diag" << std::endl;
        if(has_diag_el && k!= first)
        {   
            //Move the diagonal element as first element of the row
            std::swap(a_samg_[k], a_samg_[first]);
            std::swap(ja_samg_[k], ja_samg_[first]);
        }
        
        std::cout << " if not diag" << std::endl;
        if(!has_diag_el){
            //Add 0 diag element to the element vector in first position of row
            auto it1= a_samg_.begin();
            a_samg_.insert(it1+first,0);
            //Add the column index (=row) of the diag element
            auto it2=ja_samg_.begin();
            ja_samg_.insert(it2+first,row);
            //increment by 1 the 'starting point' of the following rows 
            for (size_type i=row+1; i<ia_samg_.size(); i++)
                ia_samg_[i]=ia_samg_[i]+1;
        }

    }// end for on the rows

	return;
}



void AMG::solve(void)
{
        
    // ===> Set primary parameters. Others can be set by access functions as shown below.
    
    int nnu = ia_samg_.size() -1; // number of unknowns;
    int nna = a_samg_.size(); // number of nnz entries
    std::vector <unsigned int> ia (ia_samg_.size());
    gmm::copy(ia_samg_, ia);
    std::vector <unsigned int> ja (ja_samg_.size());
    gmm::copy(ja_samg_, ja);
    std::vector <scalar_type> a (a_samg_.size());
    gmm::copy(a_samg_, a);
    std::vector <scalar_type> f (F_.size());
    gmm::copy(F_, f);
    std::vector <scalar_type u (F_.size()); gmm::clear(u);
    

    int nsys = 1; 
    
    // int matrix = X Y where
    //   X = 1 if A is symmetric, full matrix stored (normal case)
    //     = 2 A is not symmetric.   
    //     = 3 A is symmetric, lower triangular part stored (s. Remark below!)
    //   Y = 1 if A is a zero row sum matrix 11 . 
    //     = 2 A is not a zero row sum matrix 12 
    
    int matrix = 1; 
    int ifirst = 1;        // first approx for u : =1 (first approx = zero), =2 (first approx = 1)
    double eps  = 1.0e-8;   // required (relative) residual reduction
    std::vector<int> iscale(nsys, 0); // indicating which unknowns require scaling.
    
    std::vector<int> iu (1,0);        //  (dummy) vector iu
    std::vector<int> ndip (1,0);        // (dummy) vector ip    
    int    ndiu      = 1;        // dimension of (dummy) vector iu
    int    ndip      = 1;        // dimension of (dummy) vector ip
    
    
    // ===> Selecting the solution approach
    
      int   nsolve    = 2;        // results in scalar approach (current system is scalar)

      

    // ===> Selecting SAMG cycling process  
      
      int    ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations
      int    n_default = 20;       // select default settings for secondary parameters
                                   // CURRENTLY AVAILABLE: 10-13, 15-18, 20-23, 25-28
                                   // NOTE: the higher the respective SECOND digit, the
                                   // more aggressive the coarsening (--> lower memory at
                                   // the expense of slower convergence)
      int    iswtch    = 5100+n_default; // complete SAMG run ....
                                   // ... memory de-allocation upon return ....
                                   // ... memory extension feature activated ....
                                   // ... residuals measured in the L2-norm .....
                                   // ... secondary parameter default setting # n_default

      double a_cmplx   = 2.2;      // estimated dimensioning
      double g_cmplx   = 1.7;      // estimated dimensioning
      double w_avrge   = 2.4;      // estimated dimensioning
      double p_cmplx   = 0.0;      // estimated dimensioning (irrelevant for scalar case)

      double chktol    = -1.0;     // input checking de-activated (we know it's ok!)
      int    idump     = 0;        // minimum output during setup
      int    iout      = 2;        // display residuals per iteration and work statistics

      
        
      
        
        
        samg(nnu,nna,nsys,
            ia,ja,a,f,u,iu,ndiu,ip,ndip,matrix,iscale,
            res_in,res_out,ncyc_done,ierr,
            nsolve,ifirst,eps,ncyc,iswtch,
            a_cmplx,g_cmplx,p_cmplx,w_avrge,
            chktol,idump,iout);
        
	return;
}*/
//===============================================================
//===============================================================
//===============================================================
//===============================================================

/*
void AMG::set_pt2uk(int * dofpt , int q_dof, int l_dof, int npts){
	
		std::std::cout << std::string(100, '=') << std::std::endl;
		std::std::cout<< "AMG::set_pt2uk start"  << std::std::endl;
		std::std::cout << std::string(100, '=') << std::std::endl;
		
		
		_q_dof=q_dof; _l_dof=l_dof; _npts=npts;
		_pt2uk.resize(_q_dof + l_dof);
		for (int idof=0; idof< q_dof+l_dof; idof++)  _pt2uk[idof] = dofpt[idof] ;
			
				
				std::std::cout << std::string(100, '=') << std::std::endl;
				std::std::cout<< "AMG::set_pt2uk end "  << std::std::endl;
				std::std::cout << std::string(100, '=') << std::std::endl;
}

void AMG::set_dof(int pt , int ut, int pv, int uv){
	_Pt=pt;_Ut=ut;_Pv=pv; _Uv=uv;
}

*/


//#endif // if def with samg









