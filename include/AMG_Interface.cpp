//#ifdef WITH_SAMG
#include "AMG_Interface.hpp" 
#include <iostream>
#include <fstream>
//#include <memory>
#include "samg.h"   
AMG::AMG(std::string name)
{
	std::cout<<"Build class AMG for "<< name << std::endl;
		
}

AMG::AMG(std::string name, gmm::csr_matrix<scalar_type> A_csr, 
      std::vector <scalar_type> U, std::vector <scalar_type> F)
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
        
 
        //check if row has a nonnull diag element
        bool has_diag_el = false;  
        //TODO can I improve it with iterators? not to have so many for cycles
        size_type k=first;
        while (k<=last && ja_samg_[k] != row)
            k++;
        if (ja_samg_[k] == row)
            has_diag_el=true;
        
        
        if(has_diag_el && k!= first)
        {   
            //Move the diagonal element as first element of the row
            std::swap(a_samg_[k], a_samg_[first]);
            std::swap(ja_samg_[k], ja_samg_[first]);
        }
        
 
        if(!has_diag_el){
            //Add 0 diag element to the element vector in first position of row
            auto it1 = a_samg_.begin();
            a_samg_.insert(it1+first,0);
            //Add the column index (=row) of the diag element
            auto it2 = ja_samg_.begin();
            ja_samg_.insert(it2+first,row);
            //increment by 1 the 'starting point' of the following rows 
            for (size_type i=row+1; i<ia_samg_.size(); i++)
                ia_samg_[i] = ia_samg_[i]+1;
        }

    }// end for on the rows
    
    //Covert to Fortran: indices start from 1...
    //TODO can I improve it with iterators? not to have so many for cycles
    for(int k=0; k<ia_samg_.size(); k++) ia_samg_[k] +=1; 
    for(int k=0; k<ja_samg_.size(); k++) ja_samg_[k] +=1;

	return;
}



void AMG::solve(void)
{
  
    // ===> Set primary parameters. Others can be set by access functions as shown below.

      APPL_INT    ndiu      = 1;        // dimension of (dummy) vector iu
      //TODO Add case in which nsys !=1 and ndiu is not trivial
      APPL_INT    ndip      = 1;        // dimension of (dummy) vector ip

      APPL_INT    nsolve    = 2;        // results in scalar approach (current system is scalar)
      APPL_INT    ifirst    = 1;        // first approximation = zero
      double eps       = 1.0e-8;   // required (relative) residual reduction
      APPL_INT    ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations
      /* If we want to use Samg as a stand-alone solver (not as a preconditioner), both ncgrad (the "nd number in ncyc) and ncgrad_default must be equal to 0
 				
      APPL_INT ncgrad_default=0;
      SAMG_SET_NCGRAD_DEFAULT(&ncgrad_default);
      ncyc      = 10050;    // V-cycle as pre-conditioner for CG; at most 50 iterations
      */
   
      double a_cmplx   = 2.2;      // estimated dimensioning
      double g_cmplx   = 1.7;      // estimated dimensioning
      double w_avrge   = 2.4;      // estimated dimensioning
      double p_cmplx   = 0.0;      // estimated dimensioning (irrelevant for scalar case)

      double chktol    = -1.0;     // input checking de-activated (we know it's ok!)
      
      //============================================
      //     idump controls the matrix dumping of SAMG				
      // 1  Standard prAPPL_INT output, no matrix dump. 
      // 2‐6  Write matrices to disk: level 2 up to level idmp. 
      // 7  Write matrices to disk: level 2 up to the coarsest level. 
      // 8  Write finest‐level matrix to disk (incl. right hand side etc.). 
      // 9  Write all matrices to disk.       
      APPL_INT    idump     = 0;        // minimum output during setup

	  // iout page 44 Userguide. it controls display outpu. default 2 very verbose 43			
      APPL_INT    iout      = 43;        // display residuals per iteration and work statistics
      
            
      APPL_INT    n_default = 20;       // select default settings for secondary parameters
                                   // CURRENTLY AVAILABLE: 10-13, 15-18, 20-23, 25-28
                                   // NOTE: the higher the respective SECOND digit, the
                                   // more aggressive the coarsening (--> lower memory at
                                   // the expense of slower convergence)
      APPL_INT    iswtch    = 5100+n_default; // complete SAMG run ....
                                   // ... memory de-allocation upon return ....
                                   // ... memory extension feature activated ....
                                   // ... residuals measured in the L2-norm .....
                                   // ... secondary parameter default setting # n_default

    // ===> Secondary parameters which have to be set if n_default=0
    //      (at the same time a demonstration of how to access secondary or hidden parameters)

      APPL_INT intin=0;
      double dblin=0;

      if (n_default == 0) {
         intin=25;    SAMG_SET_LEVELX(&intin);
         intin=100;   SAMG_SET_NPTMN(&intin);
         intin=4;     SAMG_SET_NCG(&intin);
         intin=2;     SAMG_SET_NWT(&intin);
         intin=1;     SAMG_SET_NTR(&intin);
         intin=131;   SAMG_SET_NRD(&intin);
         intin=131;   SAMG_SET_NRU(&intin);
         intin=0;     SAMG_SET_NRC(&intin);
         intin=0;     SAMG_SET_NP_OPT(&intin);

         dblin=21.25; SAMG_SET_ECG(&dblin);
         dblin=0.20;  SAMG_SET_EWT(&dblin);
         dblin=12.20; SAMG_SET_ETR(&dblin);
      }

    // ===> amg declarations 

    // input:

    APPL_INT npnt,nsys,matrix,nnu,nna;

    APPL_INT * ia, * ja;
    APPL_INT * iu, * ip, * iscale;

    double * a, * u, * f;

    // output:
    APPL_INT ierr,ierrl,ncyc_done;
    double res_out,res_in;

    nnu = ia_samg_.size() -1; // number of unknowns;
    nna = a_samg_.size(); // number of nnz entries
      
    // APPL_INT matrix = XY where
    //   X = 1 if A is symmetric, full matrix stored (normal case)
    //     = 2 A is not symmetric.   
    //     = 3 A is symmetric, lower triangular part stored (s. Remark below!)
    //   Y = 1 if A is a zero row sum matrix 11 . 
    //     = 2 A is not a zero row sum matrix 12 
    
    matrix = 22;
    nsys = 1;
    npnt = 0;

    if (nsys==1) 
        iu= new int[1];
    //TODO case nsys != 1 
    /*else {
        iu  = new APPL_INT[nnu];
        ndiu   = nnu;
        for(APPL_INT iiu=0;iiu<nnu;iiu++){
            //TODO only working for nsys = 2 check it out for larger systems
            if(iiu<dof_darcy.Pt())iu[iiu]=1;
            else iu[iiu]=2;
            }//end for iiu
        }*/
        
    ip = new int[1];
   
    ia = new int[nnu+1];
    ja = new int[nna];
    a  = new double[nna];
    if (!(ia && ja && a)) {
        std::cout << " allocation failed (ia,ja,a) " << std::endl;
        ierr=1;// return ierr;
    }
    int i;
    for (i=0;i<nnu+1;i++)   ia[i] = ia_samg_[i];
    for (i=0;i<nna;i++) a[i] = a_samg_[i]; 
    for (i=0;i<nna;i++) ja[i] = ja_samg_[i];
    
    f = new double[nnu]; u = new double[nnu]; 
    if (!(f && u)) {
        std::cout << " allocation failed (f,u) " << std::endl;
        ierr=1; //return ierr;
    }
    for (i=0;i<nnu;i++) u[i]=0.0;
    for (i=0;i<nnu;i++) f[i] = F_[i];
    
    // this vector (iscale) indicates which uknowns require scaling if 0 no scaling
    iscale = new int[nsys]; for(APPL_INT i_sys=0; i_sys<nsys; i_sys++) iscale[i_sys]=0;  
 
    /*//****************************
    //checking samg vectors...
    std::ofstream s_var("samg_variables.txt");
    s_var << " a_samg_ " ;
    for (i=0; i < nna; i++) s_var << " " << a_samg_[i] ;
    s_var << "\n";
    s_var << " *a " ;
    for (i=0; i < nna; i++) s_var << " " << a[i];
    s_var << " \n " ;
    
    s_var << " ia_samg_ " ;
    for (i=0; i < nnu+1; i++) s_var << " " << ia_samg_[i] ;
    s_var << "\n";
    s_var << " *ia " ;
    for (i=0; i < nnu+1; i++) s_var << " " << ia[i];
    s_var << " \n " ;
    
    s_var << " ja_samg_ " ;
    for (i=0; i < nna; i++) s_var << " " << ja_samg_[i] ;
    s_var << "\n";
    s_var << " *ja " ;
    for (i=0; i < nna; i++) s_var << " " << ja[i];
    s_var << " \n " ;
    
    s_var << " F_ " ;
    for (i=0; i < nnu; i++) s_var << " " << F_[i] ;
    s_var << "\n";
    s_var << " *f " ;
    for (i=0; i < nnu; i++) s_var << " " << f[i];
    s_var << " \n " ;
    
    s_var << " U_ " ;
    for (i=0; i < nnu; i++) s_var << " " << U_[i] ;
    s_var << "\n";
    s_var << " *u " ;
    for (i=0; i < nnu; i++) s_var << " " << u[i];
    s_var << " \n " ;
    s_var.close();
    */ 
    //***************************************

// ===> if, e.g., the finest-level matrix shall be dumped:
//
//    char *string ="myfilename";
//    APPL_INT length=10;
//    SAMG_SET_FILNAM_DUMP(&string[0],&length,&length);
//    idump=8;

// ===> Call FORTRAN90 SAMG

//    SAMG_RESET_SECONDARY();   // necessary before a second SAMG run
	                        // if secondary parameters have to be reset
	                        // see manual.
  
    //=============================================
    // Choose the type of solver

    // Use samg direct solver
    APPL_INT levelx;
    SAMG_GET_LEVELX(&levelx); // retreive levelx=number of maximum coarsening levels 
        std::cout<<"..value of levelx======"<<levelx<<std::endl;
    levelx=1;
    SAMG_SET_LEVELX(&levelx);// change levelx=number of maximum coarsening levels 
    SAMG_GET_LEVELX(&levelx);// retreive levelx=number of maximum coarsening levels 
    std::cout<<"..check if change value of levelx======"<<levelx<<std::endl;
    
    APPL_INT nrc=11;
    APPL_INT  nrc_emergency=11;
    APPL_INT nptmax=5000;
    APPL_INT clsolver_finest=1;
    SAMG_SET_NRC(&nrc);// change  nrc=solver for coarser levels page 67 user guide
    SAMG_SET_NRC_EMERGENCY(&nrc_emergency);
    SAMG_SET_NPTMAX(&nptmax);// change  nptmax
    SAMG_SET_CLSOLVER_FINEST(&clsolver_finest);// change  nptmax
    ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations

/**    
    //AMG stand alone
    //Both ncgrad (the "nd number in ncyc) and ncgrad_default must be equal to 0 to use SAMG solver as a stand-alone solver (not as a preconditioner)

    APPL_INT ncgrad_default=0;
    SAMG_SET_NCGRAD_DEFAULT(&ncgrad_default);
    ncyc      = 10050;    // V-cycle as pre-conditioner for CG; at most 50 iterations

        
    // AMG accelerated
    ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations	
**/
        
    //=====================================================
    // printing format
    char ch[] = "f";int len=1;
    SAMG_SET_IOFORM(&ch[0],&len);
    
    float told,tnew,tamg;
    SAMG_CTIME(&told);
     
    SAMG(&nnu,&nna,&nsys,
        &ia[0],&ja[0],&a[0],&f[0],&u[0],&iu[0],&ndiu,&ip[0],&ndip,&matrix,&iscale[0],
        &res_in,&res_out,&ncyc_done,&ierr,
        &nsolve,&ifirst,&eps,&ncyc,&iswtch,
        &a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
        &chktol,&idump,&iout);

    if (ierr > 0) {
        std::cout << std::endl << " SAMG terminated with error code " 
            << ierr << " **** " << std::endl;
    }
    else if (ierr < 0) {
        std::cout << std::endl << " SAMG terminated with warning code " 
            << ierr << " **** " << std::endl;
    }

    SAMG_CTIME(&tnew);
    tamg=tnew-told;
    std::cout << std::endl << " ***** total run time: " << tamg << " ***** " << std::endl; 
      
    SAMG_LEAVE(&ierrl);
      
    if (ierrl != 0) {
        std::cout << std::endl << " error at samg_leave" 
            << ierr << " **** " << std::endl;
    }
    
    //delete[] ia,ja,a,f,u,iu,ip,iscale;
    
    //save the solution
    gmm::resize(sol, nnu);
    for (int i=0; i<nnu ; i++)
        sol[i]=u[i];
    
    delete[] ia,ja,a,f,u,iu,ip,iscale;
    
    return;
}
//===============================================================
//===============================================================
//===============================================================
//===============================================================

/*
void AMG::set_pt2uk(APPL_INT * dofpt , APPL_INT q_dof, APPL_INT l_dof, APPL_INT npts){
	
		std::std::cout << std::string(100, '=') << std::std::endl;
		std::std::cout<< "AMG::set_pt2uk start"  << std::std::endl;
		std::std::cout << std::string(100, '=') << std::std::endl;
		
		
		_q_dof=q_dof; _l_dof=l_dof; _npts=npts;
		_pt2uk.resize(_q_dof + l_dof);
		for (APPL_INT idof=0; idof< q_dof+l_dof; idof++)  _pt2uk[idof] = dofpt[idof] ;
			
				
				std::std::cout << std::string(100, '=') << std::std::endl;
				std::std::cout<< "AMG::set_pt2uk end "  << std::std::endl;
				std::std::cout << std::string(100, '=') << std::std::endl;
}

void AMG::set_dof(APPL_INT pt , APPL_INT ut, APPL_INT pv, APPL_INT uv){
	_Pt=pt;_Ut=ut;_Pv=pv; _Uv=uv;
}

*/


//#endif // if def with samg









