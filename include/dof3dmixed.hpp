/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for the number of degrees of freedom.
 */
#ifndef M3D1D_DOF3D1D_TRANSP_HPP_
#define M3D1D_DOF3D1D_TRANSP_HPP_

namespace getfem {

//! Class to store the number of degrees of freedom of used FEMs
struct dof3dmixed {


	//! Number of dof of the interstitial pressure FEM mf_Pt
	size_type Pt_;
	//! Number of dof of the vessel pressure FEM mf_Ut
	size_type Ut_;
    //! Total number of dof
    size_type tot_;
    //! Number of dof of the tissue data FEM mf_coeft
	size_type coeft_;

	
	//! Compute the number of dof of given FEM
	void set (
            const getfem::mesh_fem & mf_Ut, 
			const getfem::mesh_fem & mf_Pt, 
            const getfem::mesh_fem & mf_coeft
             )
	{   Ut_ = mf_Ut.nb_dof();
		Pt_ = mf_Pt.nb_dof();
        coeft_=mf_coeft.nb_dof();
        tot_= Ut_ + Pt_;
	}

	//! Accessor to the number of dof of mf_Pt
	// TODO why I need accessors?
	inline size_type Pt (void) { return Pt_; } const
	//! Accessor to the number of dof of mf_Pv
	inline size_type Ut (void) { return Ut_; } const
	//! Accessor to the number of dof of mf_coeft
	inline size_type coeft (void) { return coeft_; } const	
	//! Accessor to the total nummber of dof
	inline size_type tot (void) { return tot_; } const

	
	//! Overloading of the output operator
	friend std::ostream & operator << (
			std::ostream & out, const dof3dmixed & dof
			)
	{ 
		out << "--- DEGREES OF FREEDOM --- " << endl;
        out << "  nb_dof_Ut     : " 		 << dof.Ut_  << endl;
		out << "  nb_dof_Pt     : " 		 << dof.Pt_  << endl;
        out << "  nb_dof_coeft     : " 		 << dof.coeft_  << endl;
		out << "-------------------------- " << endl;
		return out;            
	}

};

} /* end of namespace */

#endif
