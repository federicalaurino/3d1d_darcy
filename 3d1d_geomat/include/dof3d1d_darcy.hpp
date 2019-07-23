/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
3d1d Darcy problem - IJGE
======================================================================*/
/*! 
  @file   AMG_Interface.cpp
  @author Federica Laurino <federica.laurino@polimi.it>
  @date   2019.
  @brief  Definition of the aux class for the number of degrees of freedom.
 */

#ifndef M3D1D_DOF3D1D_HPP_
#define M3D1D_DOF3D1D_HPP_

namespace getfem {

//! Class to store the number of degrees of freedom of used FEMs
struct dof3d1d_darcy {


	//! Number of dof of the interstitial pressure FEM mf_Pt
	size_type Pt_;
	//! Number of dof of the vessel pressure FEM mf_Pv
	size_type Pv_;
    //! Total number of dof
    size_type tot_;
    //! Number of dof of the tissue data FEM mf_coeft
	size_type coeft_;
	//! Number of dof of the vessel data FEM mf_coefv
	size_type coefv_;
	
	//! Compute the number of dof of given FEM
	void set (
			const getfem::mesh_fem & mf_Pt, const getfem::mesh_fem & mf_Pv,
            const getfem::mesh_fem & mf_coeft, const getfem::mesh_fem & mf_coefv
             )
	{
		Pt_ = mf_Pt.nb_dof();
		Pv_ = mf_Pv.nb_dof();
        coeft_=mf_coeft.nb_dof();
        coefv_=mf_coefv.nb_dof();
        tot_= Pt_ + Pv_;
	}

	//! Accessor to the number of dof of mf_Pt
	// TODO why I need accessors?
	inline size_type Pt (void) { return Pt_; } const
	//! Accessor to the number of dof of mf_Pv
	inline size_type Pv (void) { return Pv_; } const
	//! Accessor to the number of dof of mf_coeft
	inline size_type coeft (void) { return coeft_; } const
	//! Accessor to the number of dof of mf_coefv
	inline size_type coefv (void) { return coefv_; } const	
	//! Accessor to the total nummber of dof
	inline size_type tot (void) { return tot_; } const

	
	//! Overloading of the output operator
	friend std::ostream & operator << (
			std::ostream & out, const dof3d1d_darcy & dof
			)
	{ 
		out << "--- DEGREES OF FREEDOM --- " << endl;
		out << "  nb_dof_Pt     : " 		 << dof.Pt_  << endl;
		out << "  nb_dof_Pv     : " 		 << dof.Pv_  << endl;
        out << "  nb_dof_coeft     : " 		 << dof.coeft_  << endl;
		out << "  nb_dof_coefv     : " 		 << dof.coefv_  << endl;
		out << "-------------------------- " << endl;
		return out;            
	}

};

} /* end of namespace */

#endif
