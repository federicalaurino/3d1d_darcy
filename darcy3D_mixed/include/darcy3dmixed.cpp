/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
  "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
  Course on Advanced Programming for Scientific Computing
  Politecnico di Milano
  A.Y. 2016-2017


  Copyright (C) 2016 Stefano Brambilla
  ======================================================================*/
/*! 
  @file   darcy3dmixed.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
  */

#include <darcy3dmixed.hpp>

//#ifdef USE_SAMG

//SAMG
//#define CSC_INTERFACE 
//#define SPARSE_INTERFACE
//#define CSR_INTERFACE 
// ------------------------------------
// #define DIRECT_SOLVER 
//#define AMG_STAND_ALONE
//#define AMG_ACCELERATED


//#include "samg.h"

/* default 4 Byte integer types */
//#ifndef APPL_INT
//#define APPL_INT int
//#endif
/* end of integer.h */

//#endif

namespace getfem {
    double Lx = 1.0, Ly = 1.0, Lz = 1.0, kappat = 1.0; 
    //! Exact pressure 
    double sol_pt(const bgeot::base_node & x){
        return sin(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*sin(2*pi/Lz*x[2]);
    }
    //! Exact x velocity
    double sol_utx(const bgeot::base_node & x){
        return -2.0*pi*kappat/Lx*cos(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*sin(2.0*pi/Lz*x[2]);
    }
    //! Exact y velocity
    double sol_uty(const bgeot::base_node & x){
        return -2.0*pi*kappat/Ly*sin(2.0*pi/Lx*x[0])*cos(2.0*pi/Ly*x[1])*sin(2.0*pi/Lz*x[2]);
    }
    //! Exact z velocity
    double sol_utz(const bgeot::base_node & x){
        return -2.0*pi*kappat/Lz*sin(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*cos(2.0*pi/Lz*x[2]);
    }
    //! Exact velocity magnitude
    double sol_utm(const bgeot::base_node & x){
        return sqrt(sol_utx(x)*sol_utx(x)+sol_uty(x)*sol_uty(x)+sol_utz(x)*sol_utz(x));
    }
    //! Exact vectorial velocity
    std::vector<double> sol_ut(const bgeot::base_node & x){
        std::vector<double> out(3);
        out[0] = sol_utx(x); out[1] = sol_uty(x); out[2] = sol_utz(x);
        return out;
    }
    //! Exact rhs
    double sol_gt(const bgeot::base_node & x){
        return 4.0*pi*pi*kappat*(1.0/(Lx*Lx)+1.0/(Ly*Ly)+1.0/(Lz*Lz))*sin(2.0*pi/Lx*x[0])*sin(2.0*pi/Ly*x[1])*sin(2.0*pi/Lz*x[2]);
    }
    

	void darcy3dmixed::init (int argc, char *argv[]) 
	{
		std::cout << "initialize the darcy problem..."<<std::endl;
    	//1. Read the .param_darcy filename from standard input
        PARAM.read_command_line(argc, argv);
        //2. Import data (algorithm specifications, boundary conditions, ...)
        import_data();
        //3. Import mesh for tissue (3D) and vessel network (1D)
        build_mesh();
        //4. Set finite elements and integration methods
        set_im_and_fem();
        //5. Build problem parameters
        build_param();
        //6. Build the list of tissue boundary data
        build_tissue_boundary();
	} // end of init


	// Aux methods for init

	//! Import algorithm specifications
	void darcy3dmixed::import_data (void)
	{
		std::cout<<"init part 1: import data!......" <<std::endl;
        #ifdef M3D1D_VERBOSE_
            cout << "Importing descriptors for tissue and vessel problems ..." << endl;
        #endif
		descr_darcy.import(PARAM);
        #ifdef M3D1D_VERBOSE_
            cout << descr_darcy;
        #endif
	};

	//! Import mesh for tissue (3D) and vessel (1D)  
	void darcy3dmixed::build_mesh (void)
    {
		mesht.clear();
        //TODO import 3d mesh 
        /*
		bool test = 0;
		test = PARAM.int_value("TEST_GEOMETRY");
		if(test==0){
        #ifdef M3D1D_VERBOSE_
                    cout << "Importing the 3D mesh for the tissue ...  "   << endl;
        #endif
			import_msh_file(descr.MESH_FILET, mesht);
		}else{*/
        #ifdef M3D1D_VERBOSE_
			cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
        #endif
        string st("GT='" + PARAM.string_value("GT_T") + "'; " +
                "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
                "ORG=" + PARAM.string_value("ORG_T") + "; " +  
                "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
                "NOISED=" + PARAM.string_value("NOISED_T")); 
        cout << "mesht description: " << st << endl;
            regular_mesh(mesht, st);
		//}
	}
	
    //! Set finite elements methods and integration methods 
	void darcy3dmixed::set_im_and_fem (void)
	{
		std::cout<<"init part 2: set fem methods!......" <<std::endl;
        
          #ifdef M3D1D_VERBOSE_
        cout << "Setting IMs for tissue and vessel problems ..." << endl;
        #endif
        pintegration_method pim_t = int_method_descriptor(descr_darcy.IM_TYPET);
        mimt.set_integration_method(mesht.convex_index(), pim_t);
        
        #ifdef M3D1D_VERBOSE_
        cout << "Setting FEMs for tissue and vessel problems ..." << endl;
        #endif
        bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr_darcy.MESH_TYPET);
        pfem pf_Ut = fem_descriptor(descr_darcy.FEM_TYPET_U);
        pfem pf_Pt = fem_descriptor(descr_darcy.FEM_TYPET_P);
        pfem pf_coeft = fem_descriptor(descr_darcy.FEM_TYPET_DATA);
        DIMT = pgt_t->dim();	//DIMV = 1;
        mf_Ut.set_qdim(bgeot::dim_type(DIMT));         
        mf_Ut.set_finite_element(mesht.convex_index(), pf_Ut);
        std::cout<< mf_Ut.get_qdim() << " " << mf_Ut.fem_of_element(0)->target_dim() << std::endl;
        GMM_ASSERT1(mf_Ut.get_qdim() == mf_Ut.fem_of_element(0)->target_dim(), 
            "Intrinsic vectorial FEM used"); // RT0 IS INTRINSIC VECTORIAL!!!
        mf_Pt.set_finite_element(mesht.convex_index(), pf_Pt);
        mf_coeft.set_finite_element(mesht.convex_index(), pf_coeft); 
        #ifdef M3D1D_VERBOSE_
        cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
        #endif
        dof_darcy.set(mf_Ut, mf_Pt, mf_coeft);
        #ifdef M3D1D_VERBOSE_
        cout << std::scientific << dof_darcy;
        #endif

    }


    //! Build problem parameters
    void darcy3dmixed::build_param (void)
    {
        std::cout<<"init part 3: build dimensionless parameters!" <<std::endl;

        #ifdef M3D1D_VERBOSE_
            cout << "Building parameters  ..." << endl;
        #endif
        param_darcy.build(PARAM, mf_coeft);
        #ifdef M3D1D_VERBOSE_
        //cout << param_darcy ;
        #endif
    }
    
    void darcy3dmixed::build_tissue_boundary (void) 
    {
        #ifdef M3D1D_VERBOSE_
            cout << "Building tissue boundary ..." << endl;
        #endif
        BC_darcy.clear();
        BC_darcy.reserve(2*DIMT);
        // Parse BC data
        string label_in = PARAM.string_value("BClabel", "Array of tissue boundary labels");
        string value_in = PARAM.string_value("BCvalue", "Array of tissue boundary values");
        vector<string> labels = split(label_in, ' ');
        vector<string> values = split(value_in, ' ');
        GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
        GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
        for (unsigned f=0; f<2*DIMT; ++f) {
            BC_darcy.emplace_back(labels[f], std::stof(values[f]), 0, f);
            #ifdef M3D1D_VERBOSE_
                cout << "  face " << f << " : " << BC_darcy.back() << endl;
            #endif
        } 

        for (size_type bc=0; bc < BC_darcy.size(); bc++)
            cout<<BC_darcy[bc]<<endl;

        // Build mesht regions
        mesh_region border_faces;
        outer_faces_of_mesh(mesht, border_faces);

        for (mr_visitor i(border_faces); !i.finished(); ++i) {

            assert(i.is_face());

            // Unit outward normal : used to identify faces
            base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
            un /= gmm::vect_norm2(un);

            if (gmm::abs(un[0] + 1.0) < 1.0E-7)      // back
                mesht.region(0).add(i.cv(), i.f());
            else if (gmm::abs(un[0] - 1.0) < 1.0E-7) // front
                mesht.region(1).add(i.cv(), i.f());
            else if (gmm::abs(un[1] + 1.0) < 1.0E-7) // left
                mesht.region(2).add(i.cv(), i.f());
            else if (gmm::abs(un[1] - 1.0) < 1.0E-7) // right
                mesht.region(3).add(i.cv(), i.f());
            else if (gmm::abs(un[2] + 1.0) < 1.0E-7) // bottom
                mesht.region(4).add(i.cv(), i.f());
            else if (gmm::abs(un[2] - 1.0) < 1.0E-7) // top
                mesht.region(5).add(i.cv(), i.f());

        }

    }

    
    void darcy3dmixed::assembly (void)
    {
        #ifdef M3D1D_VERBOSE_
            cout << "Allocating AM_darcy, UM_darcy, FM_darcy ..." << endl;
        #endif
        gmm::resize(AM_darcy, dof_darcy.tot(), dof_darcy.tot()); gmm::clear(AM_darcy);
        gmm::resize(UM_darcy, dof_darcy.tot()); gmm::clear(UM_darcy);
        gmm::resize(FM_darcy, dof_darcy.tot()); gmm::clear(FM_darcy);
    
        #ifdef M3D1D_VERBOSE_
            cout << "Assembling the monolithic matrix AM_darcy ..." << endl;
        #endif
        
        //1 Build the monolithic matrix AM_darcy
        darcy3dmixed::assembly3d_darcy();
        //2 Add bc
        sparse_matrix_type M(dof_darcy.Ut(), dof_darcy.Ut()); gmm::clear(M);
        vector_type F(dof_darcy.Ut()); gmm::clear(F);
        darcy3dmixed::asm_tissue_bc(M, F, mimt, mf_Ut, mf_coeft, BC_darcy);
        gmm::add(M, gmm::sub_matrix(AM_darcy, gmm::sub_interval(0, dof_darcy.Ut()), gmm::sub_interval(0, dof_darcy.Ut())));
        gmm::add(F, gmm::sub_vector(FM_darcy, gmm::sub_interval(0, dof_darcy.Ut())));
        //2 Build the monolithic rhs FM_darcy
        assembly_rhs();
    } // end of assembly
    
        
    void darcy3dmixed::assembly3d_darcy(void)
    {            
        #ifdef M3D1D_VERBOSE    
            cout << "Assembling the 3d darcy (k_p grad, grad) " << endl;
        #endif
            
        // Assembling the mass matrix for velocity (u,u)
        sparse_matrix_type Mtt(dof_darcy.Ut(), dof_darcy.Ut()); gmm::clear(Mtt);
        getfem::asm_mass_matrix(Mtt, mimt, mf_Ut);
        // TODO generalize to vector Kp
        gmm::scale(Mtt, 1.0/param_darcy.Kp()[1]);
        gmm::add(Mtt, gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Ut()), 
                        gmm::sub_interval(0, dof_darcy.Ut())));
        
        // Assembling the mixed term -(p, \div u)
        sparse_matrix_type Dtt(dof_darcy.Pt(), dof_darcy.Ut()); gmm::clear(Dtt);
        getfem::asm_stokes_B (Dtt, mimt, mf_Ut, mf_Pt);
        gmm::add(gmm::transposed(Dtt), gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Ut()), 
                        gmm::sub_interval(dof_darcy.Ut(), dof_darcy.Pt())));
        
        gmm::add(gmm::scaled(Dtt, -1.0), gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(dof_darcy.Ut(), dof_darcy.Pt()), 
                        gmm::sub_interval(0, dof_darcy.Ut())));

        
        gmm::clear(Mtt);
        gmm::clear(Dtt); 
    }
    
    
    void darcy3dmixed::assembly_rhs(void)
    {   
        std::cout << "Adding Rhs..." << std::endl;
        vector_type sol_Gt(mf_coeft.nb_dof()); gmm::clear(sol_Gt);
        interpolation_function(mf_coeft, sol_Gt, sol_gt);
        vector_type Gt(dof_darcy.Pt()); gmm::clear(Gt);
        asm_source_term(Gt, mimt, mf_Pt, mf_coeft, sol_Gt);
        gmm::add(Gt, gmm::sub_vector(FM_darcy,gmm::sub_interval(dof_darcy.Ut(), dof_darcy.Pt())));
        gmm::clear(Gt);

            }
    
    template<typename MAT, typename VEC>
    void
    darcy3dmixed::asm_tissue_bc
        (MAT & M, VEC & F,
        const mesh_im & mim,
        const mesh_fem & mf_u,
        const mesh_fem & mf_data,
        const std::vector<getfem::node> & BC
        )
    {
        GMM_ASSERT1(mf_u.get_qdim()>1,  "invalid data mesh fem (Qdim>1 required)");
        GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");

        std::vector<scalar_type> G(mf_data.nb_dof());
        std::vector<scalar_type> ones(mf_data.nb_dof(), 1.0);
        

        // Define assembly for velocity bc (\Gamma_u)
        /*generic_assembly 
        assemU("g=data$1(#2);" "c=data$2(#2);" 
            "V$1(#1)+=-g(i).comp(Base(#2).vBase(#1).Normal())(i,:,k,k);"
            "M$1(#1,#1)+=c(i).comp(Base(#2).vBase(#1).Normal().vBase(#1).Normal())(i,:,j,j,:,k,k);");
        assemU.push_mi(mim);
        assemU.push_mf(mf_u);
        assemU.push_mf(mf_data);
            assemU.push_data(P0_face);
            assemU.push_data(beta_face);
        assemU.push_vec(F);
        assemU.push_mat(M);*/
        // Define assembly for pressure bc (\Gamma_p)
        generic_assembly 
        assemP("p=data$1(#2);" 
            "V$1(#1)+=-p(i).comp(Base(#2).vBase(#1).Normal())(i,:,k,k);");
        assemP.push_mi(mim);
        assemP.push_mf(mf_u);
        assemP.push_mf(mf_data);
        assemP.push_data(G);
        assemP.push_vec(F);

        for (size_type f=0; f < BC.size(); ++f) {

            GMM_ASSERT1(mf_u.linked_mesh().has_region(f), 
                    "missed mesh region" << f);
            if (BC[f].label=="DIR") { // Dirichlet BC
                gmm::copy(gmm::scaled(ones, BC[f].value), G);	
                assemP.assembly(mf_u.linked_mesh().region(BC[f].rg));
            } /*
            else if (BC[f].label=="MIX") { // Robin BC
                        //Luca MIT
                            std::vector<scalar_type> onesMIX(P0.size(), 1.0);
                            gmm::copy(gmm::scaled(onesMIX, BC[f].value), P0_face);
                            gmm::copy(gmm::scaled(onesMIX, 1/coef[f]), beta_face);
                assemU.assembly(mf_u.linked_mesh().region(BC[f].rg));
            }*/
            else if (BC[f].label=="INT") { // Internal Node
                GMM_WARNING1("internal node passed as boundary.");
            }
            else if (BC[f].label=="JUN") { // Junction Node
                GMM_WARNING1("junction node passed as boundary.");
            }
            else {
                GMM_ASSERT1(0, "Unknown Boundary Condition " << BC[f].label << endl);
            }
        }

    } /* end of asm_tissue_bc */

    
    bool darcy3dmixed::solve (void)
    {
        std::cout << "solve darcy problem" << std::endl;
        
        double time = gmm::uclock_sec();

        if ( descr_darcy.SOLVE_METHOD == "SuperLU" ) { 
            // direct solver //
            #ifdef M3D1D_VERBOSE_
                cout << "  Applying the SuperLU method ... " << endl;
            #endif	
            scalar_type cond;
            gmm::SuperLU_solve(AM_darcy, UM_darcy, FM_darcy, cond);
            cout << "  Condition number (transport problem): " << cond << endl;
        }
        else { // Iterative solver //
            
			gmm::iteration iter(descr_darcy.RES);  // iteration object with the max residu
			iter.set_noisy(1);               // output of iterations (2: sub-iteration)
			iter.set_maxiter(descr_darcy.MAXITER); // maximum number of iterations

			// Preconditioners
			//! TODO Add preconditioner choice to param file
			// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
			gmm::identity_matrix PM; // no precond
			//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
			//gmm::ilu_precond<sparse_matrix_type> PM(AM);
			// ...
			//gmm::clear(AM);
			// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15

			double time_iter = gmm::uclock_sec();
			if ( descr_darcy.SOLVE_METHOD == "CG" ) {
                #ifdef M3D1D_VERBOSE_
                    cout << "  Applying the Conjugate Gradient method ... " << endl;
                #endif
				gmm::identity_matrix PS;  // optional scalar product
				gmm::cg(AM_darcy, UM_darcy, FM_darcy, PS, PM, iter);
			}
			else if ( descr_darcy.SOLVE_METHOD == "BiCGstab" ) {
                #ifdef M3D1D_VERBOSE_
                    cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
                #endif
				gmm::bicgstab(AM_darcy, UM_darcy, FM_darcy, PM, iter);
			}
			else if ( descr_darcy.SOLVE_METHOD == "GMRES" ) {
                #ifdef M3D1D_VERBOSE_
                    cout << "  Applying the Generalized Minimum Residual method ... " << endl;
                #endif
				size_type restart = 50;
                darcy_precond<sparse_matrix_type> P (AM_darcy, mf_Ut, mf_Pt, mimt);
				gmm::gmres(AM_darcy, UM_darcy, FM_darcy, P, restart, iter);
			}
			else if ( descr_darcy.SOLVE_METHOD == "QMR" ) {
                #ifdef M3D1D_VERBOSE_
                    cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
                #endif
				gmm::qmr(AM_darcy, UM_darcy, FM_darcy, PM, iter);
			}
			else if ( descr_darcy.SOLVE_METHOD == "LSCG" ) {
                #ifdef M3D1D_VERBOSE_
                    cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
                #endif
				gmm::least_squares_cg(AM_darcy, UM_darcy, FM_darcy, iter);
			}
			// Check
			if (iter.converged())
				cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
			else if (iter.get_iteration() == descr_darcy.MAXITER)
				cerr << "  ... reached the maximum number of iterations!" << endl;
		}


        cout << endl<<"... total time to solve : " << gmm::uclock_sec() - time << " seconds\n";
        
        //extract the velocity
        vector_type Ut(dof_darcy.Ut());
        gmm::copy(gmm::sub_vector(UM_darcy, 
                    gmm::sub_interval(0, dof_darcy.Ut())), Ut);
        //exact sol
        vector_type Ut_ex (dof_darcy.Ut()); gmm::clear(Ut_ex);
        interpolation_function(mf_Ut, Ut_ex, sol_ut);
        
        scalar_type res = 0; 
        for (int i=0; i< dof_darcy.Ut(); i++)
            res += std::pow(Ut[i] - Ut_ex[i], 2);
        res = std::sqrt(res);
        
        std::cout << "res" << res << std::endl;
        
        
        return true;
        
	}; // end of solve

    
    void darcy3dmixed::export_vtk ()
    {
        #ifdef M3D1D_VERBOSE_
            cout << "Exporting the solution (vtk format) to " << descr_darcy.OUTPUT << " ..." << endl;
        #endif

        // Array of unknown dof of the interstitial pressure
        vector_type Ut(dof_darcy.Ut()); 
        // Array of unknown dof of the vessel pressure
        vector_type Pt(dof_darcy.Pt()); 

        //Copy solution
        gmm::copy(gmm::sub_vector(UM_darcy, 
                    gmm::sub_interval(0, dof_darcy.Ut())), Ut);
        gmm::copy(gmm::sub_vector(UM_darcy, 
                    gmm::sub_interval(dof_darcy.Ut(), dof_darcy.Pt())), Pt);
        
        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Ut ..." << endl;
        #endif
        pfem pf_Ut = fem_descriptor(descr_darcy.FEM_TYPET_U);
        if(pf_Ut->is_lagrange()==0){ 
            /*
                There is no built-in export for non-lagrangian FEM.
                If this is the case, we need to project before exporting.
            */
            #ifdef M3D1D_VERBOSE_
            cout << "    Projecting Ut on P1 ..." << endl;
            #endif
            mesh_fem mf_P1(mesht);
            mf_P1.set_qdim(bgeot::dim_type(DIMT)); 
            mf_P1.set_classical_finite_element(1);
            sparse_matrix_type M_RT0_P1(mf_P1.nb_dof(), dof_darcy.Ut());
            sparse_matrix_type M_P1_P1(mf_P1.nb_dof(), mf_P1.nb_dof());
            vector_type Ut_P1(mf_P1.nb_dof());
                    asm_mass_matrix(M_RT0_P1, mimt, mf_P1, mf_Ut);
            asm_mass_matrix(M_P1_P1,  mimt, mf_P1, mf_P1);
            
            vector_type Utt(mf_P1.nb_dof());
            gmm::mult(M_RT0_P1, Ut, Utt);
            double cond1;
            gmm::SuperLU_solve(M_P1_P1, Ut_P1, Utt, cond1);

            vtk_export exp1(descr_darcy.OUTPUT+"Ut.vtk");
            exp1.exporting(mf_P1);
            exp1.write_mesh();
            exp1.write_point_data(mf_P1, Ut_P1, "Ut");
            
        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Ut_ex ..." << endl;
        #endif
        vector_type Ut_ex(mf_P1.nb_dof());
        interpolation_function(mf_P1, Ut_ex, sol_ut);
        vtk_export exp_Utex(descr_darcy.OUTPUT+"Utexact.vtk");
        exp_Utex.exporting(mf_P1);
        exp_Utex.write_mesh();
        exp_Utex.write_point_data(mf_P1, Ut_ex, "Utex");
        }	
        else {
            vtk_export exp_Ut(descr_darcy.OUTPUT+"Ut.vtk");
            exp_Ut.exporting(mf_Ut);
            exp_Ut.write_mesh();
            exp_Ut.write_point_data(mf_Ut, Ut, "Ut");	
            #ifdef M3D1D_VERBOSE_
                cout << "  Exporting Ut_ex ..." << endl;
            #endif
            vector_type Ut_ex(dof_darcy.Ut());
            interpolation_function(mf_Ut, Ut_ex, sol_utm);
            vtk_export exp_Utex(descr_darcy.OUTPUT+"Utexact.vtk");
            exp_Utex.exporting(mf_Ut);
            exp_Utex.write_mesh();
            exp_Utex.write_point_data(mf_Ut, Ut_ex, "Utex");
            }
        
        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Pt ..." << endl;
        #endif
        vtk_export exp_Pt(descr_darcy.OUTPUT+"Pt.vtk");
        exp_Pt.exporting(mf_Pt);
        exp_Pt.write_mesh();
        exp_Pt.write_point_data(mf_Pt, Pt, "Pt");
        

        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Pt_ex ..." << endl;
        #endif
        vector_type Pt_ex(dof_darcy.Pt());
        interpolation_function(mf_Pt, Pt_ex, sol_pt);
        vtk_export exp_Ptex(descr_darcy.OUTPUT+"Ptexact.vtk");
        exp_Ptex.exporting(mf_Pt);
        exp_Ptex.write_mesh();
        exp_Ptex.write_point_data(mf_Pt, Pt_ex, "Ptex");


        #ifdef M3D1D_VERBOSE_
            cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
        #endif
            
    }; // end of export
    
    /*double U(const bgeot::base_node & x){ return 1.0-x[0];} //source
    
    void darcy3d1d::test ()
    {
        //Solving just the 3D problem with a fixed line source(x)=x-1
        
        //compute source
        std::vector <scalar_type> source (dof_darcy.Pv()); gmm::clear(source);
        interpolation_function(mf_Pv, source, U);
        
        //Assembling the matrix
        gmm::resize(AM_darcy, dof_darcy.Pt(), dof_darcy.Pt()); gmm::clear(AM_darcy);
        gmm::resize(UM_darcy, dof_darcy.Pt()); gmm::clear(UM_darcy);
        gmm::resize(FM_darcy, dof_darcy.Pt()); gmm::clear(FM_darcy);
        
        // 3d darcy block
        darcy3d1d::assembly3d_darcy();
        // 3d1d block
        sparse_matrix_type Btt(dof_darcy.Pt(), dof_darcy.Pt());gmm::clear(Btt);
        sparse_matrix_type Btv(dof_darcy.Pt(), dof_darcy.Pv());gmm::clear(Btv);
        sparse_matrix_type Bvt(dof_darcy.Pv(), dof_darcy.Pt());gmm::clear(Bvt);
        sparse_matrix_type Bvv(dof_darcy.Pv(), dof_darcy.Pv());gmm::clear(Bvv);

        sparse_matrix_type Mbar(dof_darcy.Pv(), dof_darcy.Pt());gmm::clear(Mbar);
        sparse_matrix_type Mlin(dof_darcy.Pv(), dof_darcy.Pt());gmm::clear(Mlin);
        
        asm_exchange_aux_mat(Mbar, Mlin, 
                mimv, mf_Pt, mf_Pv, param_darcy.R(), descr_darcy.NInt);
    
        asm_exchange_mat(Btt, Btv, Bvt, Bvv,
                mimv, mf_Pv, mf_coefv, param_darcy.kappa(), param_darcy.R(), Mbar);
        
        gmm::add(gmm::scaled(Btt, 2.0*pi*param_darcy.R(0)),			 
                    gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Pt()), 
                        gmm::sub_interval(0, dof_darcy.Pt()))); 
        
        gmm::mult(gmm::scaled(Btv, 2.0*pi*param_darcy.R(0)), source, FM_darcy);
        
        // Add boundary conditions
        scalar_type beta_tissue = PARAM.real_value("BETA_T", "Coefficient for MIX condition");
        darcy3d1d::asm_tissue_bc (FM_darcy, AM_darcy, mimt, mf_Pt, mf_coeft, BCt_darcy, beta_tissue);
        
        // Solve the problem
        darcy3d1d::solve();
         
        // Export
        vtk_export exp_Pt(descr_darcy.OUTPUT+"Pt.vtk");
        exp_Pt.exporting(mf_Pt);
        exp_Pt.write_mesh();
        exp_Pt.write_point_data(mf_Pt, UM_darcy, "Pt");        
        
        vtk_export exp_source(descr_darcy.OUTPUT+"source.vtk");
        exp_source.exporting(mf_Pv);
        exp_source.write_mesh();
        exp_source.write_point_data(mf_Pv, source, "Pv"); 
    }*/

} // end of namespace
