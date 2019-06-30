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
            cout << "Setting IMs ..." << endl;
        #endif
        pintegration_method pim_t = int_method_descriptor(descr_darcy.IM_TYPET);
        mimt.set_integration_method(mesht.convex_index(), pim_t);
        
        #ifdef M3D1D_VERBOSE_
            cout << "Setting FEMs ..." << endl;
        #endif
        bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr_darcy.MESH_TYPET);
        DIMT = pgt_t->dim();
        pfem pf_Ut = fem_descriptor(descr_darcy.FEM_TYPET_U);
		pfem pf_Pt = fem_descriptor(descr_darcy.FEM_TYPET_P);
        mf_Ut.set_qdim(bgeot::dim_type(DIMT)); 
        mf_Ut.set_finite_element(mesht.convex_index(), pf_Ut);
		mf_Pt.set_finite_element(mesht.convex_index(), pf_Pt);
		pfem pf_coeft = fem_descriptor(descr_darcy.FEM_TYPET_DATA);
		mf_coeft.set_finite_element(mesht.convex_index(), pf_coeft);
		
        #ifdef M3D1D_VERBOSE_
            cout << "Setting FEM dimensions..." << endl;
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
        //assembly_rhs();
    } // end of assembly
    
        
    void darcy3dmixed::assembly3d_darcy(void)
    {            
        #ifdef M3D1D_VERBOSE    
            cout << "Assembling the 3d darcy (k_p grad, grad) " << endl;
        #endif
            
        // Assembling the mass matrix for velocity (u,u)
        std::cout << "mass matrix" << endl;
        sparse_matrix_type Mtt(dof_darcy.Ut(), dof_darcy.Ut()); gmm::clear(Mtt);
        getfem::asm_mass_matrix(Mtt, mimt, mf_Ut);
        gmm::scale(Mtt, 1.0/param_darcy.Kp()[1]);
        std::cout << "scaled matrix" << endl;
        gmm::add(Mtt, gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Ut()), 
                        gmm::sub_interval(0, dof_darcy.Ut())));
        
        // Assembling the mixed term -(p, \div u)
        std::cout << "mixed term" << std::endl;
        sparse_matrix_type Dtt(dof_darcy.Pt(), dof_darcy.Ut()); gmm::clear(Dtt);
        generic_assembly 
        assem("M$1(#2,#1)+=comp(Base(#2).vGrad(#1))(:,:,i,i);");
        assem.push_mi(mimt);
        assem.push_mf(mf_Ut);
        assem.push_mf(mf_Pt);
        assem.push_mat(Dtt);
        assem.assembly(mesh_region::all_convexes());
               std::cout << "mixed term assembled" << std::endl;
               
        gmm::add(gmm::transposed(gmm::scaled(Dtt, -1.0)), gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Ut()), 
                        gmm::sub_interval(0, dof_darcy.Pt())));
        
        gmm::add(Dtt, gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(dof_darcy.Ut(), dof_darcy.Pt()), 
                        gmm::sub_interval(0, dof_darcy.Pt())));
        
        gmm::clear(Mtt);
        gmm::clear(Dtt); 
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
				gmm::gmres(AM_darcy, UM_darcy, FM_darcy, PM, restart, iter);
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
        
        //export solution
        std::cout << "solved! going to export..." << std::endl; 

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
        vtk_export exp_Ut(descr_darcy.OUTPUT+"Ut.vtk");
        exp_Ut.exporting(mf_Ut);
        exp_Ut.write_mesh();
        exp_Ut.write_point_data(mf_Ut, Ut, "Ut");
        
        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Pt ..." << endl;
        #endif
        vtk_export exp_Pt(descr_darcy.OUTPUT+"Pt.vtk");
        exp_Pt.exporting(mf_Pt);
        exp_Pt.write_mesh();
        exp_Pt.write_point_data(mf_Pt, Pt, "Pt");


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
