/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
  "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
  Course on Advanced Programming for Scientific Computing
  Politecnico di Milano
  A.Y. 2016-2017


  Copyright (C) 2016 Stefano Brambilla
  ======================================================================*/
/*! 
  @file   darcy3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
  */

#include <darcy3d1d.hpp>

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

	void darcy3d1d::init (int argc, char *argv[]) 
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
        //7. Build the list of tissue boundary (and junction) data
        build_vessel_boundary();
	} // end of init


	// Aux methods for init

	//! Import algorithm specifications
	void darcy3d1d::import_data (void)
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
	void darcy3d1d::build_mesh (void)
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
        
        //TODO generate a mesh 1D
            
        #ifdef M3D1D_VERBOSE_
            cout << "Importing the 1D mesh for the vessel (darcy problem)... "   << endl; 
        #endif
		std::ifstream ifs(descr_darcy.MESH_FILEV);
		meshv.clear();
		GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_darcy.MESH_FILEV);
		import_pts_file(ifs, meshv, BCv_darcy, nb_vertices, descr_darcy.MESH_TYPEV);
        nb_branches = nb_vertices.size();
		ifs.close();
	}
	
    //! Set finite elements methods and integration methods 
	void darcy3d1d::set_im_and_fem (void)
	{
		std::cout<<"init part 2: set fem methods!......" <<std::endl;
        
        #ifdef M3D1D_VERBOSE_
            cout << "Setting IMs for tissue and vessel problems ..." << endl;
        #endif
        pintegration_method pim_t = int_method_descriptor(descr_darcy.IM_TYPET);
        pintegration_method pim_v = int_method_descriptor(descr_darcy.IM_TYPEV);
        mimt.set_integration_method(mesht.convex_index(), pim_t);
        mimv.set_integration_method(meshv.convex_index(), pim_v);

        #ifdef M3D1D_VERBOSE_
            cout << "Setting FEMs for tissue and vessel problems ..." << endl;
        #endif
        bgeot::pgeometric_trans pgt_t = bgeot::geometric_trans_descriptor(descr_darcy.MESH_TYPET);
        DIMT = pgt_t->dim();
		pfem pf_Pt = fem_descriptor(descr_darcy.FEM_TYPET_P);
		pfem pf_Pv = fem_descriptor(descr_darcy.FEM_TYPEV_P);
		mf_Pt.set_finite_element(mesht.convex_index(), pf_Pt);
		mf_Pv.set_finite_element(meshv.convex_index(), pf_Pv);
		pfem pf_coeft = fem_descriptor(descr_darcy.FEM_TYPET_DATA);
		pfem pf_coefv = fem_descriptor(descr_darcy.FEM_TYPEV_DATA);
		mf_coeft.set_finite_element(mesht.convex_index(), pf_coeft);
		mf_coefv.set_finite_element(meshv.convex_index(), pf_coefv);

        #ifdef M3D1D_VERBOSE_
            cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
        #endif
		dof_darcy.set(mf_Pt, mf_Pv, mf_coeft, mf_coefv);
        #ifdef M3D1D_VERBOSE_
		cout << std::scientific << dof_darcy;
        #endif
    }


    //! Build problem parameters
    void darcy3d1d::build_param (void)
    {
        std::cout<<"init part 3: build dimensionless parameters!" <<std::endl;

        #ifdef M3D1D_VERBOSE_
            cout << "Building parameters for tissue and vessel problems ..." << endl;
        #endif
        param_darcy.build(PARAM, mf_coeft, mf_coefv);
        #ifdef M3D1D_VERBOSE_
        //cout << param_darcy ;
        #endif
    }

    
    void darcy3d1d::build_tissue_boundary (void) 
    {
        #ifdef M3D1D_VERBOSE_
            cout << "Building tissue boundary ..." << endl;
        #endif
        BCt_darcy.clear();
        BCt_darcy.reserve(2*DIMT);
        // Parse BC data
        string label_in = PARAM.string_value("BClabel", "Array of tissue boundary labels");
        string value_in = PARAM.string_value("BCvalue", "Array of tissue boundary values");
        vector<string> labels = split(label_in, ' ');
        vector<string> values = split(value_in, ' ');
        GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
        GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
        for (unsigned f=0; f<2*DIMT; ++f) {
            BCt_darcy.emplace_back(labels[f], std::stof(values[f]), 0, f);
            #ifdef M3D1D_VERBOSE_
                cout << "  face " << f << " : " << BCt_darcy.back() << endl;
            #endif
        } 

        for (size_type bc=0; bc < BCt_darcy.size(); bc++)
            cout<<BCt_darcy[bc]<<endl;

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
    
    
	void darcy3d1d::build_vessel_boundary (void)
    {
        #ifdef M3D1D_VERBOSE_
            cout << "Building vessel boundary ..." << endl;
        #endif

        dal::bit_vector junctions; // global idx of junctions vertices in meshv
        dal::bit_vector extrema;   // global idx of extreme vertices in meshv

        Jv_darcy.clear();
        nb_extrema=0; 
        nb_junctions=0;

        size_type fer = nb_branches; // first empty region
        GMM_ASSERT1(meshv.has_region(fer)==0, 
                "Overload in meshv region assembling!");
        
        // Visit all the convexes
        for (dal::bv_visitor cv(meshv.convex_index()); !cv.finished(); ++cv){
            
            // Check if cv has exactly 2 faces
            bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
            if (cvs->nb_points()>2) 
                cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
            if (cvs->nb_faces()>2)  
                cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

            // Build regions for BCs and junctions
            // Global idx of mesh vertices
            //   i1|--<-----|i0
            size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
            size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
            
            // Identify vertex type (inflow, junction, outflow)
            if (meshv.convex_to_point(i0).size() == 1){ /* inflow extremum */
                // Update information
                extrema.add(i0);
                nb_extrema++;
                // Build a new region made by a single face
                GMM_ASSERT1(meshv.has_region(fer) == 0, 
                        "Overload in meshv region assembling!");
                meshv.region(fer).add(cv, 1);
                // Store the current index and then update it
                size_type bc = 0; 
                bool found = false;
                while (!found && (bc < BCv_darcy.size())) {
                    found = (i0 == BCv_darcy[bc].idx);
                    if (!found) bc++;
                }
                GMM_ASSERT1(found == true, "Miss a boundary node in BCv list!");
                BCv_darcy[bc].rg = fer; 
                fer++;
                // Store the containing branch index
                size_type branch = 0; 
                bool contained = false;
                while (!contained && branch < nb_branches) {
                    contained = meshv.region(branch).is_in(cv);
                    if (!contained) branch++;
                }
                GMM_ASSERT1(contained == true, "No branch region contains node i0!");
                BCv_darcy[bc].branches.emplace_back(branch); 
            }
            else if (meshv.convex_to_point(i0).size() > 1){ /* junction case */
                
                if (meshv.convex_to_point(i0).size() == 2){
                #ifdef M3D1D_VERBOSE_
                    std::cout << "trivial junction at i0 = " << i0 << std::endl;
                #endif
                }
                // Search for index of containing branch (\mathcal{P}^{in}_j)
                size_type branch = 0; 
                bool contained = false;
                while (!contained && branch < nb_branches ) {
                    contained = meshv.region(branch).is_in(cv);
                    if (!contained) branch++;
                }
                GMM_ASSERT1(contained == true, "No branch region contains node i0!");
                    
                // Check if junction has been already stored, 
                // if not add to the junction list (J) and build a new region
                dal::bit_vector tmp; tmp.add(i0);
                if(!junctions.contains(tmp)){
                    // Store the junction vertex
                    junctions.add(i0);
                    nb_junctions++;
                    GMM_ASSERT1(meshv.has_region(fer) == 0, 
                            "Overload in meshv region assembling!");
                    // Build a new region with idx "first empty region"
                    meshv.region(fer).add(cv, 1); // single-face region
                    // Create a new junction node
                    Jv_darcy.emplace_back("JUN", 0, i0, fer);
                    fer++;
                    // Add the inflow branch
                    Jv_darcy.back().branches.emplace_back(-branch);
                    Jv_darcy.back().value += param_darcy.R(mimv, branch);
                    fer++;
                }
                else{
                // Add the inflow branch (to the right junction node)
                size_type jj = 0;
                bool found = false;
                while (!found && jj < nb_junctions){
                    found = (i0 == Jv_darcy[jj].idx);
                    if (!found) jj++;
                }
                Jv_darcy[jj].value += param_darcy.R(mimv, branch);
                Jv_darcy[jj].branches.emplace_back(-branch);
                }
            }

            if (meshv.convex_to_point(i1).size() == 1){ 
                size_type bc = 0; 
                bool found = false;
                while (!found && (bc < BCv_darcy.size())) {
                    found = (i1 == BCv_darcy[bc].idx);
                    if (!found) bc++;
                }
                if (found){ /* outflow extremum */
                    extrema.add(i1); 
                    nb_extrema++; 
                    // Build a new region made by a single face
                    GMM_ASSERT1(meshv.has_region(fer) == 0, 
                            "Overload in meshv region assembling!");
                    meshv.region(fer).add(cv, 0);
                    // Store the current index and then update it
                    //BCv_darcy[bc].value *= -1.0;
                    BCv_darcy[bc].rg = fer; 
                    fer++;
                    // Store the containing branch index
                    size_type branch = 0; 
                    bool contained = false;
                    while (!contained && branch<nb_branches) {
                        contained = meshv.region(branch).is_in(cv);
                        if (!contained) branch++;
                    }
                    GMM_ASSERT1(contained == true, "No branch region contains node i1!");
                    BCv_darcy[bc].branches.emplace_back(branch); 
                }
                /*else { // interior -> Mixed point
                // "MIX" label via post-processing
                // Build a new region made by a single face
                GMM_ASSERT1(meshv.has_region(fer)==0, 
                "Overload in meshv region assembling!");
                meshv.region(fer).add(cv, 0);
                BCv_darcy.emplace_back("MIX", 0.0, i1, fer);
                fer++;
                // Store the containing branch index
                size_type branch = 0; 
                bool contained = false;
                while (!contained && branch<nb_branches ) {
                contained = meshv.region(branch).is_in(cv);
                if (!contained) branch++;
                }
                GMM_ASSERT1(contained=true, "No branch region contains node i1!");
                BCv_darcy.back().branches.emplace_back(branch); 
                }*/
            }
            else if (meshv.convex_to_point(i1).size() > 1){ /* junction case */
                
                if(meshv.convex_to_point(i1).size() == 2){
                #ifdef M3D1D_VERBOSE_
                    std::cout << "Found a trivial junction at i1 = " << i1 << std::endl;
                #endif
                }
                // Search for index of containing branch (\mathcal{P}^{out}_j)
                size_type branch = 0; 
                bool contained = false;
                while (!contained && branch<nb_branches ) {
                    contained = meshv.region(branch).is_in(cv);
                    if (!contained) branch++;
                }
                GMM_ASSERT1(contained == true, "No branch region contains node i1!");

                // Check if jucntion has been already stored, 
                // if not add to the junction list (J) and build a new region
                dal::bit_vector tmp; tmp.add(i1);
                if(!junctions.contains(tmp)){
                    // Store the junction vertex
                    junctions.add(i1);
                    nb_junctions++;
                    GMM_ASSERT1(meshv.has_region(fer)==0, 
                            "Overload in meshv region assembling!");
                    // Build a new region with idx "first empty region"
                    meshv.region(fer).add(cv, 0);
                    // Create a new junction node
                    Jv_darcy.emplace_back("JUN", 0, i1, fer);
                    // Add the outflow branch
                    Jv_darcy.back().branches.emplace_back(+branch);
                    Jv_darcy.back().value += param_darcy.R(mimv, branch);
                    fer++;
                }
                else {
                    // Add the outflow branch (to the right junction node)
                    size_type jj = 0;
                    bool found = false;
                    while (!found && jj < nb_junctions){
                        found = (i1 == Jv_darcy[jj].idx);
                        if (!found) jj++;
                    }
                    Jv_darcy[jj].branches.emplace_back(+branch);
                    Jv_darcy[jj].value += param_darcy.R(mimv, branch);
                    //cout << "Branch " << branch << " added to junction " << jj << endl;
                }
            }

        } /* end of convexes loop */

        // Ckeck network assembly
        #ifdef M3D1D_VERBOSE_
            cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
            cout << "  Branches:   " << nb_branches << endl
                << "  Vertices:   " << meshv.convex_index().size()+1 << endl;
            cout << "  Extrema:    " << extrema << endl;	  
            for (size_type i=0; i<BCv_darcy.size(); ++i)
                cout << "    -  label=" << BCv_darcy[i].label 
                    << ", value=" << BCv_darcy[i].value << ", ind=" << BCv_darcy[i].idx 
                    << ", rg=" << BCv_darcy[i].rg << ", branches=" << BCv_darcy[i].branches << endl; 
            cout << "  Junctions: " << junctions << endl;
            for (size_type i=0; i<Jv_darcy.size(); ++i)
                cout << "    -  label=" << Jv_darcy[i].label 
                    << ", value=" << Jv_darcy[i].value << ", ind=" << Jv_darcy[i].idx 
                    << ", rg=" << Jv_darcy[i].rg << ", branches=" << Jv_darcy[i].branches << endl; 
            cout << "---------------------------------------- "   << endl;
            
        #endif

    } /* end of build_vessel_boundary */
    
    
    void darcy3d1d::assembly (void)
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
        darcy3d1d::assembly3d_darcy();
        darcy3d1d::assembly1d_darcy();
        darcy3d1d::assembly3d1d_darcy();
        darcy3d1d::assembly_bc();
        //2 Build the monolithic rhs FM_darcy
        //assembly_rhs();
    } // end of assembly
    
        
    void darcy3d1d::assembly3d_darcy(void)
    {            
        // Assembling the 3d darcy (k_p grad, grad)
        #ifdef M3D1D_VERBOSE    
            cout << "Assembling the 3d darcy (k_p grad, grad) " << endl;
        #endif
      
        sparse_matrix_type Dt(dof_darcy.Pt(), dof_darcy.Pt()); gmm::clear(Dt);
        getfem::asm_stiffness_matrix_for_laplacian(Dt, mimt, mf_Pt, mf_coeft, param_darcy.Kp());
        gmm::add(Dt, gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Pt()), 
                        gmm::sub_interval(0, dof_darcy.Pt())));
    
        gmm::MatrixMarket_IO::write("Dt.mm",Dt);

        gmm::clear(Dt); 
    }
    
    void darcy3d1d::assembly1d_darcy (void)
    {
        // Assembling the 1d darcy (k_v ds, ds)
        #ifdef M3D1D_VERBOSE    
            cout << "Assembling the 1d darcy (k_v ds, ds) " << endl;
        #endif
        sparse_matrix_type Dv(dof_darcy.Pv(), dof_darcy.Pv()); gmm::clear(Dv);	
        getfem::asm_stiffness_matrix_for_laplacian(Dv,mimv,mf_Pv, mf_coefv, param_darcy.Kv());	
        gmm::add(Dv, gmm::sub_matrix(AM_darcy,
                        gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()),
                        gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())));
        
        gmm::MatrixMarket_IO::write("Dv.mm",Dv);
        
        gmm::clear(Dv); 
    }
    
    void darcy3d1d::assembly3d1d_darcy (void)
    {
        // Assembling coupling terms
        #ifdef M3D1D_VERBOSE    
            cout << "Assembling coupling terms " << endl;
        #endif
        // Tissue-to-tissue exchange matrix
        sparse_matrix_type Btt(dof_darcy.Pt(), dof_darcy.Pt());gmm::clear(Btt);
        // Vessel-to-tissue exchange matrix
        sparse_matrix_type Btv(dof_darcy.Pt(), dof_darcy.Pv());gmm::clear(Btv);
        // Tissue-to-vessel exchange matrix
        sparse_matrix_type Bvt(dof_darcy.Pv(), dof_darcy.Pt());gmm::clear(Bvt);
        // Vessel-to-vessel exchange matrix
        sparse_matrix_type Bvv(dof_darcy.Pv(), dof_darcy.Pv());gmm::clear(Bvv);

        // Aux tissue-to-vessel averaging matrix
        sparse_matrix_type Mbar(dof_darcy.Pv(), dof_darcy.Pt());gmm::clear(Mbar);
        // Aux tissue-to-vessel interpolation matrix
        sparse_matrix_type Mlin(dof_darcy.Pv(), dof_darcy.Pt());gmm::clear(Mlin);
        
        asm_exchange_aux_mat(Mbar, Mlin, 
                mimv, mf_Pt, mf_Pv, param_darcy.R(), descr_darcy.NInt);

        asm_exchange_mat(Btt, Btv, Bvt, Bvv,
                mimv, mf_Pv, mf_coefv, param_darcy.kappa(), param_darcy.R(), Mbar);
        
        gmm::add(Btt,			 
                    gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Pt()), 
                        gmm::sub_interval(0, dof_darcy.Pt()))); 
        
        gmm::add(gmm::scaled(Btv, -1.0),									
                    gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(0, dof_darcy.Pt()),
                        gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()))); 

        gmm::add(gmm::scaled(Bvt, -1.0),  	
                    gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()),
                        gmm::sub_interval(0, dof_darcy.Pt())));

        gmm::add(Bvv,								
                    gmm::sub_matrix(AM_darcy, 
                        gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()), 
                        gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())));
        //***************************************
        //checkig matrices
        gmm::MatrixMarket_IO::write("Mbar.mm",Mbar);
        gmm::MatrixMarket_IO::write("Bvv.mm",gmm::scaled(Bvv, 2.0/param_darcy.R(0)));
        gmm::MatrixMarket_IO::write("Btt.mm",gmm::scaled(Btt, 2.0*pi*param_darcy.R(0)));
        gmm::MatrixMarket_IO::write("Btv.mm",gmm::scaled(Btv, -2.0*pi*param_darcy.R(0)));
        gmm::MatrixMarket_IO::write("Bvt.mm",gmm::scaled(Bvt, 2.0/param_darcy.R(0)));
        
        // De-allocate memory
        gmm::clear(Mbar);  gmm::clear(Mlin);
        gmm::clear(Btt);  gmm::clear(Btv);
        gmm::clear(Bvt);  gmm::clear(Bvv);
        
    }
    
    
    void darcy3d1d::assembly_bc (void)
    {   
        //Assembling BC for the 3d problem
        vector_type Ft_temp(dof_darcy.Pt()); gmm::clear(Ft_temp);
        gmm::copy(gmm::sub_vector(FM_darcy, gmm::sub_interval(0, dof_darcy.Pt())), Ft_temp);
        gmm::clear(gmm::sub_vector(FM_darcy, gmm::sub_interval(0, dof_darcy.Pt())));
        sparse_matrix_type At_temp(dof_darcy.Pt(),dof_darcy.Pt()); gmm::clear(At_temp);
        gmm::copy(gmm::sub_matrix(AM_darcy, 
                                  gmm::sub_interval(0, dof_darcy.Pt()),
                                  gmm::sub_interval(0, dof_darcy.Pt())
                                 ), At_temp);
        gmm::clear(gmm::sub_matrix(AM_darcy, 
                                  gmm::sub_interval(0, dof_darcy.Pt()),
                                  gmm::sub_interval(0, dof_darcy.Pt())
                                 ));
        
        scalar_type beta_tissue = PARAM.real_value("BETA_T", "Coefficient for MIX condition");
        
        darcy3d1d::asm_tissue_bc
        (Ft_temp, At_temp, mimt, mf_Pt, mf_coeft, BCt_darcy, beta_tissue);
        
        gmm::add(Ft_temp, gmm::sub_vector(FM_darcy, gmm::sub_interval(0, dof_darcy.Pt())));
        gmm::add(At_temp, gmm::sub_matrix(AM_darcy, 
                                  gmm::sub_interval(0, dof_darcy.Pt()),
                                  gmm::sub_interval(0, dof_darcy.Pt())
                                 ));
        
        gmm::clear(At_temp);
        gmm::clear(Ft_temp);
        
        
        //Assembling bc for the 1d problem
        vector_type Fv_temp(dof_darcy.Pv()); gmm::clear(Fv_temp);
        gmm::add(gmm::sub_vector(FM_darcy, gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())), Fv_temp);
        gmm::clear(gmm::sub_vector(FM_darcy, gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())));
        sparse_matrix_type Av_temp(dof_darcy.Pv(),dof_darcy.Pv()); gmm::clear(Av_temp);
        gmm::add(gmm::sub_matrix(AM_darcy, 
                                  gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()),
                                  gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())
                                 ), Av_temp);
        gmm::clear(gmm::sub_matrix(AM_darcy, 
                                  gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()),
                                  gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())
                                 ));
        
        scalar_type beta_vessel = PARAM.real_value("BETA_V", "Coefficient for MIX condition in the vessels");
        
        darcy3d1d::asm_network_bc
        (Fv_temp, Av_temp, mimv, mf_Pv, mf_coefv, BCv_darcy, beta_vessel, param_darcy.R());      
        gmm::add(Fv_temp, gmm::sub_vector(FM_darcy, gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())));
        gmm::add(Av_temp, gmm::sub_matrix(AM_darcy, 
                                  gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv()),
                                  gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())
                                 ));
        
        gmm::clear(Av_temp);
        gmm::clear(Fv_temp);   
        
    }
    
    template<typename MAT, typename VEC>
    void darcy3d1d::asm_tissue_bc
        (VEC & F,
        MAT & M,
        const mesh_im & mim,
        const mesh_fem & mf_p,
        const mesh_fem & mf_data,
        const std::vector<getfem::node> & BC,
        const scalar_type beta
        )
    {

        GMM_ASSERT1(mf_p.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
        GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");

        for (size_type bc=0; bc < BC.size(); ++bc) {
            GMM_ASSERT1(mf_p.linked_mesh().has_region(bc), "missed mesh region" << bc);
            if (BC[bc].label=="DIR") { 
                // Dirichlet BC
                VEC BC_temp(mf_p.nb_dof(), BC[bc].value);
                getfem::assembling_Dirichlet_condition(M, F, mf_p, BC[bc].rg, BC_temp);
                gmm::clear(BC_temp);				
            } 
            else if (BC[bc].label=="MIX") { 
                // Robin BC
                VEC BETA(mf_data.nb_dof(), beta);
                getfem::asm_mass_matrix_param(M, mim, mf_p, mf_data, BETA,mf_p.linked_mesh().region(BC[bc].rg) );
                
                VEC BETA_C0(mf_data.nb_dof(), beta*BC[bc].value);
                asm_source_term(F,mim, mf_p, mf_data,BETA_C0);
            }
            else if (BC[bc].label=="INT") { 
                // Internal Node
                GMM_WARNING1("internal node passed as boundary.");
            }
            else if (BC[bc].label=="JUN") { 
                // Junction Node
                GMM_WARNING1("junction node passed as boundary.");
            }
            else {
                GMM_ASSERT1(0, "Unknown Boundary Condition " << BC[bc].label << endl);
            }
        }

    } /* end of asm_tissue_bc */
    
    template<typename MAT, typename VEC>
    void darcy3d1d::asm_network_bc
        (VEC & F, MAT & M, 
        const mesh_im & mim,
        const mesh_fem & mf_p,
        const mesh_fem & mf_data,
        const std::vector<getfem::node> & BC,
        const scalar_type beta,
        const VEC & R) 
    { 
        GMM_ASSERT1(mf_p.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
        GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");
        
        for (size_type bc=0; bc < BC.size(); bc++) { 
            GMM_ASSERT1(mf_p.linked_mesh().has_region(bc), "missed mesh region" << bc);
            if (BC[bc].label=="DIR") { // Dirichlet BC
                VEC BC_temp(mf_p.nb_dof(), BC[bc].value);
                getfem::assembling_Dirichlet_condition(M, F, mf_p, BC[bc].rg, BC_temp);
                //NOTE Everything works without penalty if we don't scale the vessel eq by pi*R^2
                double penalty =1.e+9;    
                F[mf_p.basic_dof_on_region(BC[bc].rg).first()]= BC[bc].value*penalty;
                M(mf_p.basic_dof_on_region(BC[bc].rg).first(),mf_p.basic_dof_on_region(BC[bc].rg).first())=penalty;
                gmm::clear(BC_temp);			
            } 
            else if (BC[bc].label=="MIX") { // Robin BC
                VEC BETA(mf_data.nb_dof(), beta*pi);
                for (size_type i = 0; i < mf_data.nb_dof(); i++ )
                    BETA[i]= BETA[i]*R[i]*R[i];
                getfem::asm_mass_matrix_param(M, mim, mf_p, mf_data, BETA,mf_p.linked_mesh().region(BC[bc].rg) ); //int(beta*cv*bv)
                VEC BETA_C0(mf_data.nb_dof(), pi*beta*BC[bc].value);
                for (size_type j = 0; j < mf_data.nb_dof(); j++ )
                    BETA_C0[j]= BETA_C0[j]*R[j]*R[j];
                asm_source_term(F,mim, mf_p, mf_data,BETA_C0); //int(beta*c0*bv)
            }
            else if (BC[bc].label=="INT") { // Internal Node
                GMM_WARNING1("internal node passed as boundary.");
            }
            else if (BC[bc].label=="JUN") { // Junction Node
                GMM_WARNING1("junction node passed as boundary.");
            }
            else {
                GMM_ASSERT1(0, "Unknown Boundary Condition"<< BC[bc].label << endl);
            }
        }

    }
    
    bool darcy3d1d::solve_samg(void)
    {   
        
        cout << "Solving with samg.." << endl;
        
        gmm::csr_matrix <scalar_type> AM_csr;
        gmm::copy(AM_darcy, AM_csr);
        /*AMG sys("Sys_samg", AM_csr, UM_darcy, FM_darcy);
        sys.csr2samg();
        sys.solve();
        //solution
        for(int i=0; i< dof_darcy.tot(); i++)
            UM_darcy[i] = sys.getsol()[i];*/
        
        return true;
    }
    
    
    bool darcy3d1d::solve (void)
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

    
    void darcy3d1d::export_vtk ()
    {
        #ifdef M3D1D_VERBOSE_
            cout << "Exporting the solution (vtk format) to " << descr_darcy.OUTPUT << " ..." << endl;
        #endif

        // Array of unknown dof of the interstitial pressure
        vector_type Pt(dof_darcy.Pt()); 
        // Array of unknown dof of the vessel pressure
        vector_type Pv(dof_darcy.Pv()); 

        //Copy solution
        gmm::copy(gmm::sub_vector(UM_darcy, 
                    gmm::sub_interval(0, dof_darcy.Pt())), Pt);
        gmm::copy(gmm::sub_vector(UM_darcy, 
                    gmm::sub_interval(dof_darcy.Pt(), dof_darcy.Pv())), Pv);

        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Pt ..." << endl;
        #endif
        vtk_export exp_Pt(descr_darcy.OUTPUT+"Pt.vtk");
        exp_Pt.exporting(mf_Pt);
        exp_Pt.write_mesh();
        exp_Pt.write_point_data(mf_Pt, Pt, "Pt");

        #ifdef M3D1D_VERBOSE_
            cout << "  Exporting Pv ..." << endl;
        #endif
        vtk_export exp_Pv(descr_darcy.OUTPUT+"Pv.vtk");
        exp_Pv.exporting(mf_Pv);
        exp_Pv.write_mesh();
        exp_Pv.write_point_data(mf_Pv, Pv, "Pv");

        #ifdef M3D1D_VERBOSE_
            cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
        #endif
            
    }; // end of export
    
    double U(const bgeot::base_node & x){ return 1.0-x[0];} //source
    
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
    }

} // end of namespace