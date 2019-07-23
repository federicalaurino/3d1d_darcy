/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
3d1d Darcy problem - IJGE
======================================================================*/
/*! 
  @file   defines.hpp
  @author Federica Laurino <<federica.laurino@polimi.it>
  @date   2019.
  @brief  Miscellaneous definitions for the 3D/1D coupling.
 */
#ifndef M3D1D_DEFINES_HPP_
#define M3D1D_DEFINES_HPP_

#include <gmm/gmm.h>

namespace getfem {

// Some useful abbreviations

// Useful type definitions (built using the predefined types in Gmm++)
//! Special class for small (dim < 16) vectors 
using bgeot::base_small_vector;  
//! Geometrical nodes (derived from base_small_vector)
using bgeot::base_node;   	 	 
//! Double-precision FP numbers 
using bgeot::scalar_type; 	 	
//! Unsigned long integers 
using bgeot::size_type;   	 	 
//! Short integers 
using bgeot::short_type;         
//! Type for vector of integers 
typedef std::vector<size_type> vector_size_type;
//! Type for dense vector
typedef std::vector<scalar_type> vector_type;
//! Type for sparse vector (std::vector) 
typedef gmm::rsvector<scalar_type> sparse_vector_type;
//! Type for dense matrix
typedef gmm::row_matrix<vector_type> matrix_type;
//! Type for sparse matrix
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;    


//! Definition of \pi
const scalar_type pi = std::atan(1)*4; 

} /* end of namespace */

#endif
