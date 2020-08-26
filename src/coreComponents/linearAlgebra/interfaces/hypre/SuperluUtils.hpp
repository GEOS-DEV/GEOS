/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SuperluUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_SUPERLUUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_SUPERLUUTILS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"

#include <superlu_ddefs.h>

namespace geosx
{

/**
 * @brief Convert GEOSX global index value to SLUD int_t
 * @param index the input value
 * @return the converted value
 */
inline int_t toSuperlu_intT( globalIndex const index )
{
  return LvArray::integerConversion< int_t >( index );
}

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to int_t
 * @param[in] index the input array
 * @return the converted array
 */
inline int_t * toSuperlu_intT( globalIndex * const index )
{
  return reinterpret_cast< int_t * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to int_t
 * @param[in] index the input array
 * @return the converted array
 */
inline int_t const * toSuperlu_intT( globalIndex const * const index )
{
  return reinterpret_cast< int_t const * >(index);
}

/**
 * @brief Converts a matrix from Hypre to SuperLU_Dist format
 * @param[in] matrix the HypreMatrix object
 * @param[out] rowPtr the row pointers (SuperMatrix will point to it)
 * @param[out] cols the column indices (SuperMatrix will point to it)
 * @param[out] vals the values (SuperMatrix will point to it)
 * @param[out] SLUDMat the matrix in SuperLU_Dist format
 */
void ConvertToSuperMatrix( HypreMatrix const & matrix,
                           array1d< globalIndex > & rowPtr,
                           array1d< globalIndex > & cols,
                           array1d< real64 > & vals,
                           SuperMatrix & SLUDMat );

/**
 * @brief Solves a linear system with SuperLU_Dist. Time is split among factorization and solution.
 * @param[in] SLUDMat the matrix
 * @param[in] b the right-hand side in Hypre format
 * @param[out] x the solution in Hypre format
 * @param[in] comm the MPI communicator
 * @param[in,out] options SuperLU_Dist option object
 * @param[in] logLevel the output level
 * @param[out] timeFact time spent in the factorization phase
 * @param[out] timeSolve time spent in the solution phase
 * @return info error code
 */
int SolveSuperMatrix( SuperMatrix & SLUDMat,
                      HypreVector const & b,
                      HypreVector & x,
                      MPI_Comm const & comm,
                      superlu_dist_options_t & options,
                      integer const & logLevel,
                      real64 & timeFact,
                      real64 & timeSolve );

/**
 * @brief Deallocates a SuperLU_Dist matrix
 * @param[in,out] SLUDMat the matrix
 */
void DestroySuperMatrix( SuperMatrix & SLUDMat );

/**
 * @brief Converts from GEOSX to SuperLU_Dist columns permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
colperm_t const & getColPermType( string const & value );

/**
 * @brief Converts from GEOSX to SuperLU_Dist rows permutation option
 * @param[in] value the GEOSX option
 * @return the SuperLU_Dist option
 */
rowperm_t const & getRowPermType( string const & value );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_SUPERLUUTILS_HPP_*/
