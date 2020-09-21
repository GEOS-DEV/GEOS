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
 * @file HypreSuperlu.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUPERLU_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUPERLU_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <HYPRE_utilities.h>
#include <seq_mv.h>
#include <superlu_ddefs.h>

namespace geosx
{

/**
 * SuperLU_Dist data
 */
struct SuperLU_DistData
{
  hypre_CSRMatrix * localStrip;       //!< local strip of the HypreMatrix
  int_t * rowPtr;                     //!< row pointer
  SuperMatrix mat;                    //!< SuperLU_Dist matrix format
  dScalePermstruct_t ScalePermstruct; //!< data structure to scale and permute the matrix
  dLUstruct_t LUstruct;               //!< data structure to store the LU factorization
  SuperLUStat_t stat;                 //!< data structure to gather some statistics
  gridinfo_t grid;                    //!< SuperLU_Dist MPI subdivision of load
  dSOLVEstruct_t SOLVEstruct;         //!< data structure to solve the matrix
  superlu_dist_options_t options;     //!< SuperLU_Dist options
  MPI_Comm comm;                      //!< MPI communicator
};

/**
 * @brief Creates the SuperLU_Dist data structure
 * @param[in] matrix the HypreMatrix object
 * @param[in] params the linear solver parameters
 * @param[out] SLUDData the structure containing the matrix in SuperLU_Dist format
 */
void SuperLU_DistCreate( HypreMatrix const & matrix,
                         LinearSolverParameters const & params,
                         SuperLU_DistData & SLUDData );

/**
 * @brief Factorizes a linear system with SuperLU_Dist
 * @param[in,out] SLUDData the structure containing the matrix in SuperLU_Dist format
 * @param[out] time time spent in the factorization phase
 * @return info error code
 */
int SuperLU_DistSetup( SuperLU_DistData & SLUDData,
                       real64 & time );

/**
 * @brief Solves a linear system with SuperLU_Dist (matrix has already been factorized)
 * @param[in,out] SLUDData the structure containing the matrix in SuperLU_Dist format
 * @param[in] b the right-hand side in Hypre format
 * @param[out] x the solution in Hypre format
 * @param[out] time time spent in the solution phase
 * @return info error code
 */
int SuperLU_DistSolve( SuperLU_DistData & SLUDData,
                       HypreVector const & b,
                       HypreVector & x,
                       real64 & time );

/**
 * @brief Estimates the condition number of the matrix
 * @param[in] SLUDData the structure containing the matrix in SuperLU_Dist format
 * @return the estimated condition number
 */
real64 SuperLU_DistCondEst( SuperLU_DistData & SLUDData );

/**
 * @brief Deallocates a SuperLU_Dist data structure
 * @param[in,out] SLUDData the structure containing the matrix in SuperLU_Dist format
 */
void SuperLU_DistDestroy( SuperLU_DistData & SLUDData );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUPERLU_HPP_*/
