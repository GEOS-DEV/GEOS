/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_

#include "common/DataTypes.hpp"

#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

namespace geosx
{

/**
 * @brief Convert GEOSX integer value to hypre int
 * @param index the input value
 * @return the converted value
 */
inline HYPRE_Int toHYPRE_Int( integer const index )
{
  return LvArray::integerConversion< HYPRE_Int >( index );
}

/**
 * @brief Convert GEOSX global index value to hypre bigint
 * @param index the input value
 * @return the converted value
 */
inline HYPRE_BigInt toHYPRE_BigInt( globalIndex const index )
{
  return LvArray::integerConversion< HYPRE_BigInt >( index );
}

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to HYPRE_BigInt
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt * toHYPRE_BigInt( globalIndex * const index )
{
  return reinterpret_cast< HYPRE_BigInt * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to HYPRE_BigInt
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt const * toHYPRE_BigInt( globalIndex const * const index )
{
  return reinterpret_cast< HYPRE_BigInt const * >(index);
}

/**
 * @brief Container for hypre preconditioner function pointers.
 */
struct HyprePrecFuncs
{
  using SetupFunc = HYPRE_Int (*)( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector );
  using ApplyFunc = HYPRE_Int (*)( HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector );
  using DestroyFunc = HYPRE_Int (*)( HYPRE_Solver );

  SetupFunc setup{};
  ApplyFunc apply{};
  DestroyFunc destroy{};
};

struct HypreSolverFuncs
{
  using SetPrecondFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_Solver );
  using SetupFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                     HYPRE_ParCSRMatrix,
                                     HYPRE_ParVector,
                                     HYPRE_ParVector );
  using SolveFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                     HYPRE_ParCSRMatrix,
                                     HYPRE_ParVector,
                                     HYPRE_ParVector );
  using GetNumIter = HYPRE_Int ( * )( HYPRE_Solver solver,
                                      HYPRE_Int * num_iterations );
  using GetFinalNorm = HYPRE_Int ( * )( HYPRE_Solver solver,
                                        HYPRE_Real * norm );
  using DestroyFunc = HYPRE_Int ( * )( HYPRE_Solver );

  SetPrecondFunc setPrecond{};
  SetupFunc setup{};
  SolveFunc solve{};
  GetNumIter getNumIter{};
  GetFinalNorm getFinalNorm{};
  DestroyFunc destroy{};
};

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_*/
