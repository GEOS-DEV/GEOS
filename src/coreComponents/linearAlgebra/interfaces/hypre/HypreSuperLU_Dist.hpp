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
 * @file HypreSuperLU_Dist.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUPERLU_DIST_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUPERLU_DIST_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/direct/SuperLU_Dist.hpp"

#include <seq_mv.h>

namespace geosx
{

/**
 * @brief Converts a matrix from Hypre to SuperLU_Dist format
 * @param[in] matrix the HypreMatrix object
 * @param[out] localMatrix local matrix in Hypre format
 * @param[out] SLUDData the structure containing the matrix in SuperLU_Dist format
 */
void HypreConvertToSuperMatrix( HypreMatrix const & matrix,
                                hypre_CSRMatrix * & localMatrix,
                                SuperLU_Dist & SLUDData );

/**
 * @brief Destroys data needed to convert a matrix from Hypre to SuperLU_Dist format
 * @param[inout] localMatrix local matrix in Hypre format
 */
void HypreDestroyAdditionalData( hypre_CSRMatrix * & localMatrix );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPRESUPERLU_DIST_HPP_*/
