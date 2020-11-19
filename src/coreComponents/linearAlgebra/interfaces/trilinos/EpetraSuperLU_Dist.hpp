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
 * @file EpetraSuperLU_Dist.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUPERLU_DIST_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUPERLU_DIST_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/interfaces/direct/SuperLU_Dist.hpp"

//#include <seq_mv.h>

namespace geosx
{

/**
 * @brief Converts a matrix from Epetra to SuperLU_Dist format
 * @param[in] matrix the EpetraMatrix object
 * @param[out] SLUDData the structure containing the matrix in SuperLU_Dist format
 */
void EpetraConvertToSuperMatrix( EpetraMatrix const & matrix,
                                 SuperLU_Dist & SLUDData );

/**
 * @brief Computes an accurate condition number (time consuming function!!!)
 * @param[in] matrix the EpetraMatrix object
 * @param[in] SLUDData the structure containing the matrix in SuperLU_Dist format
 * @return the condition number
 */
real64 EpetraSuperLU_DistCond( EpetraMatrix const & matrix, SuperLU_Dist & SLUDData );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUPERLU_DIST_HPP_*/
