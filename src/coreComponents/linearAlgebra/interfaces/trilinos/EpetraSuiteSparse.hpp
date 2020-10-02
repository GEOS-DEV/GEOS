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
 * @file EpetraSuiteSparse.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUITESPARSE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUITESPARSE_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/interfaces/direct/SuiteSparse.hpp"

#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>

namespace geosx
{

/**
 * @brief Converts a matrix from Epetra to SuiteSparse format
 * @param[in] matrix the EpetraMatrix object
 * @param[out] SSData the structure containing the matrix in SuiteSparse format
 */
void ConvertEpetraToSuiteSparseMatrix( EpetraMatrix const & matrix,
                                       SuiteSparseData & SSData,
                                       Epetra_Map * SerialMap,
                                       Epetra_Import * ImportToSerial );

/**
 * @brief Solves a linear system with SuiteSparse (matrix has already been factorized)
 * @param[in,out] SSData the structure containing the matrix in SuiteSparse format
 * @param[in] b the right-hand side in Epetra format
 * @param[out] x the solution in Epetra format
 * @param[out] time time spent in the solution phase
 * @return info error code
 */
int SuiteSparseSolve( SuiteSparseData & SSData,
                      Epetra_Map const * SerialMap,
                      Epetra_Import const * ImportToSerial,
                      EpetraVector const & b,
                      EpetraVector & x,
                      real64 & time );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUITESPARSE_HPP_*/
