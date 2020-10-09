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
 * @param[out] serialMap Epetra map for the serial matrix
 * @param[out] importToSerial Epetra import to convert from parallel to serial matrix
 */
void ConvertEpetraToSuiteSparseMatrix( EpetraMatrix const & matrix,
                                       SuiteSparse & SSData,
                                       Epetra_Map * & serialMap,
                                       Epetra_Import * & importToSerial );

/**
 * @brief Solves a linear system with SuiteSparse (matrix has already been factorized)
 * @param[in,out] SSData the structure containing the matrix in SuiteSparse format
 * @param[in] serialMap Epetra map for the serial matrix
 * @param[in] importToSerial Epetra import to convert from parallel to serial matrix
 * @param[in] b the right-hand side in Epetra format
 * @param[out] x the solution in Epetra format
 * @param[in] transpose whether to solve for the original or the transpose matrix
 * @return info error code
 */
int SuiteSparseSolve( SuiteSparse & SSData,
                      Epetra_Map const * serialMap,
                      Epetra_Import const * importToSerial,
                      EpetraVector const & b,
                      EpetraVector & x,
                      bool transpose = false );

/**
 * @brief Computes an accurate condition number (time consuming function!!!)
 * @param[in] matrix the EpetraMatrix object
 * @param[in] serialMap Epetra map for the serial matrix
 * @param[in] importToSerial Epetra import to convert from parallel to serial matrix
 * @param[in] SSData the structure containing the matrix in SuiteSparse format
 * @return the condition number
 */
real64 EpetraSuiteSparseCond( EpetraMatrix const & matrix,
                              Epetra_Map const * serialMap,
                              Epetra_Import const * importToSerial,
                              SuiteSparse & SSData );

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRASUITESPARSE_HPP_*/
