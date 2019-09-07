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
* @file TrilinosSolver.hpp
*
*  Created on: Aug 9, 2018
*      Author: Matthias Cremon
*/

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSSOLVER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSSOLVER_HPP_

#include <AztecOO.h>
#include <Amesos.h>
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "EpetraMatrix.hpp"
#include "EpetraVector.hpp"
#include "LinearSolverParameters.hpp"

namespace geosx
{

/**
* \class TrilinosSolver
* \brief This class creates and provides basic support for AztecOO, Amesos and ML libraries.
*/

class TrilinosSolver
{
public:

/**
* @brief Solver constructor, with parameter list reference
*
*/
TrilinosSolver( LinearSolverParameters const & parameters );

/**
* @brief Virtual destructor.
*
*/
virtual ~TrilinosSolver() = default;

/**
* @brief Solve system with an iterative solver (HARD CODED PARAMETERS, GMRES).
*
* Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
*/

void solve( EpetraMatrix & mat,
EpetraVector & sol,
EpetraVector & rhs );

private:

LinearSolverParameters const & m_parameters;

void solve_direct( EpetraMatrix & mat,
EpetraVector & sol,
EpetraVector & rhs );

void solve_krylov( EpetraMatrix & mat,
EpetraVector & sol,
EpetraVector & rhs );

};

} // end geosx namespace

#endif /* TRILINOSSOLVER_HPP_ */
