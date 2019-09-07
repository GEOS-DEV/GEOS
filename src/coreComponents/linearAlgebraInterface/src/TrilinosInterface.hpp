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
* @file TrilinosInterface.hpp
*/

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_

#include "TrilinosSolver.hpp"
#include "EpetraMatrix.hpp"
#include "EpetraVector.hpp"

namespace geosx
{

/**
* \class TrilinosInterface
* \brief This class holds aliases based on the Trilinos library.
*/

class TrilinosInterface
{
public:

//using lid = trilinosTypes::lid; // no longer necessary
//using gid = trilinosTypes::gid; // no longer necessary

// Epetra matrix and vector wrappers
using ParallelMatrix = EpetraMatrix;
using ParallelVector = EpetraVector;
using LinearSolver   = TrilinosSolver;

//! @name Constructor/Destructor Methods
//@{
/**
* @brief Empty constructor.
*/
TrilinosInterface() = default;

/**
* @brief Destructor.
*
*/
~TrilinosInterface() = default;
//@}

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_TRILINOSINTERFACE_HPP_ */
