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
* @file HypreInterface.hpp
*/

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_

//#include "HypreSparseMatrix.hpp"
//#include "HypreVector.hpp"
//#include "HypreSolver.hpp"

namespace geosx
{

/**
* \class HypreInterface
* \brief This class holds aliases based on the Hypre library.
*/

class HypreInterface
{
public:

using laiLID = hypreTypes::lid;
using laiGID = hypreTypes::gid;

// Epetra matrix and vector wrappers
using ParallelMatrix = HypreSparseMatrix;
using ParallelVector = HypreVector;

using LinearSolver = HypreSolver;

//! @name Constructor/Destructor Methods
//@{
/**
* @brief Empty constructor.
*/
HypreInterface() = default;
/**
* @brief Destructor.
*/
~HypreInterface() = default;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_HYPREINTERFACE_HPP_ */
