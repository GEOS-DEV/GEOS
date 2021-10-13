/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_VEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_VEM_HPP_

#include "physicsSolvers/simplePDE/LaplaceBaseH1.hpp"  // a base class shared by all Laplace solvers

namespace geosx
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

class LaplaceVEM : public LaplaceBaseH1
{
public:
  // The default nullary constructor is disabled to avoid compiler auto-generation:
  LaplaceVEM() = delete;

  // The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  LaplaceVEM( const string & name,
              Group * const parent );

  // Destructor
  virtual ~LaplaceVEM() override;

  // "CatalogName()" return the string used as XML tag in the input file.
  // It ties the XML tag with this C++ classes. This is important.
  static string catalogName() { return "LaplaceVEM"; }

// /**
//  * @defgroup Solver Interface Functions
//  *
//  * These functions provide the primary interface that is required for derived classes
//  */
// /**@{*/

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  /**@}*/

private:

  static constexpr localIndex m_maxCellNodes = 10;
  static constexpr localIndex m_maxFaceNodes = 10;

};
} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_VEM_HPP_ */
