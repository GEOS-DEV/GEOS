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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_DIETERICH_SEISMICITY_RATE_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_DIETERICH_SEISMICITY_RATE_HPP_

#include "physicsSolvers/inducedSeismicity/SeismicityRateBase.hpp"  // a base class shared by all Laplace solvers

namespace geos
{

// Like most physics solvers, the Laplace solver derives from a generic SolverBase class.
// The base class is densely Doxygen-commented and worth a look if you have not done so already.
// Most important system assembly steps, linear and non-linear resolutions, and time-stepping mechanisms
// are implemented at the SolverBase class level and can thus be used in Laplace without needing reimplementation.

//START_SPHINX_INCLUDE_BEGINCLASS
class DieterichSeismicityRate : public SeismicityRateBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  DieterichSeismicityRate() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the
  /// tree structure of classes)
  DieterichSeismicityRate( const string & name,
                           Group * const parent );

  /// Destructor
  virtual ~DieterichSeismicityRate() override;

  /// "CatalogName()" return the string used as XML tag in the input file.  It ties the XML tag with
  /// this C++ classes. This is important.
  static string catalogName() { return "DieterichSeismicityRate"; }

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override; 

//END_SPHINX_INCLUDE_BEGINCLASS
// /**
//  * @defgroup Solver Interface Functions
//  *
//  * These functions provide the primary interface that is required for derived classes
//  */
// /**@{*/

//START_SPHINX_INCLUDE_SOLVERINTERFACE
  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  void odeSolverStep( real64 const & time_n,
                      real64 const & dt,
                      integer const cycleNumber,
                      DomainPartition & domain );

  void integralSolverStep( real64 const & time_n,
                      real64 const & dt,
                      integer const cycleNumber,
                      DomainPartition & domain );

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               ParallelVector & rhs,
               ParallelVector & solution,
               bool const setSparsity = false ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void initializePreSubGroups() override;

  struct viewKeyStruct : public SeismicityRateBase::viewKeyStruct
  {
    static constexpr char const * directEffect() { return "directEffect"; }
    static constexpr char const * bStressRate() { return "bStressRate"; }
  };

private:
  real64 m_directEffect;
  real64 m_bStressRate;
  real64 m_initialSigma;
  real64 m_initialTau;


//END_SPHINX_INCLUDE_SOLVERINTERFACE

  /**@}*/

};
} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SIMPLEPDE_LAPLACE_FEM_HPP_ */
