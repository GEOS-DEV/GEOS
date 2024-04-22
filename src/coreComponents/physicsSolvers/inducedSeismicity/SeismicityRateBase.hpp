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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITY_RATE_BASE_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITY_RATE_BASE_HPP

#include "codingUtilities/EnumStrings.hpp"   // facilities for enum-string conversion (for reading enum values from XML input)
#include "physicsSolvers/SolverBase.hpp"  // an abstraction class shared by all physics solvers
#include "fieldSpecification/FieldSpecificationManager.hpp" // a manager that can access and set values on the discretized domain

#include "physicsSolvers/inducedSeismicity/inducedSeismicityFields.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

namespace geos
{

class SeismicityRateBase : public SolverBase
{
public:
  /// The default nullary constructor is disabled to avoid compiler auto-generation:
  SeismicityRateBase() = delete;

  /// The constructor needs a user-defined "name" and a parent Group (to place this instance in the tree structure of classes)
  SeismicityRateBase( const string & name,
                      Group * const parent );

  /// Destructor
  virtual ~SeismicityRateBase() override;

  /// This method ties properties with their supporting mesh
  virtual void registerDataOnMesh( Group & meshBodies ) override;

  struct viewKeyStruct : public SolverBase::viewKeyStruct
  {
    static constexpr char const * stressSolverNameString() { return "stressSolverName"; }
    static constexpr char const * initialFaultNormalTractionString() { return "initialFaultNormalTraction"; }
    static constexpr char const * initialFaultShearTractionString() { return "initialFaultShearTraction"; }
    static constexpr char const * faultNormalDirectionString() { return "faultNormalDirection"; }
    static constexpr char const * faultShearDirectionString() { return "faultShearDirection"; }
  };

protected:

  /**
   * @brief called in SolverStep after member stress solver is called to
   *  project the stress state to pre-defined fault orientations
   * @param subRegion The ElementSubRegionBase that will have the stress information
   */
  void updateFaultTraction( ElementSubRegionBase & subRegion );

  /**
   * @brief called in SolverStep before member stress solver is called to
   *  project the initial stress state to pre-defined fault orientations,
   *  only when cycleNumber == 0
   * @param cycleNumber current cycle number
   * @param domain The DomainPartition of the problem
   */
  void initializeFaultTraction( real64 const time_n, integer const cycleNumber, DomainPartition & domain );

  void constructFaultStressProjectionTensors(
    real64 ( &faultNormalProjectionTensor )[6], real64 ( &faultShearProjectionTensor )[6] );

  virtual void postProcessInput() override;

  /// pointer to stress solver
  SolverBase * m_stressSolver;

  /// stress solver name string
  string m_stressSolverName;

  /// intial stress conditions
  real64 m_initialFaultNormalTraction;
  real64 m_initialFaultShearTraction;

  /// fault orientation
  R1Tensor m_faultNormalDirection;
  R1Tensor m_faultShearDirection;
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITY_RATE_BASE_HPP */
