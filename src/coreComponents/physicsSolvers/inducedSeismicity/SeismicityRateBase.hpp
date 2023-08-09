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
    static constexpr char const * faultNormalString() { return "faultNormal"; }
    static constexpr char const * faultShearString() { return "faultShear"; }
    static constexpr char const * initialSigmaString() { return "initialSigma"; }
    static constexpr char const * initialTauString() { return "initialTau"; }
  };

  virtual void initializePreSubGroups() override;

  void updateMeanSolidStress( ElementSubRegionBase & subRegion );

  void initializeMeanSolidStress( integer const cycleNumber, DomainPartition & domain );
  
protected:

  virtual void postProcessInput() override;

  void initializeFaultOrientation();

  /// pointer to stress solver
  SolverBase * m_stressSolver;

  /// stress solver name string
  string m_stressSolverName;

  /// intial stress conditions
  real64 m_initialSigma;
  real64 m_initialTau;

  /// fault orientation
  array1d< real64 > m_faultNormal;
  array1d< real64 > m_faultShear;
  real64 m_faultNormalVoigt[6];
  real64 m_faultShearVoigt[6];
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITY_RATE_BASE_HPP */
