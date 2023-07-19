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
  
};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_SEISMICITY_RATE_BASE_HPP */
