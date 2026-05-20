/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FieldApplicator.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDAPPLICATOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDAPPLICATOR_HPP_

#include "events/tasks/TaskBase.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "functions/FunctionManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

/**
 * @class FieldApplicator
 * @brief Task to apply field specifications to elements during runtime.
 *
 * This task applies field specifications to elements that are created dynamically during
 * the simulation, such as fracture elements created by the SurfaceGenerator. It handles
 * both simple field specifications and complex initialization like HydrostaticEquilibrium
 * for compositional flow, including the computation of derived quantities like component
 * densities.
 */
class FieldApplicator : public TaskBase
{
public:

  /**
   * @brief Constructor for the FieldApplicator class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  FieldApplicator( const string & name,
                   Group * const parent );

  /// Destructor for the class
  ~FieldApplicator() override;

  /// Accessor for the catalog name
  static string catalogName()
  {
    return "FieldApplicator";
  }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/

private:

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * fieldSpecificationNamesString() { return "fieldSpecificationNames"; }
    constexpr static char const * solverNameString() { return "solverName"; }
    constexpr static char const * targetRegionsString() { return "targetRegions"; }
  };

  void postInputInitialization() override;

  /**
   * @brief Return the target field name.
   * @param[in] fieldName Name of the field specification.
   * @return Name of the target field being updated.
   */
  string getTargetFieldName( string const & fieldName ) const;


  /**
   * @brief Initialize the fluid state for a subregion after applying equilibrium fields.
   * @param[in] domain The domain partition.
   * @param[in,out] subRegion The element subregion to initialize.
   *
   * This function is called after applying the HydrostaticEquilibrium fields to properly
   * initialize all derived compositional variables. The full initialization sequence is:
   * 1. Porosity and permeability (fracPorosity_porosity, fracPorosity_initialPorosity,
   *    fracPorosity_referencePorosity, fracPerm_permeability) via updatePorosityAndPermeability.
   *    For fracture elements (SurfaceElementSubRegion), this uses the hydraulic aperture
   *    (defaultAperture from the XML) to compute initial fracture permeability and porosity.
   * 2. Component densities (globalCompDensity) from component fractions and total fluid density.
   * 3. Phase volume fractions (phaseVolumeFraction) via updateFluidState.
   * 4. Relative permeability values (RelPerm_phaseRelPerm) via updateFluidState.
   * 5. Phase mobilities (phaseMobility) via updateFluidState.
   * 6. All '_n' fields (converged state) via saveConvergedState.
   */
  void initializeSubRegionFluidState( DomainPartition & domain, ElementSubRegionBase & subRegion );

  /// Array of field specification names to apply
  stdVector< string > m_fieldSpecificationNames;

  /// Name of the flow solver to use for fluid state initialization
  string m_solverName;

  /// Target regions to apply field specifications to (optional filter)
  stdVector< string > m_targetRegions;
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FIELDAPPLICATOR_HPP_ */
