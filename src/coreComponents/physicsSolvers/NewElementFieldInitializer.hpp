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
 * @file NewElementFieldInitializer.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_NEWELEMENTFIELDINITIALIZER_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_NEWELEMENTFIELDINITIALIZER_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{

class EquilibriumInitialCondition;
class FunctionManager;
class FlowSolverBase;
class ElementSubRegionBase;


/**
 * @class NewElementFieldInitializer
 * @brief Task to initialize fields on newly created elements (e.g., from SurfaceGenerator).
 *
 * This task applies field specifications to elements that are created dynamically during
 * the simulation, such as fracture elements created by the SurfaceGenerator. It handles
 * both simple field specifications and complex initialization like HydrostaticEquilibrium
 * for compositional flow, including the computation of derived quantities like component
 * densities.
 */
class NewElementFieldInitializer : public TaskBase
{
public:

  /**
   * @brief Constructor for the initialization class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  NewElementFieldInitializer( const string & name,
                              Group * const parent );

  /// Destructor for the class
  ~NewElementFieldInitializer() override;

  /// Accessor for the catalog name
  static string catalogName()
  {
    return "NewElementFieldInitializer";
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
  };

  void postInputInitialization() override;

  /**
   * @brief Return the target field name.
   * @param[in] fieldName Name of the field specification.
   * @return Name of the target field being updated.
   */
  string getTargetFieldName( string const & fieldName ) const;

  /**
   * @brief Apply pressure, temperature and component fraction fields for HydrostaticEquilibrium.
   * @param[in] equilIC The EquilibriumInitialCondition (HydrostaticEquilibrium) to apply.
   * @param[in] targetSet The set of indices to apply the fields to.
   * @param[in,out] targetGroup The group containing the fields to update.
   * @param[in] functionManager The function manager containing the table functions.
   *
   * For compositional multiphase flow, HydrostaticEquilibrium sets pressure, temperature,
   * and component fractions. This function handles the temperature and component fraction
   * fields using the elevation-based tables specified in the EquilibriumInitialCondition.
   */
  void applyEquilibriumInitialConditionFields( EquilibriumInitialCondition const & equilIC,
                                               SortedArrayView< localIndex const > const & targetSet,
                                               dataRepository::Group & targetGroup,
                                               FunctionManager & functionManager );

  /**
   * @brief Initialize the fluid state for a subregion.
   * @param[in] domain The domain partition.
   * @param[in,out] subRegion The element subregion to initialize.
   *
   * This function is called after applying the HydrostaticEquilibrium fields to properly
   * initialize all the derived compositional variables (component densities, phase fractions, etc.).
   */
  void initializeSubRegionFluidState( DomainPartition & domain, ElementSubRegionBase & subRegion );

  /// Array of field specification names to apply
  stdVector< string > m_fieldSpecificationNames;

  /// Name of the flow solver to use for fluid state initialization
  string m_solverName;
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_NEWELEMENTFIELDINITIALIZER_HPP_ */
