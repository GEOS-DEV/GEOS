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
 * @file PermeabilitySpecification.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATION_HPP
#define GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATION_HPP


#include "common/DataTypes.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/MeshObjectPath.hpp"
#include "FieldSpecification.hpp"
#include "FieldSpecificationABC.hpp"
#include "FieldSpecificationFactory.hpp"

namespace geos
{

/**
 * @class PermeabilitySpecification
 *
 * Data class representing a permeability field specification
 */
class PermeabilitySpecification : public FieldSpecificationABC
{
public:

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "PermeabilitySpecification"; }

  /**
   * @brief return the catalog name
   * @return the catalog name
   */
  virtual const string getCatalogName() const
  {
    return PermeabilitySpecification::catalogName();
  }


  /**
   * @brief constructor
   * @param name the name of the PermeabilitySpecification in the data repository
   * @param parent the parent group of this group.
   */
  PermeabilitySpecification( string const & name, dataRepository::Group * parent );

  /**
   * destructor
   */
  virtual ~PermeabilitySpecification() override;


  /// Deleted copy constructor
  PermeabilitySpecification( PermeabilitySpecification const & ) = delete;

  /// Defaulted move constructor
  PermeabilitySpecification( PermeabilitySpecification && ) = default;

  /// deleted copy assignment
  PermeabilitySpecification & operator=( PermeabilitySpecification const & ) = delete;

  /// deleted move assignement
  PermeabilitySpecification & operator=( PermeabilitySpecification && ) = delete;

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * setNamesString() { return "setNames"; }
    /// @return The key for regionNames
    constexpr static char const * regionNamesString() { return "regionNames"; }
    /// @return The key for fieldName
    constexpr static char const * fieldNameString() { return "fieldName"; }
    /// @return The key for initialCondition
    constexpr static char const * initialConditionString() { return "initialCondition"; }
    /// @return The key for functionNames
    constexpr static char const * functionNamesString() { return "functionNames"; }
    /// @return The key for scales
    constexpr static char const * scalesString() { return "scales"; }
    /// @return The key for beginTime
    constexpr static char const * beginTimeString() { return "beginTime"; }
    /// @return The key for endTime
    constexpr static char const * endTimeString() { return "endTime"; }
    /// @return The key errorSetMode
    constexpr static char const * errorSetModeString() { return "errorSetMode"; }
  };

  /**
   * Accessor
   * @return const reference to m_functionNames
   */
  string_array const & getFunctionNames() const
  { return m_functionNames; }

  /**
   * Accessor
   * @return const reference to m_regionNames
   */
  string_array const & getRegionNames() const
  { return m_regionNames; }

  /**
   * Accessor
   * @return const reference to m_fieldName
   */
  virtual const string & getFieldName() const
  { return m_fieldName; }

  /**
   * Accessor
   * @return const reference to m_setNames
   */
  string_array const & getSetNames() const
  { return m_setNames; }

  /**
   * Accessor
   * @return const m_scales
   */
  array1d< real64 > getScales() const
  { return m_scales; }

  /**
   * Accessor
   * @return const m_initialCondition
   */
  int initialCondition() const
  { return m_initialCondition; }

  /**
   * Accessor
   * @return const m_beginTime
   */
  real64 getStartTime() const
  { return m_beginTime; }

  /**
   * Accessor
   * @return const m_endTime
   */
  real64 getEndTime() const
  { return m_endTime; }

  /**
   * Accessor
   * @return const m_emptySetErrorMode
   */
  FieldSpecification::SetErrorMode getErrorSetMode() const
  { return m_emptySetErrorMode; }

protected:

  virtual void postInputInitialization() override;

private:


  /// the names of the sets that the boundary condition is applied to
  string_array m_setNames;

  /// the names of the regions that the boundary condition is applied to
  string_array m_regionNames;

  /// the name of the field the boundary condition is applied to or a key string to use for
  /// determining whether or not to apply the boundary condition.
  string m_fieldName;

  /// Whether or not the boundary condition is an initial condition.
  int m_initialCondition;

  /// Name(s) of the function used to generate values for application.
  string_array m_functionNames;

  /// Scale factor(s) to use on the value of the boundary condition.
  array1d< real64 > m_scales;

  /// Time after which the bc is allowed to be applied
  real64 m_beginTime;

  /// Time after which the bc will no longer be applied.
  real64 m_endTime;

  /// Enum containing the possible output modes when an error occur
  FieldSpecification::SetErrorMode m_emptySetErrorMode;

};

/// @copydoc geos::FieldSpecificationFactory::generateFieldSpecifications
template<>
void
generateFieldSpecifications< PermeabilitySpecification >( PermeabilitySpecification const & specification,
                                                          dataRepository::Group & manager );

}

#endif //GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATION_HPP
