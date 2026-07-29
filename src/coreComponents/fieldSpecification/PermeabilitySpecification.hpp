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


  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationABC::viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * setNamesString() { return "setNames"; }
    /// @return The key for regionNames
    constexpr static char const * regionNamesString() { return "regionNames"; }
    /// @return The key for fieldName
    constexpr static char const * fieldNameString() { return "fieldName"; }
    /// @return The key for component
    constexpr static char const * componentString() { return "component"; }
  };

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
   * Accessing the considered component.
   * @return The component axis or a special value.
   */
  virtual int getComponent() const
  { return m_component; }

  /**
   * Accessor
   * @return const reference to m_setNames
   */
  string_array const & getSetNames() const
  { return m_setNames; }

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

  /// The component the boundary condition acts on. Not used if field is a scalar.
  int m_component;

};

/// @copydoc geos::FieldSpecificationFactory::generateFieldSpecifications
template<>
void expandFieldSpecification< PermeabilitySpecification >( PermeabilitySpecification const & fs,
                                                            dataRepository::Group & manager );

}

#endif //GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATION_HPP
