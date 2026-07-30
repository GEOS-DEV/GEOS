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
#include "constitutive/permeability/PermeabilityFields.hpp"
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
 *
 * @todo Currently the PermeabilitySpecification only supports cells.
 *       A future goal is to support faces too.
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
    /// @return The key for permeabilityModelName
    constexpr static char const * permeabilityModelNameString() { return "permeabilityModelName"; }
  };

  /**
   * Accessor
   * @return const reference to m_regionNames
   */
  string_array const & getRegionNames() const
  { return m_regionNames; }

  /**
   * Accessor
   * @return const reference to m_permeabilityModelName
   */
  string const & getPermeabilityModelName() const
  { return m_permeabilityModelName; }

  /**
   * Accessor
   * @return const reference to m_fieldName
   */
  virtual string const & getFieldName() const override
  {
    static string const fieldName = fields::permeability::permeability::key();
    return fieldName;
  }

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

  /// the name of the constitutive permeability model
  string m_permeabilityModelName;

};

/// @copydoc geos::FieldSpecificationFactory::generateFieldSpecifications
template<>
void expandFieldSpecification< PermeabilitySpecification >( PermeabilitySpecification const & fs,
                                                            dataRepository::Group & manager );

}

#endif //GEOS_FIELDSPECIFICATION_PERMEABILITYSPECIFICATION_HPP
