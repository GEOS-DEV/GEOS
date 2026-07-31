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
 * @file FieldSpecification.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATION_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATION_HPP


#include "common/DataTypes.hpp"
#include "fieldSpecification/FieldSpecificationABC.hpp"
#include "mesh/MeshObjectPath.hpp"

namespace geos
{
class Function;


/**
 * @class FieldSpecification
 * A class to hold values for and administer a single boundary condition
 */
class FieldSpecification : public FieldSpecificationABC
{
public:

  /**
   * @defgroup alias and functions to defined statically initialized catalog
   * @{
   */

  /**
   * alias to define the catalog type for this base type
   */
  using CatalogInterface = dataRepository::CatalogInterface< FieldSpecification,
                                                             string const &,
                                                             dataRepository::Group * const >;

  /**
   * @brief static function to return static catalog.
   * @return the static catalog to create derived types through the static factory methods.
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string catalogName() { return "FieldSpecification"; }

  /**
   * @brief return the catalog name
   * @return the catalog name
   */
  virtual const string getCatalogName() const
  { return FieldSpecification::catalogName(); }

  /**
   * @}
   */


  /**
   * @brief constructor
   * @param name the name of the FieldSpecification in the data repository
   * @param parent the parent group of this group.
   */
  FieldSpecification( string const & name, dataRepository::Group * parent );

  /**
   * destructor
   */
  virtual ~FieldSpecification() override;


  /// Deleted copy constructor
  FieldSpecification( FieldSpecification const & ) = delete;

  /// Defaulted move constructor
  FieldSpecification( FieldSpecification && ) = default;

  /// deleted copy assignment
  FieldSpecification & operator=( FieldSpecification const & ) = delete;

  /// deleted move assignement
  FieldSpecification & operator=( FieldSpecification && ) = delete;

  /**
   * @brief View keys
   */
  struct viewKeyStruct : public FieldSpecificationABC::viewKeyStruct
  {
    /// @return The key for constitutivePath
    constexpr static char const * constitutivePathString() { return "constitutivePath"; }
    /// @return The key for objectPath
    constexpr static char const * objectPathString() { return "objectPath"; }
    /// @return The key for fieldName
    constexpr static char const * fieldNameString() { return "fieldName"; }
    /// @return The key for dataType
    constexpr static char const * dataTypeString() { return "dataType"; }
    /// @return The key for component
    constexpr static char const * componentString() { return "component"; }
    /// @return The key for direction
    constexpr static char const * directionString() { return "direction"; }
  };

  /**
   * Accessor
   * @return const reference to m_objectPath
   */
  virtual const string & getObjectPath() const
  { return m_objectPath; }

  /**
   * Accessor
   * @return const reference to m_fieldName
   */
  virtual const string & getFieldName() const override
  { return m_fieldName; }

  /**
   * Accessing the considered component.
   * @return The component axis or a special value.
   */
  virtual int getComponent() const
  { return m_component; }

  /**
   * Accessor
   * @return const reference to m_direction
   */
  virtual R1Tensor const & getDirection() const
  { return m_direction; }

  /**
   * Mutator
   * @param[in] fieldName The name of the field
   */
  void setFieldName( string const & fieldName )
  { m_fieldName = fieldName; }

  /**
   * Mutator
   * @param[in] component The component axis or a special value.
   */
  void setComponent( int component )
  { m_component = component; }

  /**
   * Mutator
   * @param[in] objectPath The path for the object
   */
  void setObjectPath( string const & objectPath )
  { m_objectPath = objectPath; }

  /**
   * @brief Set the Mesh Object Path object
   *
   * @param meshBodies The group containing all the MeshBody objects
   */
  void setMeshObjectPath( Group const & meshBodies );

  /**
   * @brief Get the Mesh Object Paths object
   *
   * @return reference to const m_meshObjectPaths
   */
  MeshObjectPath const & getMeshObjectPaths() const
  { return *(m_meshObjectPaths.get()); }

  /**
   * @brief Validate that the size of @p m_scales and @p m_functionNames correspond to the
   *        size of the targeted field or expand them by duplicating values if possible.
   *
   * Validate that @p m_scales has the same size as the targeted field.
   * If @p m_scales as a single value and the targeted field expect multiple, @p m_scales will
   * be resized to the size of the field and its values be duplicated.
   * Else, if there is a size mismatch and @p m_scales has more than one value, it throws.
   * (The same applies for @p m_functionNames)
   *
   * @note This method can mutate the FieldSpecification by resizing its @p m_scales and
   *       its @p m_functionNames arrays
   */
  void validateNumArrayComp( localIndex numComp );

protected:

  virtual void postInputInitialization() override;


private:


  /// the path to the object which contains the fields that the boundary condition is applied to
  string m_objectPath;

  std::unique_ptr< MeshObjectPath > m_meshObjectPaths;

  /// the name of the field the boundary condition is applied to or a key string to use for
  /// determining whether or not to apply the boundary condition.
  string m_fieldName;

  /// The component the boundary condition acts on. Not used if field is a scalar.
  int m_component;

  /// The direction the boundary condition acts in.
  R1Tensor m_direction;
};

}

#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATION_HPP
