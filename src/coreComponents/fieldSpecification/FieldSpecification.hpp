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
#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
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
   * @enum  SetErrorMode
   * @brief Indicate the error handling mode.
   */
  enum class SetErrorMode : integer
  {
    silent,
    error,
    warning
  };

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
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * setNamesString() { return "setNames"; }
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
    /// @return The key for scale
    constexpr static char const * scalesString() { return "scale"; }
    /// @return The key for functionName
    constexpr static char const * functionNamesString() { return "functionName"; }
    /// @return The key for initialCondition
    constexpr static char const * initialConditionString() { return "initialCondition"; }
    /// @return The key for beginTime
    constexpr static char const * beginTimeString() { return "beginTime"; }
    /// @return The key for endTime
    constexpr static char const * endTimeString() { return "endTime"; }
    /// @return The key errorSetMode
    constexpr static char const * errorSetModeString() { return "errorSetMode"; }
  };

  /**
   * Accessor
   * @return first entry of m_functionNames, or an empty string if empty
   *
   * @note Legacy scalar accessor.
   *       Use getFunctionNames() to access the full list of function names when using non-scalar
   *       field specifications (eg. functionName="{ f1, f2, f3 }")
   */
  string const & getFunctionName() const
  {
    static string const emptyName;
    return m_functionNames.empty() ? emptyName : m_functionNames.front();
  }

  /**
   * Accessor
   * @return const reference to m_functionNames
   */
  string_array const & getFunctionNames() const
  {
    return m_functionNames;
  }

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
   * @return const reference to m_direction
   */
  virtual R1Tensor const & getDirection() const
  { return m_direction; }

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
   * @return const reference to m_setNames
   */
  string_array const & getSetNames() const
  { return m_setNames; }

  /**
   * Accessor
   * @return const m_initialCondition
   */
  int initialCondition() const
  { return m_initialCondition; }

  /**
   * Accessor
   * @return first entry of m_scales, or 0 if m_scales is empty
   *
   * @note Legacy scalar accessor.
   *       Use getScales() to access the full list of scales when using non-scalar
   *       field specifications (eg. scales="{ 1, 2, 3 }")
   */
  real64 getScale() const
  { return m_scales.empty() ? 0.0 : m_scales.front(); }

  /**
   * Accessor
   * @return const m_scales
   */
  arrayView1d< real64 const > getScales() const
  { return m_scales.toViewConst(); }

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
   * Mutator
   * @param[in] scale Scaling factor
   */
  void setScale( real64 const & scale )
  {
    m_scales.resize( 1 );
    m_scales[ 0 ] = scale;
  }

  /**
   * Mutator
   * @brief Set the per-component scale factors
   * @param[in] scales The tensor-valued scale
   */
  void setScales( array1d< real64 > const & scales )
  { m_scales = scales; }

  /**
   * Mutator
   * @brief Set the per-component function names
   * @param[in] functionNames The per-component function names. Must either be empty,
   *                          have a single entry, or be sized exactly as @p m_scales
   */
  void setFunctionNames( string_array const & functionNames )
  { m_functionNames = functionNames; }

  /**
   * Mutator
   * @param[in] isInitialCondition Logical value to indicate if it is an initial condition
   */
  void initialCondition( bool isInitialCondition )
  { m_initialCondition = isInitialCondition; }

  /**
   * Mutator
   * @param[in] setName The name of the set
   */
  void addSetName( string const & setName )
  { m_setNames.emplace_back( setName ); }

  /**
   * @brief Set the Mesh Object Path object
   *
   * @param meshBodies The group containing all the MeshBody objects
   */
  void setMeshObjectPath( Group const & meshBodies );

  /**
   * Mutator
   * @param[in] beginTime Time after which the bc is allowed to be applied
   */
  void setStartTime( real64 beginTime )
  { m_beginTime = beginTime; }

  /**
   * Mutator
   * @param[in] endTime Time after which the bc will no longer be applied.
   */
  void setEndTime( real64 endTime )
  { m_endTime = endTime; }

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


  /// the names of the sets that the boundary condition is applied to
  string_array m_setNames;

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
  SetErrorMode m_emptySetErrorMode;
};

/**
 * @brief Indicate the error handling mode
 */
ENUM_STRINGS( FieldSpecification::SetErrorMode,
              "silent",
              "error",
              "warning" );

}

#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATION_HPP
