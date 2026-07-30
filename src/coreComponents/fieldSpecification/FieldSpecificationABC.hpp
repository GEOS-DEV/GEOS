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
 * @file FieldSpecificationABC.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONABC_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONABC_HPP


#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/FunctionBase.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{
class Function;


/**
 * @class FieldSpecificationABC
 *
 * Abstract Base Class for concrete field modifiers (`FieldSpecification`) and high-level user-defined field specification.
 */
class FieldSpecificationABC : public dataRepository::Group
{
public:

  /**
   * @defgroup alias and functions to define statically initialized catalog
   * @{
   */

  /**
   * alias to define the catalog type for this abstract type
   */
  using CatalogInterface = dataRepository::CatalogInterface< FieldSpecificationABC,
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
   * @brief return the catalog name
   * @return the catalog name
   */
  virtual const string getCatalogName() const = 0;

  /**
   * @}
   */


  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationABC in the data repository
   * @param parent the parent group of this group.
   */
  FieldSpecificationABC( string const & name, dataRepository::Group * parent );

  /**
   * destructor
   */
  virtual ~FieldSpecificationABC() override;


  /// Deleted copy constructor
  FieldSpecificationABC( FieldSpecificationABC const & ) = delete;

  /// Defaulted move constructor
  FieldSpecificationABC( FieldSpecificationABC && ) = default;

  /// deleted copy assignment
  FieldSpecificationABC & operator=( FieldSpecificationABC const & ) = delete;

  /// deleted move assignement
  FieldSpecificationABC & operator=( FieldSpecificationABC && ) = delete;

  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr static char const * setNamesString() { return "setNames"; }
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
   * @return const reference to the field name
   */
  virtual const string & getFieldName() const = 0;

  /**
   * Accessor
   * @return const reference to m_setNames
   */
  string_array const & getSetNames() const
  { return m_setNames; }

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
  array1d< real64 > getScales() const
  { return m_scales; }

  /**
   * Accessor
   * @return const m_emptySetErrorMode
   */
  SetErrorMode getErrorSetMode() const
  { return m_emptySetErrorMode; }

  /**
   * Mutator
   * @param[in] setName The name of the set
   */
  void addSetName( string const & setName )
  { m_setNames.emplace_back( setName ); }

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
   * Mutator
   * @param[in] errorSetMode Time after which the bc will no longer be applied.
   */
  void setErrorSetMode( SetErrorMode const & errorSetMode )
  { m_emptySetErrorMode = errorSetMode; }

protected:

  virtual void postInputInitialization() override;

  /// the names of the sets that the boundary condition is applied to
  string_array m_setNames;

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
ENUM_STRINGS( FieldSpecificationABC::SetErrorMode,
              "silent",
              "error",
              "warning" );

}

#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONABC_HPP
