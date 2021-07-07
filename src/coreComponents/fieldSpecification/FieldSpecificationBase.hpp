/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FieldSpecificationBase.hpp
 */

#ifndef GEOSX_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP
#define GEOSX_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP

#include "common/DataTypes.hpp"
#include "codingUtilities/traits.hpp"
#include "codingUtilities/Utilities.hpp"
#include "dataRepository/Group.hpp"
#include "functions/FunctionBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "functions/FunctionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{
class Function;


/**
 * @class FieldSpecificationBase
 * A class to hold values for and administer a single boundary condition
 */
class FieldSpecificationBase : public dataRepository::Group
{
public:

  /**
   * @defgroup alias and functions to defined statically initialized catalog
   * @{
   */

  /**
   * alias to define the catalog type for this base type
   */
  using CatalogInterface = dataRepository::CatalogInterface< FieldSpecificationBase,
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
  {
    return FieldSpecificationBase::catalogName();
  }

  /**
   * @}
   */


  /**
   * @brief constructor
   * @param name the name of the FieldSpecificationBase in the data repository
   * @param parent the parent group of this group.
   */
  FieldSpecificationBase( string const & name, dataRepository::Group * parent );

  /**
   * destructor
   */
  virtual ~FieldSpecificationBase() override;


  /// Deleted copy constructor
  FieldSpecificationBase( FieldSpecificationBase const & ) = delete;

  /// Defaulted move constructor
  FieldSpecificationBase( FieldSpecificationBase && ) = default;

  /// deleted copy assignment
  FieldSpecificationBase & operator=( FieldSpecificationBase const & ) = delete;

  /// deleted move assignement
  FieldSpecificationBase & operator=( FieldSpecificationBase && ) = delete;



  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] field the field to apply the value to.
   * @param[in] targetSet the set of indices which the value will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the value.
   * @param[in] dataGroup the Group that contains the field to apply the value to.
   *
   * This function applies the value to a field variable.
   */
  template< typename FIELD_OP, typename POLICY, typename T, int N, int USD >
  void applyFieldValueKernel( ArrayView< T, N, USD > const & field,
                              SortedArrayView< localIndex const > const & targetSet,
                              real64 const time,
                              Group & dataGroup ) const;

  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] targetSet the set of indices which the value will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the value.
   * @param[in] dataGroup the Group that contains the field to apply the value to.
   * @param[in] fieldname the name of the field to apply the value to.
   *
   * This function applies the value to a field variable. This function is typically
   * called from within the lambda to a call to FieldSpecificationManager::applyFieldValue().
   */
  template< typename FIELD_OP, typename POLICY=parallelHostPolicy >
  void applyFieldValue( SortedArrayView< localIndex const > const & targetSet,
                        real64 const time,
                        dataRepository::Group & dataGroup,
                        string const & fieldname ) const;

  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] fieldName The name of the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
   * @param[in] dofDim The number of degrees of freedom per index of the primary field. For instance
   *                   this will be 1 for a pressure degree of freedom, and 3 for a displacement
   *                   degree of freedom.
   * @param[inout] matrix A ParallelMatrix object: the system matrix.
   * @param[inout] rhs A ParallelVector object: the right-hand side.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename LAI >
  void applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                       real64 const time,
                                       dataRepository::Group const & dataGroup,
                                       string const & fieldName,
                                       string const & dofMapName,
                                       integer const & dofDim,
                                       typename LAI::ParallelMatrix & matrix,
                                       typename LAI::ParallelVector & rhs ) const;


  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *               Either \ref OpEqual or \ref OpAdd.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofDim The number of degrees of freedom per index of the primary field. For instance
   *                   this will be 1 for a pressure degree of freedom, and 3 for a displacement
   *                   degree of freedom.
   * @param[inout] matrix A ParallelMatrix object: the system matrix.
   * @param[inout] rhs A ParallelVector object: the right-hand side
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                provided by the FIELD_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename LAI, typename LAMBDA >
  void
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  integer const & dofDim,
                                  typename LAI::ParallelMatrix & matrix,
                                  typename LAI::ParallelVector & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *               Either \ref OpEqual or \ref OpAdd.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dt time step size which is applied as a factor to bc values
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofDim The number of degrees of freedom per index of the primary field. For instance
   *                   this will be 1 for a pressure degree of freedom, and 3 for a displacement
   *                   degree of freedom.
   * @param[inout] matrix A ParallelMatrix object: the system matrix.
   * @param[inout] rhs A ParallelVector object: the right-hand side
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                provided by the FIELD_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename LAI, typename LAMBDA >
  void
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  integer const & dofDim,
                                  typename LAI::ParallelMatrix & matrix,
                                  typename LAI::ParallelVector & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Zero matrix rows to apply boundary conditions.
   * @tparam LAI The linear algebra interface
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[inout] matrix A ParallelMatrix object: the system matrix.
   *
   * This function zeroes the rows of the matrix that correspond to boundary conditions.
   */
  template< typename LAI >
  void zeroSystemRowsForBoundaryCondition( SortedArrayView< localIndex const > const & targetSet,
                                           arrayView1d< globalIndex const > const & dofMap,
                                           typename LAI::ParallelMatrix & matrix ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over the target set.
   * @tparam T Data type of the field.
   * @tparam NDIM Number of dimensions in the field array.
   * @tparam USD Unit stride dimension of the field array.
   * @param targetSet The set of indices which the boundary condition will be applied.
   * @param time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param dataGroup The Group that contains the field to apply the boundary condition to.
   * @param dofMap The map from the local index of the primary field to the global degree of freedom number.
   * @param dofRankOffset Offset of dof indices on current rank.
   * @param matrix Local part of the system matrix.
   * @param rhs Local part of the system rhs vector.
   * @param fieldView Array view of the field data.
   *
   * @note This function is rarely used directly. More often it is called by other ApplyBoundaryCondition functions.
   */
  template< typename FIELD_OP, typename POLICY, typename T, int NDIM, int USD >
  void applyBoundaryConditionToSystemKernel( SortedArrayView< localIndex const > const & targetSet,
                                             real64 const time,
                                             dataRepository::Group const & dataGroup,
                                             arrayView1d< globalIndex const > const & dofMap,
                                             globalIndex const dofRankOffset,
                                             CRSMatrixView< real64, globalIndex const > const & matrix,
                                             arrayView1d< real64 > const & rhs,
                                             ArrayView< T const, NDIM, USD > const & fieldView ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *            application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] fieldName The name of the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[in,out] matrix Local part of the system matrix.
   * @param[in,out] rhs Local part of the system rhs vector.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename POLICY >
  void applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                       real64 const time,
                                       dataRepository::Group const & dataGroup,
                                       string const & fieldName,
                                       string const & dofMapName,
                                       globalIndex const dofRankOffset,
                                       CRSMatrixView< real64, globalIndex const > const & matrix,
                                       arrayView1d< real64 > const & rhs ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[inout] matrix Local part of the system matrix.
   * @param[inout] rhs Local part of the system rhs vector.
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                   provided by the FIELD_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename POLICY, typename LAMBDA >
  void
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dt time step size which is applied as a factor to bc values
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[inout] matrix Local part of the system matrix.
   * @param[inout] rhs Local part of the system rhs vector.
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                   provided by the FIELD_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename FIELD_OP, typename POLICY, typename LAMBDA >
  void
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Function to zero matrix rows to apply boundary conditions
   * @tparam POLICY the execution policy to use when zeroing rows
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[inout] matrix the local system matrix
   *
   * This function zeroes the rows of the matrix that correspond to boundary conditions.
   */
  template< typename POLICY >
  void zeroSystemRowsForBoundaryCondition( SortedArrayView< localIndex const > const & targetSet,
                                           arrayView1d< globalIndex const > const & dofMap,
                                           CRSMatrixView< real64, globalIndex const > const & matrix ) const;

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
    /// @return The key for bcApplicationTableName
    constexpr static char const * bcApplicationTableNameString() { return "bcApplicationTableName"; }
    /// @return The key for scale
    constexpr static char const * scaleString() { return "scale"; }
    /// @return The key for functionName
    constexpr static char const * functionNameString() { return "functionName"; }
    /// @return The key for initialCondition
    constexpr static char const * initialConditionString() { return "initialCondition"; }
    /// @return The key for beginTime
    constexpr static char const * beginTimeString() { return "beginTime"; }
    /// @return The key for endTime
    constexpr static char const * endTimeString() { return "endTime"; }
    /// @return The key for fluxBoundaryCondition
    constexpr static char const * fluxBoundaryConditionString() { return "fluxBoundaryCondition"; }
  };

  /**
   * Accessor
   * @return const reference to m_function
   */
  string const & getFunctionName() const
  {
    return m_functionName;
  }

  /**
   * Accessor
   * @return const reference to m_objectPath
   */
  virtual const string & getObjectPath() const
  {
    return m_objectPath;
  }

  /**
   * Accessor
   * @return const reference to m_fieldName
   */
  virtual const string & getFieldName() const
  {
    return m_fieldName;
  }

  /**
   * Accessing the considered component.
   * @return The component axis or a special value.
   */
  virtual int getComponent() const
  {
    return m_component;
  }

  /**
   * Accessor
   * @return const reference to m_direction
   */
  virtual R1Tensor const & getDirection() const
  {
    GEOSX_UNUSED_VAR( time );
    return m_direction;
  }

  /**
   * Accessor
   * @return const m_beginTime
   */
  real64 getStartTime() const
  {
    return m_beginTime;
  }

  /**
   * Accessor
   * @return const m_endTime
   */
  real64 getEndTime() const
  {
    return m_endTime;
  }

  /**
   * Accessor
   * @return const reference to m_setNames
   */
  string_array const & getSetNames() const
  {
    return m_setNames;
  }

  /**
   * Accessor
   * @return const m_initialCondition
   */
  int initialCondition() const
  {
    return m_initialCondition;
  }

  /**
   * Accessor
   * @return const m_scale
   */
  real64 getScale() const
  {
    return m_scale;
  }

  /**
   * Mutator
   * @param[in] fieldName The name of the field
   */
  void setFieldName( string const & fieldName )
  {
    m_fieldName = fieldName;
  }

  /**
   * Mutator
   * @param[in] objectPath The path for the object
   */
  void setObjectPath( string const & objectPath )
  {
    m_objectPath = objectPath;
  }

  /**
   * Mutator
   * @param[in] scale Scaling factor
   */
  void setScale( real64 const & scale )
  {
    m_scale = scale;
  }

  /**
   * Mutator
   * @param[in] isInitialCondition Logical value to indicate if it is an initial condition
   */
  void initialCondition( bool isInitialCondition )
  {
    m_initialCondition = isInitialCondition;
  }

  /**
   * Mutator
   * @param[in] setName The name of the set
   */
  void addSetName( string const & setName )
  {
    m_setNames.emplace_back( setName );
  }


protected:
  virtual void postProcessInput() override;

  /// The flag used to decide if the BC value is normalized by the size of the set on which it is applied
  bool m_normalizeBySetSize;

private:


  /// the names of the sets that the boundary condition is applied to
  string_array m_setNames;

  /// the path to the object which contains the fields that the boundary condition is applied to
  string m_objectPath;

  /// the name of the field the boundary condition is applied to or a key string to use for
  /// determining whether or not to apply the boundary condition.
  string m_fieldName;


//  string m_dataType;

  /// The component the boundary condition acts on. Not used if field is a scalar.
  int m_component;

  /// The direction the boundary condition acts in.
  R1Tensor m_direction;

  /// Whether or not the boundary condition is an initial condition.
  int m_initialCondition;

  /// The name of the function used to generate values for application.
  string m_functionName;

  /// The scale factor to use on the value of the boundary condition.
  real64 m_scale;

  /// Time after which the bc is allowed to be applied
  real64 m_beginTime;

  /// Time after which the bc will no longer be applied.
  real64 m_endTime;

  /// The name of a function used to turn on and off the boundary condition.
  string m_bcApplicationFunctionName;

  /// The factor used to normalize the boundary flux by the size of the set it is applied to
  //real64 m_setSizeScalingFactor;

};


template< typename FIELD_OP, typename POLICY, typename T, int N, int USD >
void FieldSpecificationBase::applyFieldValueKernel( ArrayView< T, N, USD > const & field,
                                                    SortedArrayView< localIndex const > const & targetSet,
                                                    real64 const time,
                                                    Group & dataGroup ) const
{
  integer const component = getComponent();
  FunctionManager & functionManager = FunctionManager::getInstance();

  if( m_functionName.empty() )
  {
    real64 const value = m_scale;
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      FIELD_OP::SpecifyFieldValue( field, a, component, value );
    } );
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( m_functionName );

    if( function.isFunctionOfTime()==2 )
    {
      real64 const value = m_scale * function.evaluate( &time );
      forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, value );
      } );
    }
    else
    {
      real64_array result( static_cast< localIndex >( targetSet.size() ) );
      function.evaluate( dataGroup, time, targetSet, result );
      arrayView1d< real64 const > const & resultView = result.toViewConst();
      real64 const scale = m_scale;
      forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, scale * resultView[i] );
      } );
    }
  }
}


template< typename FIELD_OP, typename POLICY >
void FieldSpecificationBase::applyFieldValue( SortedArrayView< localIndex const > const & targetSet,
                                              real64 const time,
                                              dataRepository::Group & dataGroup,
                                              string const & fieldName ) const
{
  dataRepository::WrapperBase & wrapper = dataGroup.getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( wrapper.getTypeId());

  rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                  true,
                                  [&]( auto arrayInstance, auto GEOSX_UNUSED_PARAM( dataTypeInstance ) )
  {
    using ArrayType = decltype(arrayInstance);
    dataRepository::Wrapper< ArrayType > & view = dynamicCast< dataRepository::Wrapper< ArrayType > & >( wrapper );

    auto const & field = view.reference().toView();
    applyFieldValueKernel< FIELD_OP, POLICY >( field, targetSet, time, dataGroup );
  } );
}

template< typename FIELD_OP, typename LAI >
void FieldSpecificationBase::applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                                             real64 const time,
                                                             dataRepository::Group const & dataGroup,
                                                             string const & fieldName,
                                                             string const & dofMapName,
                                                             integer const & dofDim,
                                                             typename LAI::ParallelMatrix & matrix,
                                                             typename LAI::ParallelVector & rhs ) const
{
  dataRepository::WrapperBase const & wrapperBase = dataGroup.getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( wrapperBase.getTypeId());
  arrayView1d< globalIndex const > const & dofMap = dataGroup.getReference< array1d< globalIndex > >( dofMapName );
  integer const component = getComponent();

  rtTypes::applyArrayTypeLambda1( rtTypes::typeID( typeIndex ), [&]( auto type )
  {
    using FieldType = decltype( type );
    dataRepository::Wrapper< FieldType > const & wrapper = dynamicCast< dataRepository::Wrapper< FieldType > const & >( wrapperBase );
    traits::ViewTypeConst< FieldType > fieldView = wrapper.reference();

    this->applyBoundaryConditionToSystem< FIELD_OP, LAI >( targetSet, time, dataGroup, dofMap, dofDim, matrix, rhs,
                                                           [&]( localIndex const a )
    {
      real64 value = 0.0;
      FieldSpecificationEqual::readFieldValue( fieldView, a, component, value );
      return value;
    } );
  } );
}

template< typename FIELD_OP, typename LAI, typename LAMBDA >
void
FieldSpecificationBase::
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  integer const & dofDim,
                                  typename LAI::ParallelMatrix & matrix,
                                  typename LAI::ParallelVector & rhs,
                                  LAMBDA && lambda ) const
{
  GEOSX_UNUSED_VAR( dofDim );

  integer const component = getComponent();
  FunctionManager & functionManager = FunctionManager::getInstance();

  globalIndex_array dof( targetSet.size() );
  real64_array rhsContribution( targetSet.size() );

  real64 sizeScalingFactor = 0;
  if( m_normalizeBySetSize )
  {
    // note: this assumes that the ghost elements have been filtered out

    // recompute the set size here to make sure that topology changes are accounted for
    integer const localSetSize = targetSet.size();
    integer globalSetSize = 0;

    // synchronize
    MpiWrapper::allReduce( &localSetSize, &globalSetSize, 1, MPI_SUM, MPI_COMM_GEOSX );

    // set the scaling factor
    sizeScalingFactor = globalSetSize >= 1 ? 1.0 / globalSetSize : 1;
  }
  else
  {
    sizeScalingFactor = 1;
  }

  if( m_functionName.empty() )
  {

    integer counter=0;
    for( auto a : targetSet )
    {
      dof( counter ) = dofMap[a]+component;
      FIELD_OP::template SpecifyFieldValue< LAI >( dof( counter ),
                                                   matrix,
                                                   rhsContribution( counter ),
                                                   m_scale * sizeScalingFactor,
                                                   lambda( a ) );
      ++counter;
    }
    FIELD_OP::template PrescribeRhsValues< LAI >( rhs, counter, dof.data(), rhsContribution.data() );
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( m_functionName );

    if( function.isFunctionOfTime() == 2 )
    {
      real64 value = m_scale * function.evaluate( &time ) * sizeScalingFactor;
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue< LAI >( dof( counter ),
                                                     matrix,
                                                     rhsContribution( counter ),
                                                     value,
                                                     lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues< LAI >( rhs, counter, dof.data(), rhsContribution.data() );
    }
    else
    {
      real64_array result;
      result.resize( LvArray::integerConversion< localIndex >( targetSet.size()));
      function.evaluate( dataGroup, time, targetSet, result );
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue< LAI >( dof( counter ),
                                                     matrix,
                                                     rhsContribution( counter ),
                                                     m_scale * result[counter] * sizeScalingFactor,
                                                     lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues< LAI >( rhs, counter, dof.data(), rhsContribution.data() );
    }
  }
}

template< typename FIELD_OP, typename LAI, typename LAMBDA >
void
FieldSpecificationBase::
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  integer const & dofDim,
                                  typename LAI::ParallelMatrix & matrix,
                                  typename LAI::ParallelVector & rhs,
                                  LAMBDA && lambda ) const
{
  GEOSX_UNUSED_VAR( dofDim );

  integer const component = getComponent();
  FunctionManager & functionManager = FunctionManager::getInstance();

  globalIndex_array dof( targetSet.size() );
  real64_array rhsContribution( targetSet.size() );

  real64 sizeScalingFactor = 0.0;
  if( m_normalizeBySetSize )
  {
    // note: this assumes that the ghost elements have been filtered out

    // recompute the set size here to make sure that topology changes are accounted for
    integer const localSetSize = targetSet.size();
    integer globalSetSize = 0;

    // synchronize
    MpiWrapper::allReduce( &localSetSize, &globalSetSize, 1, MPI_SUM, MPI_COMM_GEOSX );

    // set the scaling factor
    sizeScalingFactor = globalSetSize >= 1 ? 1.0 / globalSetSize : 1;
  }
  else
  {
    sizeScalingFactor = 1;
  }


  if( m_functionName.empty() )
  {

    integer counter=0;
    for( auto a : targetSet )
    {
      dof( counter ) = dofMap[a]+component;
      FIELD_OP::template SpecifyFieldValue< LAI >( dof( counter ),
                                                   matrix,
                                                   rhsContribution( counter ),
                                                   m_scale * dt * sizeScalingFactor,
                                                   lambda( a ) );
      ++counter;
    }
    FIELD_OP::template PrescribeRhsValues< LAI >( rhs, counter, dof.data(), rhsContribution.data() );
  }
  else
  {
    FunctionBase const & function  = functionManager.getGroup< FunctionBase >( m_functionName );

    if( function.isFunctionOfTime() == 2 )
    {
      real64 value = m_scale * dt * function.evaluate( &time ) * sizeScalingFactor;
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue< LAI >( dof( counter ),
                                                     matrix,
                                                     rhsContribution( counter ),
                                                     value,
                                                     lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues< LAI >( rhs, counter, dof.data(), rhsContribution.data() );
    }
    else
    {
      real64_array result;
      result.resize( LvArray::integerConversion< localIndex >( targetSet.size()));
      function.evaluate( dataGroup, time, targetSet, result );
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue< LAI >( dof( counter ),
                                                     matrix,
                                                     rhsContribution( counter ),
                                                     m_scale * dt * result[counter] * sizeScalingFactor,
                                                     lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues< LAI >( rhs, counter, dof.data(), rhsContribution.data() );
    }
  }
}

template< typename LAI >
void FieldSpecificationBase::zeroSystemRowsForBoundaryCondition( SortedArrayView< localIndex const > const & targetSet,
                                                                 arrayView1d< globalIndex const > const & dofMap,
                                                                 typename LAI::ParallelMatrix & matrix ) const

{
  integer const component = getComponent();
  for( auto a : targetSet )
  {
    globalIndex const dof = dofMap[a]+component;
    matrix.clearRow( dof );
  }
}

template< typename FIELD_OP, typename POLICY, typename T, int NDIM, int USD >
void FieldSpecificationBase::applyBoundaryConditionToSystemKernel( SortedArrayView< localIndex const > const & targetSet,
                                                                   real64 const time,
                                                                   dataRepository::Group const & dataGroup,
                                                                   arrayView1d< globalIndex const > const & dofMap,
                                                                   globalIndex const dofRankOffset,
                                                                   CRSMatrixView< real64, globalIndex const > const & matrix,
                                                                   arrayView1d< real64 > const & rhs,
                                                                   ArrayView< T const, NDIM, USD > const & fieldView ) const
{
  integer const component = getComponent();
  this->applyBoundaryConditionToSystem< FIELD_OP, POLICY >( targetSet, time, dataGroup, dofMap, dofRankOffset, matrix, rhs,
                                                            [fieldView, component] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    real64 value = 0.0;
    FieldSpecificationEqual::readFieldValue( fieldView, a, component, value );
    return value;
  } );
}

template< typename FIELD_OP, typename POLICY >
void FieldSpecificationBase::applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                                             real64 const time,
                                                             dataRepository::Group const & dataGroup,
                                                             string const & fieldName,
                                                             string const & dofMapName,
                                                             globalIndex const dofRankOffset,
                                                             CRSMatrixView< real64, globalIndex const > const & matrix,
                                                             arrayView1d< real64 > const & rhs ) const
{
  dataRepository::WrapperBase const & wrapperBase = dataGroup.getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( wrapperBase.getTypeId());
  arrayView1d< globalIndex const > const & dofMap = dataGroup.getReference< array1d< globalIndex > >( dofMapName );

  rtTypes::applyArrayTypeLambda1( rtTypes::typeID( typeIndex ), [&]( auto type )
  {
    using FieldType = decltype( type );
    dataRepository::Wrapper< FieldType > const & wrapper = dynamicCast< dataRepository::Wrapper< FieldType > const & >( wrapperBase );
    applyBoundaryConditionToSystemKernel< FIELD_OP, POLICY >( targetSet, time, dataGroup, dofMap, dofRankOffset, matrix, rhs, wrapper.reference() );
  } );
}

template< typename FIELD_OP, typename POLICY, typename LAMBDA >
void
FieldSpecificationBase::
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const
{
  return applyBoundaryConditionToSystem< FIELD_OP, POLICY >( targetSet,
                                                             time,
                                                             1.0,
                                                             dataGroup,
                                                             dofMap,
                                                             dofRankOffset,
                                                             matrix,
                                                             rhs,
                                                             std::forward< LAMBDA >( lambda ) );
}

template< typename FIELD_OP, typename POLICY, typename LAMBDA >
void
FieldSpecificationBase::
  applyBoundaryConditionToSystem( SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda ) const
{
  integer const component = getComponent();
  string const & functionName = getReference< string >( viewKeyStruct::functionNameString() );
  FunctionManager & functionManager = FunctionManager::getInstance();

  array1d< globalIndex > dofArray( targetSet.size() );
  arrayView1d< globalIndex > const & dof = dofArray.toView();

  array1d< real64 > rhsContributionArray( targetSet.size() );
  arrayView1d< real64 > const & rhsContribution = rhsContributionArray.toView();

  real64 sizeScalingFactor = 1.0;
  if( m_normalizeBySetSize )
  {
    // note: this assumes that the ghost elements have been filtered out
    // recompute the set size here to make sure that topology changes are accounted for
    globalIndex const globalSetSize = MpiWrapper::sum( LvArray::integerConversion< globalIndex >( targetSet.size() ), MPI_COMM_GEOSX );
    sizeScalingFactor = ( globalSetSize > 0 ) ? ( 1.0 / globalSetSize ) : 1.0;
  }


  if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
  {
    real64 value = m_scale * dt * sizeScalingFactor;
    if( !functionName.empty() )
    {
      FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
      value *= function.evaluate( &time );
    }

    forAll< POLICY >( targetSet.size(),
                      [targetSet, dof, dofMap, dofRankOffset, component, matrix, rhsContribution, value, lambda] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      dof[ i ] = dofMap[ a ] + component;
      FIELD_OP::SpecifyFieldValue( dof[ i ],
                                   dofRankOffset,
                                   matrix,
                                   rhsContribution[ i ],
                                   value,
                                   lambda( a ) );
    } );
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );

    real64_array resultsArray( targetSet.size() );
    function.evaluate( dataGroup, time, targetSet, resultsArray );
    arrayView1d< real64 const > const & results = resultsArray.toViewConst();
    real64 const value = m_scale * dt * sizeScalingFactor;

    forAll< POLICY >( targetSet.size(),
                      [targetSet, dof, dofMap, dofRankOffset, component, matrix, rhsContribution, results, value, lambda] GEOSX_HOST_DEVICE (
                        localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      dof[ i ] = dofMap[ a ] + component;
      FIELD_OP::SpecifyFieldValue( dof[ i ],
                                   dofRankOffset,
                                   matrix,
                                   rhsContribution[ i ],
                                   value * results[ i ],
                                   lambda( a ) );
    } );
  }

  FIELD_OP::template prescribeRhsValues< POLICY >( rhs, dof, dofRankOffset, rhsContribution );
}

template< typename POLICY >
void FieldSpecificationBase::zeroSystemRowsForBoundaryCondition( SortedArrayView< localIndex const > const & targetSet,
                                                                 arrayView1d< globalIndex const > const & dofMap,
                                                                 CRSMatrixView< real64, globalIndex const > const & matrix ) const

{
  integer const component = getComponent();
  forAll< POLICY >( targetSet.size(), [targetSet, dofMap, matrix, component] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    globalIndex const dof = dofMap[ a ] + component;

    arraySlice1d< real64 > const entries = matrix.getEntries( dof );
    localIndex const numEntries = matrix.numNonZeros( dof );

    for( localIndex j = 0; j < numEntries; ++j )
    {
      entries[ j ] = 0;
    }
  } );
}

}

#endif //GEOSX_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP
