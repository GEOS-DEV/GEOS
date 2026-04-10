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
 * @file FieldSpecificationImpl.hpp
 */

#ifndef GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONIMPL_HPP
#define GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONIMPL_HPP

#include "FieldSpecification.hpp"
#include "common/DataTypes.hpp"
#include "common/TypeDispatch.hpp"
#include "codingUtilities/traits.hpp"
#include "codingUtilities/Utilities.hpp"
#include "dataRepository/Group.hpp"
#include "functions/FunctionBase.hpp"
#include "common/FieldSpecificationOps.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/MeshObjectPath.hpp"
#include "functions/FunctionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{
class Function;


/**
 * @class FieldSpecificationImpl
 * @brief Methods to apply field specifications
 */
class FieldSpecificationImpl
{
public:

  /**
   * @brief Apply this field specification to the discretization
   *
   * @tparam OBJECT_TYPE The type of discretization/mesh object that the
   *   specification is being applied to.
   * @tparam BC_TYPE The type of BC being applied
   * @tparam LAMBDA
   * @param fs The field specification data object
   * @param mesh The MeshLevel that the specification is applied to
   * @param lambda The being executed
   */
  template< typename OBJECT_TYPE,
            typename BC_TYPE = FieldSpecification,
            typename LAMBDA >
  static void apply( BC_TYPE const & fs,
                     MeshLevel & mesh,
                     LAMBDA && lambda )
  {
    MeshObjectPath const & meshObjectPaths = fs.getMeshObjectPaths();
    meshObjectPaths.forObjectsInPath< OBJECT_TYPE >( mesh,
                                                     [&] ( OBJECT_TYPE & object )
    {
      {
        dataRepository::Group const & setGroup = object.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
        string_array setNames = fs.getSetNames();
        for( auto & setName : setNames )
        {
          if( setGroup.hasWrapper( setName ) )
          {
            SortedArrayView< localIndex const > const & targetSet = setGroup.getReference< SortedArray< localIndex > >( setName );
            lambda( fs, setName, targetSet, object, fs.getFieldName() );
          }
        }
      }
    } );
  }

  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] fs the field specfication data object.
   * @param[in] field the field to apply the value to.
   * @param[in] targetSet the set of indices which the value will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the value.
   * @param[in] dataGroup the Group that contains the field to apply the value to.
   *
   * This function applies the value to a field variable.
   */
  template< typename FIELD_OP, typename POLICY, typename T, int N, int USD >
  static void applyFieldValueKernel( FieldSpecification const & fs,
                                     ArrayView< T, N, USD > const & field,
                                     SortedArrayView< localIndex const > const & targetSet,
                                     real64 const time,
                                     dataRepository::Group & dataGroup )
  {
    integer const component = fs.getComponent();
    string const & functionName = fs.getFunctionName();
    FunctionManager & functionManager = FunctionManager::getInstance();

    if( functionName.empty() )
    {
      real64 const value = fs.getScale();
      forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, value );
      } );
    }
    else
    {
      FunctionBase const & function = [&]() -> FunctionBase const &
      {
        try
        {
          return functionManager.getGroup< FunctionBase >( functionName );
        }
        catch( std::exception const & e )
        {
          string const errorMsg = GEOS_FMT( "Error while reading {}:\n",
                                            fs.getWrapperDataContext( FieldSpecification::viewKeyStruct::functionNameString() ) );
          ErrorLogger::global().modifyCurrentExceptionMessage()
            .addToMsg( errorMsg )
            .addContextInfo( fs.getWrapperDataContext( FieldSpecification::viewKeyStruct::functionNameString() )
                               .getContextInfo()
                               .setPriority( 1 ) );
          throw InputError( e, errorMsg );
        }
      }();

      if( function.isFunctionOfTime()==2 )
      {
        real64 const value = fs.getScale() * function.evaluate( &time );
        forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
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
        real64 const scale = fs.getScale();
        forAll< POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
        {
          localIndex const a = targetSet[ i ];
          FIELD_OP::SpecifyFieldValue( field, a, component, scale * resultView[i] );
        } );
      }
    }
  }


  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] fs the field specification data object
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
  static void applyFieldValue( FieldSpecification const & fs,
                               SortedArrayView< localIndex const > const & targetSet,
                               real64 const time,
                               dataRepository::Group & dataGroup,
                               string const & fieldName )
  {
    dataRepository::WrapperBase & wrapper = dataGroup.getWrapperBase( fieldName );

    // This function is used in setting boundary/initial conditions on simulation fields.
    // This is meaningful for 1/2/3D real arrays and sometimes 1D integer (indicator) arrays.
    using FieldTypes = types::ListofTypeList< types::Join< types::ArrayTypes< types::RealTypes, types::DimsUpTo< 3 > >,
                                                           types::ArrayTypes< types::TypeList< integer >, types::DimsSingle< 1 > > > >;


    types::dispatch( FieldTypes{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;
      auto & wrapperT = dataRepository::Wrapper< ArrayType >::cast( wrapper );
      applyFieldValueKernel< FIELD_OP, POLICY >( fs, wrapperT.reference().toView(), targetSet, time, dataGroup );
    }, wrapper );
  }

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over the target set.
   * @tparam T Data type of the field.
   * @tparam NDIM Number of dimensions in the field array.
   * @tparam USD Unit stride dimension of the field array.
   * @param fs The field specification data object.
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
  static void applyBoundaryConditionToSystemKernel( FieldSpecification const & fs,
                                                    SortedArrayView< localIndex const > const & targetSet,
                                                    real64 const time,
                                                    dataRepository::Group const & dataGroup,
                                                    arrayView1d< globalIndex const > const & dofMap,
                                                    globalIndex const dofRankOffset,
                                                    CRSMatrixView< real64, globalIndex const > const & matrix,
                                                    arrayView1d< real64 > const & rhs,
                                                    ArrayView< T const, NDIM, USD > const & fieldView )
  {
    integer const component = fs.getComponent();
    applyBoundaryConditionToSystem< FIELD_OP, POLICY >( fs, targetSet, time, dataGroup, dofMap, dofRankOffset, matrix, rhs,
                                                        [fieldView, component] GEOS_HOST_DEVICE ( localIndex const a )
    {
      real64 value = 0.0;
      FieldSpecificationEqual::readFieldValue( fieldView, a, component, value );
      return value;
    } );
  }

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @param[in] fs The field specification data object
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
  static void applyBoundaryConditionToSystem( FieldSpecification const & fs,
                                              SortedArrayView< localIndex const > const & targetSet,
                                              real64 const time,
                                              dataRepository::Group const & dataGroup,
                                              string const & fieldName,
                                              string const & dofMapName,
                                              globalIndex const dofRankOffset,
                                              CRSMatrixView< real64, globalIndex const > const & matrix,
                                              arrayView1d< real64 > const & rhs )
  {
    dataRepository::WrapperBase const & wrapper = dataGroup.getWrapperBase( fieldName );
    arrayView1d< globalIndex const > const & dofMap = dataGroup.getReference< array1d< globalIndex > >( dofMapName );

    // We're reading values from a field, which is only well-defined for dims 1 and 2
    using FieldTypes = types::ListofTypeList< types::ArrayTypes< types::RealTypes, types::DimsUpTo< 2 > > >;
    types::dispatch( FieldTypes{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;
      auto const & wrapperT = dataRepository::Wrapper< ArrayType >::cast( wrapper );
      applyBoundaryConditionToSystemKernel< FIELD_OP, POLICY >( fs,
                                                                targetSet,
                                                                time,
                                                                dataGroup,
                                                                dofMap,
                                                                dofRankOffset,
                                                                matrix,
                                                                rhs,
                                                                wrapperT.reference() );
    }, wrapper );
  }

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] fs The field specification data object.
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
  static void
  applyBoundaryConditionToSystem( FieldSpecification const & fs,
                                  SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda )
  {
    return applyBoundaryConditionToSystem< FIELD_OP, POLICY >( fs,
                                                               targetSet,
                                                               time,
                                                               1.0,
                                                               dataGroup,
                                                               dofMap,
                                                               dofRankOffset,
                                                               matrix,
                                                               rhs,
                                                               std::forward< LAMBDA >( lambda ) );
  }

  /**
   * @brief Apply a boundary condition to a system of equations.
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] fs The field specification data object.
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
  static void
  applyBoundaryConditionToSystem( FieldSpecification const & fs,
                                  SortedArrayView< localIndex const > const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group const & dataGroup,
                                  arrayView1d< globalIndex const > const & dofMap,
                                  globalIndex const dofRankOffset,
                                  CRSMatrixView< real64, globalIndex const > const & matrix,
                                  arrayView1d< real64 > const & rhs,
                                  LAMBDA && lambda )
  {
    array1d< globalIndex > dofArray( targetSet.size() );
    arrayView1d< globalIndex > const & dof = dofArray.toView();

    array1d< real64 > rhsContributionArray( targetSet.size() );
    arrayView1d< real64 > const & rhsContribution = rhsContributionArray.toView();

    computeRhsContribution< FIELD_OP, POLICY, LAMBDA >( fs,
                                                        targetSet,
                                                        time,
                                                        dt,
                                                        dataGroup,
                                                        dofMap,
                                                        dofRankOffset,
                                                        matrix,
                                                        dof,
                                                        rhsContribution,
                                                        std::forward< LAMBDA >( lambda ) );

    FIELD_OP::template prescribeRhsValues< POLICY >( rhs, dof, dofRankOffset, rhsContribution );
  }

  /**
   * @brief Compute the contributions that will be added/enforced to the right-hand side, and collect the corresponding dof numbers
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *                  Either \ref OpEqual or \ref OpAdd.
   * @tparam POLICY Execution policy to use when iterating over target set.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] fs The field specification data object.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dt time step size which is applied as a factor to bc values
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[in] dofRankOffset Offset of dof indices on current rank.
   * @param[inout] matrix Local part of the system matrix.
   * @param[inout] dof array storing the degrees of freedom of the rhsContribution, to know where in the rhs they will be added/enforced
   * @param[inout] rhsContribution array storing the values that will be added/enforced to the right-hand side
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                   provided by the FIELD_OP templated type.
   *
   * Note that this function only computes the rhs contributions, but does not apply them to the right-hand side.
   * The application of these rhs contributions is done in applyBoundaryConditionToSystem.
   *
   * Why did we have to extract the computation of the rhs contributions from applyBoundaryConditionToSystem?
   * Because applyBoundaryConditionToSystem is not very well suited to apply the rhsContributions to the equation layout used in the
   * compositional solvers.
   * Therefore, the compositional solvers do not call applyBoundaryConditionToSystem, but instead call computeRhsContribution directly, and
   * apply these rhs contributions "manually" according to the equation layout used in the solver
   */
  template< typename FIELD_OP, typename POLICY, typename LAMBDA >
  static void
  computeRhsContribution( FieldSpecification const & fs,
                          SortedArrayView< localIndex const > const & targetSet,
                          real64 const time,
                          real64 const dt,
                          dataRepository::Group const & dataGroup,
                          arrayView1d< globalIndex const > const & dofMap,
                          globalIndex const dofRankOffset,
                          CRSMatrixView< real64, globalIndex const > const & matrix,
                          arrayView1d< globalIndex > const & dof,
                          arrayView1d< real64 > const & rhsContribution,
                          LAMBDA && lambda )
  {
    integer const component = ( fs.getComponent() >= 0 ) ? fs.getComponent() : 0;
    // string const & functionName = fs.getReference< string >( fs.viewKeyStruct::functionNameString() );
    string const & functionName = fs.getFunctionName();
    FunctionManager & functionManager = FunctionManager::getInstance();

    // Compute the value of the rhs terms, and collect the dof numbers
    // The rhs terms will be assembled in applyBoundaryConditionToSystem (or in the solver for CompositionalMultiphaseBase)

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      real64 value = fs.getScale() * dt;
      if( !functionName.empty() )
      {
        FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
        value *= function.evaluate( &time );
      }

      forAll< POLICY >( targetSet.size(),
                        [targetSet, dof, dofMap, dofRankOffset, component, matrix, rhsContribution, value, lambda] GEOS_HOST_DEVICE ( localIndex const i )
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
      real64 const value = fs.getScale() * dt;

      forAll< POLICY >( targetSet.size(),
                        [targetSet, dof, dofMap, dofRankOffset, component, matrix, rhsContribution, results, value, lambda] GEOS_HOST_DEVICE ( localIndex const i )
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
  }

  /**
   * @brief Function to zero matrix rows to apply boundary conditions
   * @tparam POLICY the execution policy to use when zeroing rows
   * @param[in] fs The field specification data object
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] dofMap The map from the local index of the primary field to the global degree of
   *                   freedom number.
   * @param[inout] matrix the local system matrix
   *
   * This function zeroes the rows of the matrix that correspond to boundary conditions.
   */
  template< typename POLICY >
  static void zeroSystemRowsForBoundaryCondition( FieldSpecification const & fs,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  arrayView1d< globalIndex const > const & dofMap,
                                                  CRSMatrixView< real64, globalIndex const > const & matrix )
  {
    integer const component = ( fs.getComponent() >= 0 ) ? fs.getComponent() : 0;
    forAll< POLICY >( targetSet.size(), [targetSet, dofMap, matrix, component] GEOS_HOST_DEVICE ( localIndex const i )
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

};

} /* namespace geos */

#endif //GEOS_FIELDSPECIFICATION_FIELDSPECIFICATIONIMPL_HPP
