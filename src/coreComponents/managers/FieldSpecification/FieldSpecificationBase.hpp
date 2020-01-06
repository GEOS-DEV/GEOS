/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FieldBase.hpp
 */

#ifndef GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP
#define GEOSX_MANAGERS_FIELDSPECIFICATION_FIELDSPECIFICATIONBASE_HPP

#include "common/DataTypes.hpp"
#include "codingUtilities/traits.hpp"
#include "codingUtilities/Utilities.hpp"
#include "dataRepository/Group.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "managers/FieldSpecification/FieldSpecificationOps.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "managers/ObjectManagerBase.hpp"


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
  static CatalogInterface::CatalogType& GetCatalog();

  static string CatalogName() { return "FieldSpecification"; }

  virtual const string getCatalogName() const
  {
    return FieldSpecificationBase::CatalogName();
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

  template < typename FIELD_OP, typename POLICY, typename T, int N, int UNIT_STRIDE_DIM >
  void ApplyFieldValueKernel( LvArray::ArrayView< T, N, UNIT_STRIDE_DIM, localIndex > const & field,
                              SortedArrayView< localIndex const > const & targetSet,
                              real64 const time,
                              Group * dataGroup ) const;

  /**
   * @tparam FIELD_OP type that contains static functions to apply the value to the field
   * @param[in] targetSet the set of indices which the value will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the value.
   * @param[in] dataGroup the Group that contains the field to apply the value to.
   * @param[in] fieldname the name of the field to apply the value to.
   *
   * This function applies the value to a field variable. This function is typically
   * called from within the lambda to a call to FieldSpecificationManager::ApplyFieldValue().
   */
  template< typename FIELD_OP, typename POLICY=parallelHostPolicy >
  void ApplyFieldValue( set<localIndex> const & targetSet,
                        real64 const time,
                        dataRepository::Group * dataGroup,
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
  void ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                       real64 const time,
                                       dataRepository::Group * dataGroup,
                                       string const & fieldName,
                                       string const & dofMapName,
                                       integer const & dofDim,
                                       typename LAI::ParallelMatrix & matrix,
                                       typename LAI::ParallelVector & rhs ) const;


  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *               Either \ref BcEqual or \ref BcAdd.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
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
  ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                  real64 const time,
                                  dataRepository::Group * dataGroup,
                                  arrayView1d<globalIndex const> const & dofMap,
                                  integer const & dofDim,
                                  typename LAI::ParallelMatrix & matrix,
                                  typename LAI::ParallelVector & rhs,
                                  LAMBDA && lambda ) const;

  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @tparam FIELD_OP A wrapper struct to define how the boundary condition operates on the variables.
   *               Either \ref BcEqual or \ref BcAdd.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dt time step size which is applied as a factor to bc values
   * @param[in] dataGroup The Group that contains the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
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
  ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                  real64 const time,
                                  real64 const dt,
                                  dataRepository::Group * dataGroup,
                                  arrayView1d<globalIndex const> const & dofMap,
                                  integer const & dofDim,
                                  typename LAI::ParallelMatrix & matrix,
                                  typename LAI::ParallelVector & rhs,
                                  LAMBDA && lambda ) const;

  
  
  struct viewKeyStruct
  {
    constexpr static auto setNamesString = "setNames";
    constexpr static auto constitutivePathString = "constitutivePath";
    constexpr static auto objectPathString = "objectPath";
    constexpr static auto fieldNameString = "fieldName";
    constexpr static auto dataTypeString = "dataType";
    constexpr static auto componentString = "component";
    constexpr static auto directionString = "direction";
    constexpr static auto bcApplicationTableNameString = "bcApplicationTableName";
    constexpr static auto scaleString = "scale";
    constexpr static auto functionNameString = "functionName";
    constexpr static auto initialConditionString = "initialCondition";
    constexpr static auto beginTimeString = "beginTime";
    constexpr static auto endTimeString = "endTime";
    constexpr static auto fluxBoundaryConditionString = "fluxBoundaryConditionString"; 
  } viewKeys;

  struct groupKeyStruct
  {} groupKeys;


  /**
   * Accessor
   * @return const reference to m_function
   */
  string const & GetFunctionName() const
  {
    return m_functionName;
  }

  virtual const string& GetObjectPath() const
  {
    return m_objectPath;
  }

  virtual const string& GetFieldName() const
  {
    return m_fieldName;
  }

  virtual int GetComponent() const
  {
    return m_component;
  }

  virtual const R1Tensor& GetDirection( realT GEOSX_UNUSED_ARG( time ) )
  {
    return m_direction;
  }

  real64 GetStartTime() const
  {
    return m_beginTime;
  }

  real64 GetEndTime() const
  {
    return m_endTime;
  }

  string_array const & GetSetNames() const
  {
    return m_setNames;
  }

  int initialCondition() const
  {
    return m_initialCondition;
  }

  real64 GetScale() const
  { return m_scale; }

  void SetFieldName( string const & fieldName )
  {
    m_fieldName = fieldName;
  }

  void SetObjectPath( string const & objectPath )
  {
    m_objectPath = objectPath;
  }

  void SetScale( real64 const & scale )
  {
    m_scale = scale;
  }

  void InitialCondition( bool isInitialCondition)
  {
    m_initialCondition = isInitialCondition;
  }

  void AddSetName( string const & setName )
  {
    m_setNames.push_back( setName );
  }


protected:
  void PostProcessInput() override final;

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


template < typename FIELD_OP, typename POLICY, typename T, int N, int UNIT_STRIDE_DIM >
void FieldSpecificationBase::ApplyFieldValueKernel( LvArray::ArrayView< T, N, UNIT_STRIDE_DIM, localIndex > const & field,
                                                    SortedArrayView< localIndex const > const & targetSet,
                                                    real64 const time,
                                                    Group * dataGroup ) const
{
  integer const component = GetComponent();
  string const & functionName = getReference<string>( viewKeyStruct::functionNameString );
  FunctionManager & functionManager = FunctionManager::Instance();

  if( functionName.empty() )
  {
    forall_in_range< POLICY >( 0, targetSet.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      FIELD_OP::SpecifyFieldValue( field, a, component, m_scale );
    });
  }
  else
  {
    FunctionBase const * const function  = functionManager.GetGroup<FunctionBase>( functionName );

    GEOSX_ERROR_IF( function == nullptr, "Function '" << functionName << "' not found" );

    if( function->isFunctionOfTime()==2 )
    {
      real64 value = m_scale * function->Evaluate( &time );
      forall_in_range< POLICY >( 0, targetSet.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, value );
      });
    }
    else
    {
      real64_array result( static_cast<localIndex>( targetSet.size() ) );
      function->Evaluate( dataGroup, time, targetSet, result );
      arrayView1d<real64 const> const & resultView = result;
      forall_in_range< POLICY >( 0, targetSet.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        FIELD_OP::SpecifyFieldValue( field, a, component, m_scale*resultView[i] );
      });
    }
  }
}


template< typename FIELD_OP, typename POLICY >
void FieldSpecificationBase::ApplyFieldValue( set<localIndex> const & targetSet,
                                              real64 const time,
                                              Group * dataGroup,
                                              string const & fieldName ) const
{
  dataRepository::WrapperBase * wrapper = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( wrapper->get_typeid());

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( typeIndex ),
                                 false,
                                 [&]( auto arrayInstance, auto GEOSX_UNUSED_ARG( dataTypeInstance ) )
  {
    using ArrayType = decltype(arrayInstance);
    dataRepository::Wrapper<ArrayType> & view = dataRepository::Wrapper<ArrayType>::cast( *wrapper );

    typename ArrayType::ViewType const & field = view.referenceAsView();
    ApplyFieldValueKernel< FIELD_OP, POLICY >( field, targetSet, time, dataGroup );
  });
}

template< typename FIELD_OP, typename LAI >
void FieldSpecificationBase::ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                                             real64 const time,
                                                             dataRepository::Group * dataGroup,
                                                             string const & fieldName,
                                                             string const & dofMapName,
                                                             integer const & dofDim,
                                                             typename LAI::ParallelMatrix & matrix,
                                                             typename LAI::ParallelVector & rhs ) const
{
  dataRepository::WrapperBase * wrapperBase = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( wrapperBase->get_typeid());
  arrayView1d<globalIndex> const & dofMap = dataGroup->getReference<array1d<globalIndex>>( dofMapName );
  integer const component = GetComponent();

  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ),
    [&]( auto type ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::Wrapper<fieldType> & wrapper = dynamic_cast< dataRepository::Wrapper<fieldType> & >(*wrapperBase);
      typename dataRepository::Wrapper<fieldType>::ViewTypeConst fieldView = wrapper.referenceAsView();

      this->ApplyBoundaryConditionToSystem<FIELD_OP, LAI>( targetSet, time, dataGroup, dofMap, dofDim, matrix, rhs,
        [&]( localIndex const a )->real64
        {
          real64 value = 0.0;
          FieldSpecificationEqual::ReadFieldValue( fieldView, a, component, value );
          return value;
        }
      );
    }
  );
}

template< typename FIELD_OP, typename LAI, typename LAMBDA >
void
FieldSpecificationBase::
ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                real64 const time,
                                dataRepository::Group * dataGroup,
                                arrayView1d<globalIndex const> const & dofMap,
                                integer const & GEOSX_UNUSED_ARG( dofDim ),
                                typename LAI::ParallelMatrix & matrix,
                                typename LAI::ParallelVector & rhs,
                                LAMBDA && lambda ) const
{
  integer const component = GetComponent();
  string const & functionName = getReference<string>( viewKeyStruct::functionNameString );
  FunctionManager & functionManager = FunctionManager::Instance();

  globalIndex_array  dof( targetSet.size() );
  real64_array rhsContribution( targetSet.size() );

  real64 sizeScalingFactor = 0;
  if (m_normalizeBySetSize)
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
  
  if( functionName.empty() )
  {

    integer counter=0;
    for( auto a : targetSet )
    {
      dof( counter ) = dofMap[a]+component;
      FIELD_OP::template SpecifyFieldValue<LAI>( dof( counter ),
                                                 matrix,
                                                 rhsContribution( counter ),
                                                 m_scale * sizeScalingFactor,
                                                 lambda( a ) );
      ++counter;
    }
    FIELD_OP::template PrescribeRhsValues<LAI>( rhs, counter, dof.data(), rhsContribution.data() );
  }
  else
  {
    FunctionBase const * const function  = functionManager.GetGroup<FunctionBase>( functionName );

    GEOSX_ERROR_IF( function == nullptr, "Function '" << functionName << "' not found" );

    if( function->isFunctionOfTime()==2 )
    {
      real64 value = m_scale * function->Evaluate( &time ) * sizeScalingFactor;
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue<LAI>( dof( counter ),
                                                   matrix,
                                                   rhsContribution( counter ),
                                                   value,
                                                   lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues<LAI>( rhs, counter, dof.data(), rhsContribution.data() );
    }
    else
    {
      real64_array result;
      result.resize( integer_conversion<localIndex>( targetSet.size()));
      function->Evaluate( dataGroup, time, targetSet, result );
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue<LAI>( dof( counter ),
                                                   matrix,
                                                   rhsContribution( counter ),
                                                   m_scale * result[counter] * sizeScalingFactor,
                                                   lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues<LAI>( rhs, counter, dof.data(), rhsContribution.data() );
    }
  }
}

template< typename FIELD_OP, typename LAI, typename LAMBDA >
void
FieldSpecificationBase::
ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                real64 const time,
                                real64 const dt,
                                dataRepository::Group * dataGroup,
                                arrayView1d<globalIndex const> const & dofMap,
                                integer const & GEOSX_UNUSED_ARG( dofDim ),
                                typename LAI::ParallelMatrix & matrix,
                                typename LAI::ParallelVector & rhs,
                                LAMBDA && lambda ) const
{
  integer const component = GetComponent();
  string const & functionName = getReference<string>( viewKeyStruct::functionNameString );
  FunctionManager & functionManager = FunctionManager::Instance();

  globalIndex_array  dof( targetSet.size() );
  real64_array rhsContribution( targetSet.size() );

  real64 sizeScalingFactor = 0.0;
  if (m_normalizeBySetSize)
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

  
  if( functionName.empty() )
  {

    integer counter=0;
    for( auto a : targetSet )
    {
      dof( counter ) = dofMap[a]+component;
      FIELD_OP::template SpecifyFieldValue<LAI>( dof( counter ),
                                                 matrix,
                                                 rhsContribution( counter ),
                                                 m_scale * dt * sizeScalingFactor,
                                                 lambda( a ) );
      ++counter;
    }
    FIELD_OP::template PrescribeRhsValues<LAI>( rhs, counter, dof.data(), rhsContribution.data() );
  }
  else
  {
    FunctionBase const * const function  = functionManager.GetGroup<FunctionBase>( functionName );

    GEOSX_ERROR_IF( function == nullptr, "Function '" << functionName << "' not found" );

    if( function->isFunctionOfTime()==2 )
    {
      real64 value = m_scale * dt * function->Evaluate( &time ) * sizeScalingFactor;
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue<LAI>( dof( counter ),
                                                   matrix,
                                                   rhsContribution( counter ),
                                                   value,
                                                   lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues<LAI>( rhs, counter, dof.data(), rhsContribution.data() );
    }
    else
    {
      real64_array result;
      result.resize( integer_conversion<localIndex>( targetSet.size()));
      function->Evaluate( dataGroup, time, targetSet, result );
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofMap[a] + component;
        FIELD_OP::template SpecifyFieldValue<LAI>( dof( counter ),
                                                   matrix,
                                                   rhsContribution( counter ),
                                                   m_scale * dt * result[counter] * sizeScalingFactor,
                                                   lambda( a ) );
        ++counter;
      }
      FIELD_OP::template PrescribeRhsValues<LAI>( rhs, counter, dof.data(), rhsContribution.data() );
    }
  }
}

}
#endif
