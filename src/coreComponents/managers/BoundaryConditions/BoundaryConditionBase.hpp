/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file BoundaryConditionBase.hpp
 */

#ifndef BOUNDARYCONDITIONBASE_H
#define BOUNDARYCONDITIONBASE_H

#include "common/DataTypes.hpp"
#include "codingUtilities/GeosxTraits.hpp"
#include "codingUtilities/Utilities.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "fileIO/vtm/VtmFile.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{
class Function;


/**
 * @struct BcEqual
 * this struct a collection of static functions which adhere to an assumed interface for overwriting
 * a value for a boundary condition.
 */
struct BcEqual
{
  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component not used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index] = value.
   */
  template< typename T >
  static inline typename std::enable_if< !traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array1d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    field[index] = static_cast<T>(value);
  }

  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index][component] = value.
   */
  template< typename T >
  static inline typename std::enable_if< traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array1d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    field[index].Data()[component] = value;
  }


  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component not used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index] = value for all values of field[index].
   */
  template< typename T >
  static inline typename std::enable_if< !traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array2d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    for( localIndex a=0 ; a<field.size( 1 ) ; ++a )
    {
      field[index][a] = static_cast<T>(value);
    }
  }

  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index][component] = value for all values of field[index].
   */
  template< typename T >
  static inline typename std::enable_if< traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array2d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    for( localIndex a=0 ; a<field.size( 1 ) ; ++a )
    {
      field[index][a] = value;
    }
  }

  /**
   * @brief Function to apply a Dirichlet like boundary condition to a single dof in a system of
   *        equations.
   * @param[in] dof The degree of freedom that is to be set.
   * @param[in] blockSystem A pointer to the block system object.
   * @param[in] blockID The value of the blockID that contains the specific \p dof being set.
   * @param[out] rhs The rhs contribution resulting from the application of the BC.
   * @param[in] bcValue The target value of the Boundary Condition
   * @param[in] fieldValue The current value of the variable to be set.
   *
   * This function clears the rows in all blocks for the specified \p dof, sets the diagonal to some
   * appropriate scaled value, and sets \p rhs to the product of the scaled value of the diagonal and
   * the difference between \p bcValue and \p fieldValue.
   */
  static inline void ApplyBcValue( globalIndex const dof,
                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            systemSolverInterface::BlockIDs const blockID,
                            real64 & rhs,
                            real64 const & bcValue,
                            real64 const fieldValue )
  {

    if( true )//node_is_ghost[*nd] < 0 )
    {
      real64 LARGE = blockSystem->ClearSystemRow( blockID, static_cast< int >( dof ), 1.0 );
      rhs = -LARGE*( bcValue - fieldValue );
    }
    else
    {
      blockSystem->ClearSystemRow( blockID, static_cast< int >( dof ), 0.0 );
      rhs = 0.0;
    }
  }


  /**
   * @brief Function to replace some values of a vector.
   * @param rhs A pointer to the global vector
   * @param num The number of values in \p rhs to replace
   * @param dof A pointer to the global DOF to be replaced
   * @param values A pointer to the values corresponding to \p dof that will be replaced in \p rhs.
   */
  static inline void ReplaceGlobalValues( Epetra_FEVector * const rhs,
                                          int const num,
                                          globalIndex const * const dof,
                                          real64 const * const values )
  {
    rhs->ReplaceGlobalValues( num, dof, values );
  }
};

/**
 * @struct BcAdd
 * this struct a collection of static functions which adhere to an assumed interface for adding
 * a value for a boundary condition.
 */
struct BcAdd
{
  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component not used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index] += value.
   */
  template< typename T >
  static inline typename std::enable_if< !traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array1d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    field[index] += static_cast<T>(value);
  }

  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array1d field variable specified in @p field.
   * @param[in] field The array1d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index][component] += value.
   */
  template< typename T >
  static inline typename std::enable_if< traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array1d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    field[index].Data()[component] += value;
  }

  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component not used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index] += value for all values of field[index].
   */
  template< typename T >
  static inline typename std::enable_if< !traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array2d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    for( localIndex a=0 ; a<field.size( 1 ) ; ++a )
    {
      field[index][a] += static_cast<T>(value);
    }
  }

  /**
   * @brief Pointwise application of a boundary condition to a field variable.
   * @tparam T The type of the array2d field variable specified in @p field.
   * @param[in] field The array2d field variable to apply @p value to.
   * @param[in] index The index in field to apply @p value to.
   * @param[in] component The component of @p field to apply @p value to. If @p T is a scalar type,
   *                      this will not be used.
   * @param[in] value The value of the boundary condition to apply to @p field.
   *
   * This function performs field[index][component] += value for all values of field[index].
   */
  template< typename T >
  static inline typename std::enable_if< traits::is_tensorT<T>::value, void>::type
  ApplyBcValue( array2d<T> & field,
                localIndex const index,
                int const component,
                real64 const & value )
  {
    for( localIndex a=0 ; a<field.size( 1 ) ; ++a )
    {
      field[index][a] += value;
    }
  }

  /**
   * @brief Function to apply a boundary condition to a vector for a single dof.
   * @param[in] dof The degree of freedom that is to be modified.
   * @param[in] blockSystem A pointer to the block system object.
   * @param[in] blockID The value of the blockID that contains the specific \p dof being modified.
   * @param[out] rhs The rhs contribution to be modified
   * @param[in] bcValue The value to add to rhs
   * @param[in] fieldValue unused.
   *
   */
  static inline void ApplyBcValue( globalIndex const dof,
                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                            systemSolverInterface::BlockIDs const blockID,
                            real64 & rhs,
                            real64 const & bcValue,
                            real64 const fieldValue )
  {
    if( true )//node_is_ghost[*nd] < 0 )
    {
      rhs += bcValue;
    }
  }

  /**
   * @brief Function to add some values of a vector.
   * @param rhs A pointer to the global vector
   * @param num The number of values in \p rhs to replace
   * @param dof A pointer to the global DOF to be replaced
   * @param values A pointer to the values corresponding to \p dof that will be added to \p rhs.
   */
  static inline void ReplaceGlobalValues( Epetra_FEVector * const rhs,
                                          int const num,
                                          globalIndex * const dof,
                                          real64 * const values )
  {
    rhs->SumIntoGlobalValues( num, dof, values );
  }


};


/**
 * @class BoundaryConditionBase
 * A class to hold values for and administer a single boundary condition
 */
class BoundaryConditionBase : public dataRepository::ManagedGroup
{
public:

  /**
   * @defgroup alias and functions to defined statically initialized catalog
   * @{
   */

  /**
   * alias to define the catalog type for this base type
   */
  using CatalogInterface = cxx_utilities::CatalogInterface< BoundaryConditionBase,
                                                            string const &,
                                                            dataRepository::ManagedGroup * const >;

  /**
   * @brief static function to return static catalog.
   * @return the static catalog to create derived types through the static factory methods.
   */
  static CatalogInterface::CatalogType& GetCatalog();

  /**
   * @}
   */


  /**
   * @brief constructor
   * @param name the name of the BoundaryConditionBase in the data repository
   * @param parent the parent group of this group.
   */
  BoundaryConditionBase( string const & name, dataRepository::ManagedGroup * parent );

  /**
   * destructor
   */
  virtual ~BoundaryConditionBase() override;


  void FillDocumentationNode() override;

  void ReadXML_PostProcess() override final;

  /**
   * @tparam BC_OP type that contains static functions to apply the boundary condition to the field
   * @param[in] targetSet the set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup the ManagedGroup that contains the field to apply the boundary condition to.
   * @param[in] fieldname the name of the field to apply the boundary condition to.
   *
   * This function applies the boundary condition to a field variable. This function is typically
   * called from within the lambda to a call to BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename BC_OP >
  void ApplyBoundaryConditionToField( set<localIndex> const & targetSet,
                                      real64 const time,
                                      dataRepository::ManagedGroup * dataGroup,
                                      string const & fieldname ) const;

  // calls user-provided lambda to apply computed boundary value
//  template<typename LAMBDA>
//  void ApplyBoundaryCondition( set<localIndex> const & targetSet,
//                               real64 const time,
//                               dataRepository::ManagedGroup * dataGroup,
//                               LAMBDA && lambda );

  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The ManagedGroup that contains the field to apply the boundary condition to.
   * @param[in] fieldName The name of the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
   * @param[in] dofDim The number of degrees of freedom per index of the primary field. For instance
   *                   this will be 1 for a pressure degree of freedom, and 3 for a displacement
   *                   degree of freedom.
   * @param[in] blockSystem A pointer to a blockSystem object.
   * @param[in] blockID The blockID in the linear system where the global rows of the global degrees
   *                    of freedom are stored.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename BC_OP >
  void ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                       real64 const time,
                                       dataRepository::ManagedGroup * dataGroup,
                                       string const & fieldName,
                                       string const & dofMapName,
                                       integer const & dofDim,
                                       systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                       systemSolverInterface::BlockIDs const blockID ) const;


  /**
   * @brief Function to apply a boundary condition to a system of equations
   * @tparam BC_OP A wrapper struct to define how the boundary condition operates on the variables.
   *               Either \ref BcEqual or \ref BcAdd.
   * @tparam LAMBDA The type of lambda function passed into the parameter list.
   * @param[in] targetSet The set of indices which the boundary condition will be applied.
   * @param[in] time The time at which any time dependent functions are to be evaluated as part of the
   *             application of the boundary condition.
   * @param[in] dataGroup The ManagedGroup that contains the field to apply the boundary condition to.
   * @param[in] dofMapName The name of the map from the local index of the primary field to the
   *                       global degree of freedom number.
   * @param[in] dofDim The number of degrees of freedom per index of the primary field. For instance
   *                   this will be 1 for a pressure degree of freedom, and 3 for a displacement
   *                   degree of freedom.
   * @param[in] blockSystem A pointer to a blockSystem object.
   * @param[in] blockID The blockID in the linear system where the global rows of the global degrees
   *                    of freedom are stored.
   * @param[in] lambda A lambda function which defines how the value that is passed into the functions
   *                provided by the BC_OP templated type.
   *
   * This function applies the boundary condition to a linear system of equations. This function is
   * typically called from within the lambda to a call to
   * BoundaryConditionManager::ApplyBoundaryCondition().
   */
  template< typename BC_OP, typename LAMBDA >
  void
  ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                  real64 const time,
                                  dataRepository::ManagedGroup * dataGroup,
                                  globalIndex_array const & dofMap,
                                  integer const & dofDim,
                                  systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                  systemSolverInterface::BlockIDs const blockID,
                                  LAMBDA && lambda ) const;

  struct viewKeyStruct
  {
    constexpr static auto setNamesString = "setNames";
    constexpr static auto constitutivePathString = "constitutivePath";
    constexpr static auto objectPathString = "objectPath";
    constexpr static auto fieldNameString = "fieldName";
    constexpr static auto fieldNameInFileString = "fieldNameInFile";
    constexpr static auto propertyFileNameString = "propertyFileName";
    constexpr static auto dataTypeString = "dataType";
    constexpr static auto componentString = "component";
    constexpr static auto directionString = "direction";
    constexpr static auto bcApplicationTableNameString = "bcApplicationTableName";
    constexpr static auto scaleString = "scale";
    constexpr static auto functionNameString = "functionName";
    constexpr static auto initialConditionString = "initialCondition";
    constexpr static auto beginTimeString = "beginTime";
    constexpr static auto endTimeString = "endTime";


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

  virtual const R1Tensor& GetDirection( realT time )
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

  string const & GetFieldNameInFile() const {
      return m_fieldNameInFile;
  }

  string const & GetPropertyFileName() const {
      return m_propertyFileName;
  }

  int initialCondition() const
  {
    return m_initialCondition;
  }

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

  /// If not empty, link the boundary condition with an existing property supported in a file
  string m_fieldNameInFile;

  /// Path to the file containing the property
  string m_propertyFileName;

  /// The name of the function used to generate values for application.
  string m_functionName;

  /// The scale factor to use on the value of the boundary condition.
  real64 m_scale;

  /// Time after which the bc is allowed to be applied
  real64 m_beginTime;

  /// Time after which the bc will no longer be applied.
  real64 m_endTime;

  /// the name of a function used to turn on and off the boundary condition.
  string m_bcApplicationFunctionName;


};



template< typename BC_OP >
void BoundaryConditionBase::ApplyBoundaryConditionToField( set<localIndex> const & targetSet,
                                                           real64 const time,
                                                           ManagedGroup * dataGroup,
                                                           string const & fieldName ) const
{
  integer const component = GetComponent();
  string const functionName = getData<string>( viewKeyStruct::functionNameString );
  string const fieldNameInFile = GetFieldNameInFile();
  string const propertyFileName = GetPropertyFileName();
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  dataRepository::ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( vw->get_typeid());

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID( typeIndex ), [&]( auto type, auto baseType ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::ViewWrapper<fieldType> & view = dataRepository::ViewWrapper<fieldType>::cast( *vw );
      fieldType & field = view.reference();
      if( functionName.empty() )
      {
        if( !fieldNameInFile.empty() && !propertyFileName.empty()) {
        // TODO work only for properties which are handled by a vtm file. We should
        // abstract this for other file format.
        VtmFile vtmFile;
        vtmFile.Load(propertyFileName,false,true);
        const auto & propertyMap = vtmFile.GetRankBlock(0).GetMeshBlock(0).PropertyMap();
        const auto & property = propertyMap.find(fieldNameInFile)->second;
        for( auto a : targetSet )
        {
          if( fieldNameInFile == "PERMX" || fieldNameInFile == "PERMY" || fieldNameInFile == "PERMZ") {
          BC_OP::ApplyBcValue( field, a, component, property[a]/1000 );
          }
          else {
          BC_OP::ApplyBcValue( field, a, component, property[a] );
          }
        }
        }
        else {
        for( auto a : targetSet )
        {
          BC_OP::ApplyBcValue( field, a, component, m_scale );
        }
        }
      }
      else
      {
        FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>( functionName );

        GEOS_ASSERT( function != nullptr, "Function '" << functionName << "' not found" );

        if( function->isFunctionOfTime()==2 )
        {
          real64 value = m_scale * function->Evaluate( &time );
          for( auto a : targetSet )
          {
            BC_OP::ApplyBcValue( field, a, component, value );
          }
        }
        else
        {
          real64_array result( static_cast<localIndex>(targetSet.size()) );
          function->Evaluate( dataGroup, time, targetSet, result );
          integer count=0;
          for( auto a : targetSet )
          {
            BC_OP::ApplyBcValue( field, a, component, m_scale*result[count] );
            ++count;
          }
        }
      }
    } );
}



template< typename BC_OP >
void BoundaryConditionBase::ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                                            real64 const time,
                                                            dataRepository::ManagedGroup * dataGroup,
                                                            string const & fieldName,
                                                            string const & dofMapName,
                                                            integer const & dofDim,
                                                            systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                            systemSolverInterface::BlockIDs const blockID ) const
{

#if 1
  dataRepository::ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( vw->get_typeid());
  globalIndex_array const & dofMap = dataGroup->getReference<globalIndex_array>( dofMapName );
  integer const component = GetComponent();

  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ), [&]( auto type ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::ViewWrapper<fieldType> & view = dynamic_cast< dataRepository::ViewWrapper<fieldType> & >(*vw);
      dataRepository::view_rtype<fieldType> field = view.data();

      this->ApplyBoundaryConditionToSystem<BC_OP>( targetSet,
                                                   time,
                                                   dataGroup,
                                                   dofMap,
                                                   dofDim,
                                                   blockSystem,
                                                   blockID,
                                                   [&]( localIndex const a )->real64
      {
        return static_cast<real64>(rtTypes::value( field[a], component ));
      });
    });
#else
  dataRepository::ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index( vw->get_typeid());

  integer const component = GetComponent();
  string const functionName = getData<string>( viewKeyStruct::functionNameString );
  NewFunctionManager * functionManager = NewFunctionManager::Instance();


  integer const numBlocks = blockSystem->numBlocks();
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( blockID );

  array1d<globalIndex>  dof( targetSet.size() );
  array1d<real64>       rhsContribution( targetSet.size() );


  dataRepository::view_rtype_const<globalIndex_array> dofMap = dataGroup->getData<globalIndex_array>( dofMapName );


  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID( typeIndex ), [&]( auto type ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::ViewWrapper<fieldType> & view = dynamic_cast< dataRepository::ViewWrapper<fieldType> & >(*vw);
      dataRepository::view_rtype<fieldType> field = view.data();
      if( functionName.empty() )
      {

        integer counter=0;
        for( auto a : targetSet )
        {
          dof( counter ) = dofDim*dofMap[a]+component;
          BC_OP::ApplyBcValue( dof( counter ),
                               blockSystem,
                               blockID,
                               rhsContribution( counter ),
                               m_scale,
                               static_cast<real64>(rtTypes::value( field[a], component )));
          ++counter;
        }

        BC_OP::ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
      }
      else
      {
        FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>( functionName );
        if( function!=nullptr )
        {
          if( function->isFunctionOfTime()==2 )
          {
            real64 value = m_scale * function->Evaluate( &time );
            integer counter=0;
            for( auto a : targetSet )
            {
              dof( counter ) = dofDim*dofMap[a]+component;
              BC_OP::ApplyBcValue( dof( counter ),
                                   blockSystem,
                                   blockID,
                                   rhsContribution( counter ),
                                   value,
                                   rtTypes::value( field[a], component ));
              ++counter;
            }
            BC_OP::ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
          }
          else
          {
            real64_array result;
            result.resize( integer_conversion<localIndex>( targetSet.size()));
            function->Evaluate( dataGroup, time, targetSet, result );
            integer counter=0;
            for( auto a : targetSet )
            {
              dof( counter ) = dofDim*dofMap[a]+component;
              BC_OP::ApplyBcValue( dof( counter ),
                                   blockSystem,
                                   blockID,
                                   rhsContribution( counter ),
                                   result[counter],
                                   rtTypes::value( field[a], component ));
              ++counter;
            }
            BC_OP::ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
          }
        }
      }
    } );
#endif
}



template< typename BC_OP, typename LAMBDA >
void
BoundaryConditionBase::
ApplyBoundaryConditionToSystem( set<localIndex> const & targetSet,
                                real64 const time,
                                dataRepository::ManagedGroup * dataGroup,
                                globalIndex_array const & dofMap,
                                integer const & dofDim,
                                systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                systemSolverInterface::BlockIDs const blockID,
                                LAMBDA && lambda ) const
{
  integer const component = GetComponent();
  string const functionName = getData<string>( viewKeyStruct::functionNameString );
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  integer const numBlocks = blockSystem->numBlocks();
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( blockID );

  globalIndex_array  dof( targetSet.size() );
  real64_array     rhsContribution( targetSet.size() );

  if( functionName.empty() )
  {

    integer counter=0;
    for( auto a : targetSet )
    {
      dof( counter ) = dofDim*dofMap[a]+component;
      BC_OP::ApplyBcValue( dof( counter ),
                           blockSystem,
                           blockID,
                           rhsContribution( counter ),
                           m_scale,
                           lambda( a ) );
      ++counter;
    }
    BC_OP::ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
  }
  else
  {
    FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>( functionName );

    GEOS_ASSERT( function != nullptr, "Function '" << functionName << "' not found" );

    if( function->isFunctionOfTime()==2 )
    {
      real64 value = m_scale * function->Evaluate( &time );
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofDim*integer_conversion<int>( dofMap[a] )+component;
        BC_OP::ApplyBcValue( dof( counter ),
                             blockSystem,
                             blockID,
                             rhsContribution( counter ),
                             value,
                             lambda( a ) );
        ++counter;
      }
      BC_OP::ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
    }
    else
    {
      real64_array result;
      result.resize( integer_conversion<localIndex>( targetSet.size()));
      function->Evaluate( dataGroup, time, targetSet, result );
      integer counter=0;
      for( auto a : targetSet )
      {
        dof( counter ) = dofDim*integer_conversion<int>( dofMap[a] )+component;
        BC_OP::ApplyBcValue( dof( counter ),
                             blockSystem,
                             blockID,
                             rhsContribution( counter ),
                             m_scale*result[counter],
                             lambda( a ) );
        ++counter;
      }
      BC_OP::ReplaceGlobalValues( rhs, counter, dof.data(), rhsContribution.data() );
    }
  }
}

//template<typename LAMBDA>
//void BoundaryConditionBase::ApplyBoundaryCondition( set<localIndex> const & targetSet,
//                                                    real64 const time,
//                                                    dataRepository::ManagedGroup * dataGroup,
//                                                    LAMBDA && lambda )
//{
//  integer const component = GetComponent();
//  string const functionName = getData<string>( viewKeyStruct::functionNameString );
//  NewFunctionManager * functionManager = NewFunctionManager::Instance();
//
//  if( functionName.empty())
//  {
//    real64 const value = m_scale;
//    integer counter = 0;
//    for( auto a : targetSet )
//    {
//      lambda( dataGroup, a, counter, value );
//      ++counter;
//    }
//  }
//  else
//  {
//    FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>( functionName );
//    if( function!=nullptr )
//    {
//      if( function->isFunctionOfTime() == 2 )
//      {
//        real64 const value = m_scale * function->Evaluate( &time );
//        integer counter = 0;
//        for( auto a : targetSet )
//        {
//          lambda( dataGroup, a, counter, value );
//          ++counter;
//        }
//      }
//      else
//      {
//        real64_array result;
//        result.resize( integer_conversion<localIndex>( targetSet.size()));
//        function->Evaluate( dataGroup, time, targetSet, result );
//        integer counter = 0;
//        for( auto a : targetSet )
//        {
//          real64 const value = m_scale * result[counter];
//          lambda( dataGroup, a, counter, value );
//          ++counter;
//        }
//      }
//    }
//  }
//}


}
#endif
