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



#ifndef BOUNDARYCONDITIONBASE_H
#define BOUNDARYCONDITIONBASE_H

#include "common/DataTypes.hpp"
#include "dataRepository/ManagedGroup.hpp"
//#include "managers/TableManager.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{
class Function;


class BoundaryConditionBase : public dataRepository::ManagedGroup
{
public:

  using CatalogInterface = cxx_utilities::CatalogInterface< BoundaryConditionBase, string const &, dataRepository::ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  BoundaryConditionBase( string const & name, dataRepository::ManagedGroup * parent );

  virtual ~BoundaryConditionBase() override;



  void FillDocumentationNode() override;

  void ReadXML_PostProcess() override final;



//  real64 GetValue( realT time ) const;


//  template< typename T >
//  void ApplyBounaryConditionDefaultMethod( lSet const & set,
//                                           real64 const time,
//                                           array<R1Tensor> const & X,
//                                           array<T> & field );

//  void ApplyBounaryConditionDefaultMethod( lSet const & set,
//                                           real64 const time,
//                                           array<R1Tensor> const & X,
//                                           array<R1Tensor> & field );

  template< typename OPERATION >
  void ApplyBounaryConditionDefaultMethod( lSet const & set,
                                           real64 const time,
                                           dataRepository::ManagedGroup * dataGroup,
                                           string const & fieldname ) const;

  // calls user-provided lambda to apply computed boundary value
  template<typename LAMBDA>
  void ApplyBoundaryCondition(lSet const & set,
                              real64 const time,
                              dataRepository::ManagedGroup * dataGroup,
                              LAMBDA && lambda);

  template< int OPERATION >
  void ApplyDirichletBounaryConditionDefaultMethod( lSet const & set,
                                                    real64 const time,
                                                    dataRepository::ManagedGroup * dataGroup,
                                                    string const & fieldName,
                                                    string const & dofMapName,
                                                    integer const & dofDim,
                                                    systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                    systemSolverInterface::BlockIDs const blockID ) const;


  template< int OPERATION, typename LAMBDA >
  void
  ApplyDirichletBounaryConditionDefaultMethod( lSet const & set,
                                               real64 const time,
                                               dataRepository::ManagedGroup * dataGroup,
                                               globalIndex_array const & dofMap,
                                               integer const & dofDim,
                                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                               systemSolverInterface::BlockIDs const blockID,
                                               LAMBDA && lambda ) const;

  template< int OPERATION >
  inline void ApplyBounaryConditionDefaultMethodPoint( globalIndex const dof,
                                                       systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                       systemSolverInterface::BlockIDs const blockID,
                                                       real64 & rhs,
                                                       real64 const & bcValue,
                                                       real64 const fieldValue ) const;


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

  } viewKeys;

  struct groupKeyStruct
  {
  } groupKeys;



  string const & GetFunctionName() const
  {
    return m_functionName;
  }

  virtual const string& GetConstitutivePath() const
  {
    return m_constitutivePath;
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

  virtual const R1Tensor& GetDirection(realT time)
  {
    return m_direction;
  }

  real64 GetStartTime() const
  {
    return -1;
  }

  real64 GetEndTime() const
  {
    return 1.0e9;
  }

  string_array const & GetSetNames() const
  {
    return m_setNames;
  }

  int initialCondition() const
  {
    return m_initialCondition;
  }

private:

  string_array m_setNames; // sets the boundary condition is applied to

  string m_constitutivePath;
  string m_objectPath;

  string m_fieldName;    // the name of the field the boundary condition is
                         // applied to or a description of the boundary
                         // condition.


  string m_dataType;
  // TODO get rid of components. Replace with direction only.

  int m_component;       // the component the boundary condition acts on (-ve
                         // indicates that direction should be used).
  R1Tensor m_direction;  // the direction the boundary condition acts in.

  int m_initialCondition;

  string m_functionName;
  string m_bcApplicationFunctionName;

  real64 m_scale;


};


template< typename OPERATION >
void BoundaryConditionBase::ApplyBounaryConditionDefaultMethod( lSet const & set,
                                                                real64 const time,
                                                                ManagedGroup * dataGroup,
                                                                string const & fieldName ) const
{

  integer const component = GetComponent();
  string const functionName = getData<string>(viewKeyStruct::functionNameString);
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  dataRepository::ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index(vw->get_typeid());

  rtTypes::ApplyArrayTypeLambda2( rtTypes::typeID(typeIndex), [&]( auto type, auto baseType ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::ViewWrapper<fieldType> & view = dynamic_cast< dataRepository::ViewWrapper<fieldType> & >(*vw);
      dataRepository::view_rtype<fieldType> field = view.data();
      if( functionName.empty() )
      {
        for( auto a : set )
        {
          OPERATION::f( field[a], component, (m_scale) );
//        OPERATION::f( field[a], component,
// static_cast<decltype(baseType)>(m_scale) );
        }
      }
      else
      {
        FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
        if( function!=nullptr)
        {
          if( function->isFunctionOfTime()==2 )
          {
            real64 value = m_scale * function->Evaluate( &time );
            for( auto a : set )
            {
              OPERATION::f( field[a], component, (value) );
            }
          }
          else
          {
            real64_array result(static_cast<localIndex>(set.size()));
            function->Evaluate( dataGroup, time, set, result );
            integer count=0;
            for( auto a : set )
            {
              OPERATION::f( field[a], component, (result[count]) );
              ++count;
            }
          }
        }
      }
    });
}



template<>
inline void BoundaryConditionBase::ApplyBounaryConditionDefaultMethodPoint<0>( globalIndex const dof,
                                                                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                                               systemSolverInterface::BlockIDs const blockID,
                                                                               real64 & rhs,
                                                                               real64 const & bcValue,
                                                                               real64 const fieldValue ) const
{

  if( true )//node_is_ghost[*nd] < 0 )
  {
    real64 LARGE = blockSystem->ClearSystemRow( blockID, dof, 1.0 );
    rhs = -LARGE*( bcValue - fieldValue );
  }
  else
  {
    blockSystem->ClearSystemRow( blockID, dof, 0.0 );
    rhs = 0.0;
  }
}


template<>
inline void BoundaryConditionBase::ApplyBounaryConditionDefaultMethodPoint<1>( globalIndex const dof,
                                                                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                                               systemSolverInterface::BlockIDs const blockID,
                                                                               real64 & rhs,
                                                                               real64 const & bcValue,
                                                                               real64 const fieldValue ) const
{
  if( true )//node_is_ghost[*nd] < 0 )
  {
    rhs += bcValue;
  }
}





template< int OPERATION >
void BoundaryConditionBase::ApplyDirichletBounaryConditionDefaultMethod( lSet const & set,
                                                                         real64 const time,
                                                                         dataRepository::ManagedGroup * dataGroup,
                                                                         string const & fieldName,
                                                                         string const & dofMapName,
                                                                         integer const & dofDim,
                                                                         systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                                         systemSolverInterface::BlockIDs const blockID ) const
{
  integer const component = GetComponent();
  string const functionName = getData<string>(viewKeyStruct::functionNameString);
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  dataRepository::ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
  std::type_index typeIndex = std::type_index(vw->get_typeid());

  integer const numBlocks = blockSystem->numBlocks();
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( blockID );

  Epetra_IntSerialDenseVector  node_dof( integer_conversion<int>( set.size() ) );
  Epetra_SerialDenseVector     node_rhs( integer_conversion<int>( set.size() ) );


  dataRepository::view_rtype_const<localIndex_array> dofMap = dataGroup->getData<localIndex_array>(dofMapName);


  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID(typeIndex), [&]( auto type ) -> void
    {
      using fieldType = decltype(type);
      dataRepository::ViewWrapper<fieldType> & view = dynamic_cast< dataRepository::ViewWrapper<fieldType> & >(*vw);
      dataRepository::view_rtype<fieldType> field = view.data();
      if( functionName.empty() )
      {

        integer counter=0;
        for( auto a : set )
        {
          node_dof(counter) = dofDim*integer_conversion<int>(dofMap[a])+component;
          this->ApplyBounaryConditionDefaultMethodPoint<OPERATION>( node_dof(counter),
                                                                    blockSystem,
                                                                    blockID,
                                                                    node_rhs(counter),
                                                                    m_scale,
                                                                    static_cast<real64>(rtTypes::value(field[a],component)));
          ++counter;
        }
        if( OPERATION==0 )
        {
          rhs->ReplaceGlobalValues(node_dof, node_rhs);
        }
        else if( OPERATION==1 )
        {
          rhs->SumIntoGlobalValues(node_dof, node_rhs);
        }
      }
      else
      {
        FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
        if( function!=nullptr)
        {
          if( function->isFunctionOfTime()==2 )
          {
            real64 value = m_scale * function->Evaluate( &time );
            integer counter=0;
            for( auto a : set )
            {
              node_dof(counter) = dofDim*integer_conversion<int>(dofMap[a])+component;
              this->ApplyBounaryConditionDefaultMethodPoint<OPERATION>( node_dof(counter),
                                                                        blockSystem,
                                                                        blockID,
                                                                        node_rhs(counter),
                                                                        value,
                                                                        rtTypes::value(field[a],component));
              ++counter;
            }
            if( OPERATION==0 )
            {
              rhs->ReplaceGlobalValues(node_dof, node_rhs);
            }
            else if( OPERATION==1 )
            {
              rhs->SumIntoGlobalValues(node_dof, node_rhs);
            }
          }
          else
          {
            real64_array result;
            result.resize( integer_conversion<localIndex>(set.size()));
            function->Evaluate( dataGroup, time, set, result );
            integer counter=0;
            for( auto a : set )
            {
              node_dof(counter) = dofDim*integer_conversion<int>(dofMap[a])+component;
              this->ApplyBounaryConditionDefaultMethodPoint<OPERATION>( node_dof(counter),
                                                                        blockSystem,
                                                                        blockID,
                                                                        node_rhs(counter),
                                                                        result[counter],
                                                                        rtTypes::value(field[a],component));
              ++counter;
            }
            if( OPERATION==0 )
            {
              rhs->ReplaceGlobalValues(node_dof, node_rhs);
            }
            else if( OPERATION==1 )
            {
              rhs->SumIntoGlobalValues(node_dof, node_rhs);
            }

          }
        }
      }
    });
}




template< int OPERATION, typename LAMBDA >
void
BoundaryConditionBase::
ApplyDirichletBounaryConditionDefaultMethod( lSet const & set,
                                             real64 const time,
                                             dataRepository::ManagedGroup * dataGroup,
                                             globalIndex_array const & dofMap,
                                             integer const & dofDim,
                                             systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                             systemSolverInterface::BlockIDs const blockID,
                                             LAMBDA && lambda ) const
{
  integer const component = GetComponent();
  string const functionName = getData<string>(viewKeyStruct::functionNameString);
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  integer const numBlocks = blockSystem->numBlocks();
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( blockID );

  globalIndex_array  node_dof( set.size() );
  real64_array     node_rhs( set.size() );

  if( functionName.empty() )
  {

    integer counter=0;
    for( auto a : set )
    {
      node_dof(counter) = dofDim*dofMap[a]+component;
      this->ApplyBounaryConditionDefaultMethodPoint<OPERATION>( node_dof(counter),
                                                                blockSystem,
                                                                blockID,
                                                                node_rhs(counter),
                                                                m_scale,
                                                                lambda(a) );
      ++counter;
    }
    if( OPERATION==0 )
    {
      rhs->ReplaceGlobalValues( node_dof.size(), node_dof.data(), node_rhs.data() );
    }
    else if( OPERATION==1 )
    {
      rhs->SumIntoGlobalValues( node_dof.size(), node_dof.data(), node_rhs.data() );
    }
  }
  else
  {
    FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
    if( function!=nullptr)
    {
      if( function->isFunctionOfTime()==2 )
      {
        real64 value = m_scale * function->Evaluate( &time );
        integer counter=0;
        for( auto a : set )
        {
          node_dof(counter) = dofDim*integer_conversion<int>(dofMap[a])+component;
          this->ApplyBounaryConditionDefaultMethodPoint<OPERATION>( node_dof(counter),
                                                                    blockSystem,
                                                                    blockID,
                                                                    node_rhs(counter),
                                                                    value,
                                                                    lambda(a) );
          ++counter;
        }
        if( OPERATION==0 )
        {
          rhs->ReplaceGlobalValues( node_dof.size(), node_dof.data(), node_rhs.data() );
        }
        else if( OPERATION==1 )
        {
          rhs->SumIntoGlobalValues( node_dof.size(), node_dof.data(), node_rhs.data() );
        }
      }
      else
      {
        real64_array result;
        result.resize( integer_conversion<localIndex>(set.size()));
        function->Evaluate( dataGroup, time, set, result );
        integer counter=0;
        for( auto a : set )
        {
          node_dof(counter) = dofDim*integer_conversion<int>(dofMap[a])+component;
          this->ApplyBounaryConditionDefaultMethodPoint<OPERATION>( node_dof(counter),
                                                                    blockSystem,
                                                                    blockID,
                                                                    node_rhs(counter),
                                                                    m_scale*result[counter],
                                                                    lambda(a) );
          ++counter;
        }
        if( OPERATION==0 )
        {
          rhs->ReplaceGlobalValues( node_dof.size(), node_dof.data(), node_rhs.data() );
        }
        else if( OPERATION==1 )
        {
          rhs->SumIntoGlobalValues( node_dof.size(), node_dof.data(), node_rhs.data() );
        }

      }
    }
  }
}

template<typename LAMBDA>
void BoundaryConditionBase::ApplyBoundaryCondition(lSet const & set,
                                                   real64 const time,
                                                   dataRepository::ManagedGroup * dataGroup,
                                                   LAMBDA && lambda)
{
  integer const component = GetComponent();
  string const functionName = getData<string>(viewKeyStruct::functionNameString);
  NewFunctionManager * functionManager = NewFunctionManager::Instance();

  if (functionName.empty())
  {
    real64 const value = m_scale;
    integer counter = 0;
    for (auto a : set)
    {
      lambda(dataGroup, a, counter, value);
      ++counter;
    }
  }
  else
  {
    FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
    if (function!=nullptr)
    {
      if (function->isFunctionOfTime() == 2)
      {
        real64 const value = m_scale * function->Evaluate(&time);
        integer counter = 0;
        for (auto a : set)
        {
          lambda(dataGroup, a, counter, value);
          ++counter;
        }
      }
      else
      {
        real64_array result;
        result.resize(integer_conversion<localIndex>(set.size()));
        function->Evaluate(dataGroup, time, set, result);
        integer counter = 0;
        for (auto a : set)
        {
          real64 const value = m_scale * result[counter];
          lambda(dataGroup, a, counter, value);
          ++counter;
        }
      }
    }
  }
}


}
#endif
