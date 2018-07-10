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

#include "BoundaryConditionBase.hpp"

namespace geosx
{
using namespace dataRepository;

BoundaryConditionBase::BoundaryConditionBase( string const & name, ManagedGroup * parent ):
  ManagedGroup(name,parent)
{
  RegisterViewWrapper( viewKeyStruct::setNamesString, &m_setNames, 0 );
  RegisterViewWrapper( viewKeyStruct::constitutivePathString, &m_constitutivePath, 0 );
  RegisterViewWrapper( viewKeyStruct::objectPathString, &m_objectPath, 0 );
  RegisterViewWrapper( viewKeyStruct::fieldNameString, &m_fieldName, 0 );
  RegisterViewWrapper( viewKeyStruct::componentString, &m_component, 0 );
  RegisterViewWrapper( viewKeyStruct::directionString, &m_direction, 0 );
  RegisterViewWrapper( viewKeyStruct::functionNameString, &m_functionName, 0 );
  RegisterViewWrapper( viewKeyStruct::bcApplicationTableNameString, &m_bcApplicationFunctionName, 0 );
  RegisterViewWrapper( viewKeyStruct::scaleString, &m_scale, 0 );
  RegisterViewWrapper( viewKeyStruct::initialConditionString, &m_initialCondition, 0 );
}


BoundaryConditionBase::~BoundaryConditionBase()
{}

BoundaryConditionBase::CatalogInterface::CatalogType& BoundaryConditionBase::GetCatalog()
{
  static BoundaryConditionBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void BoundaryConditionBase::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->AllocateChildNode( viewKeyStruct::initialConditionString,
                              viewKeyStruct::initialConditionString,
                              -1,
                              "integer",
                              "integer",
                              "BC is applied as an initial condition.",
                              "",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::setNamesString,
                              viewKeyStruct::setNamesString,
                              -1,
                              "string_array",
                              "string_array",
                              "Name of sets that boundary condition is applied to.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );


  docNode->AllocateChildNode( viewKeyStruct::constitutivePathString,
                              viewKeyStruct::constitutivePathString,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::objectPathString,
                              viewKeyStruct::objectPathString,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::fieldNameString,
                              viewKeyStruct::fieldNameString,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::dataTypeString,
                              viewKeyStruct::dataTypeString,
                              -1,
                              "string",
                              "string",
                              "Name of field that boundary condition is applied to.",
                              "",
                              "REQUIRED",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::componentString,
                              viewKeyStruct::componentString,
                              -1,
                              "integer",
                              "integer",
                              "Component of field (if tensor) to apply boundary condition to",
                              "",
                              "0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::directionString,
                              viewKeyStruct::directionString,
                              -1,
                              "R1Tensor",
                              "R1Tensor",
                              "Direction to apply boundary condition to",
                              "",
                              "{0,0,0}",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::functionNameString,
                              viewKeyStruct::functionNameString,
                              -1,
                              "string",
                              "string",
                              "Name of table that specifies variation of BC with time.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::bcApplicationTableNameString,
                              viewKeyStruct::bcApplicationTableNameString,
                              -1,
                              "string",
                              "string",
                              "Name of table that specifies the on/off application of the bc.",
                              "",
                              "",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( viewKeyStruct::scaleString,
                              viewKeyStruct::scaleString,
                              -1,
                              "real64",
                              "real64",
                              "Component of field (if tensor) to apply boundary condition to",
                              "",
                              "-1",
                              "",
                              0,
                              1,
                              0 );
}

void BoundaryConditionBase::ReadXML_PostProcess()
{
}

//real64 BoundaryConditionBase::GetValue( realT time ) const
//{
//
//  real64 rval = m_scale;
//  if (!(m_functionName.empty()))
//  {
//    array<real64> t(1);
//    t[0] = time;
//    real64 const tableval =
// TableManager::Instance().LookupTable<1>(m_functionName, t);
//    rval = m_scale * tableval;
//  }
////  else if (!(m_spaceFunctionName.empty()))
////  {
//////    rval = m_scale * (*m_function)(time);
////  }
//  return rval;
//}
//
//void BoundaryConditionBase::ApplyBounaryConditionDefaultMethod( lSet const &
// set,
//                                                                real64 const
// time,
//                                                                array<R1Tensor>
// const & X,
//                                                                array<R1Tensor>
// & field )
//{
//
//  integer const component = GetComponent();
//  string const functionName =
// getData<string>(dataRepository::keys::functionName);
//  NewFunctionManager * functionManager = NewFunctionManager::Instance();
//
//  FunctionBase const * const function  =
// functionManager->GetGroup<FunctionBase>(functionName);
//
//  if( function!=nullptr )
//  {
//    real64 const factor = m_scale * ( timeFunction->Evaluate( &time ) );
//    for( auto a : set )
//    {
//      field[a][component] = factor;
//    }
//  }
//  else
//  {
//    for( auto a : set )
//    {
//      field[a][component] = m_scale;
//    }
//  }
//}
//
//void BoundaryConditionBase::ApplyBounaryConditionDefaultMethod( lSet const &
// set,
//                                                                real64 const
// time,
//                                                                ManagedGroup *
// dataGroup,
//                                                                string const &
// fieldName ) const
//{
//
//  integer const component = GetComponent();
//  string const functionName =
// getData<string>(dataRepository::keys::functionName);
//  NewFunctionManager * functionManager = NewFunctionManager::Instance();
//
//  ViewWrapperBase * vw = dataGroup->getWrapperBase( fieldName );
//  std::type_index typeIndex = std::type_index(vw->get_typeid());
//
//  rtTypes::ApplyArrayTypeLambda1( rtTypes::typeID(typeIndex) , [&]( auto type
// ) -> void
//  {
//    using fieldType = decltype(type);
//    ViewWrapper<fieldType> & view = dynamic_cast< ViewWrapper<fieldType> &
// >(*vw);
//    view_rtype<fieldType> field = view.data();
//    if( functionName.empty() )
//    {
//      for( auto a : set )
//      {
//        rtTypes::equate( field[a], component, m_scale );
//      }
//    }
//    else
//    {
//      FunctionBase const * const function  =
// functionManager->GetGroup<FunctionBase>(functionName);
//      if( function!=nullptr)
//      {
//        if( function->isFunctionOfTime()==2 )
//        {
//          real64 value = m_scale * function->Evaluate( &time );
//          for( auto a : set )
//          {
//            rtTypes::equate( field[a], component, value );
//          }
//        }
//        else
//        {
//          real64_array result(set.size());
//          function->Evaluate( dataGroup, time, set, result );
//          integer count=0;
//          for( auto a : set )
//          {
//            rtTypes::equate( field[a], component, result[count] );
//            ++count;
//          }
//        }
//      }
//    }
//  });
//}



}
