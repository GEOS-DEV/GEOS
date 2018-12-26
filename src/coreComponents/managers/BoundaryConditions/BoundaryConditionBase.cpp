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
  ManagedGroup( name, parent )
//  m_setNames(),
//  m_objectPath(),
//  m_fieldName(),
//  m_component(-1),
//  m_direction(-1),
//  m_initialCondition(0),
//  m_functionName(),
//  m_scale(0.0),
//  m_beginTime(0.0),
//  m_endTime(1e9),
//  m_bcApplicationFunctionName()
{
  RegisterViewWrapper( viewKeyStruct::setNamesString, &m_setNames, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of sets that boundary condition is applied to.");

  RegisterViewWrapper( viewKeyStruct::objectPathString, &m_objectPath, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of field that boundary condition is applied to.");

  RegisterViewWrapper( viewKeyStruct::fieldNameString, &m_fieldName, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of field that boundary condition is applied to.");

  RegisterViewWrapper( viewKeyStruct::componentString, &m_component, false )->
    setDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Component of field (if tensor) to apply boundary condition to");

  RegisterViewWrapper( viewKeyStruct::directionString, &m_direction, false )->
    setDefaultValue(R1Tensor())->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Direction to apply boundary condition to");

  RegisterViewWrapper( viewKeyStruct::functionNameString, &m_functionName, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of function that specifies variation of BC with time.");

  RegisterViewWrapper( viewKeyStruct::bcApplicationTableNameString, &m_bcApplicationFunctionName, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of table that specifies the on/off application of the bc.");

  RegisterViewWrapper( viewKeyStruct::scaleString, &m_scale, false )->
    setDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Scale multiplier of the BC.");

  RegisterViewWrapper( viewKeyStruct::initialConditionString, &m_initialCondition, false )->
    setDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("BC is applied as an initial condition.");

  RegisterViewWrapper( viewKeyStruct::beginTimeString, &m_beginTime, false )->
    setDefaultValue(-1.0e99)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Time at which BC will start being applied");

  RegisterViewWrapper( viewKeyStruct::endTimeString, &m_endTime, false )->
    setDefaultValue(1.0e99)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Time at which bc will stop being applied");
}


BoundaryConditionBase::~BoundaryConditionBase()
{}

BoundaryConditionBase::CatalogInterface::CatalogType&
BoundaryConditionBase::GetCatalog()
{
  static BoundaryConditionBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void BoundaryConditionBase::ReadXML_PostProcess()
{}



}
