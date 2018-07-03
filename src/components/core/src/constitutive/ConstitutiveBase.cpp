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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/**
 * @file ConstitutiveBase.cpp
 */



#include "ConstitutiveBase.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ConstitutiveBase::ConstitutiveBase( std::string const & name,
                                    ManagedGroup * const parent ):
  ManagedGroup(name,parent),
  m_parameterData(groupKeys().ParameterData.Key(),this),
  m_stateData(groupKeys().StateData.Key(),this)
{
  RegisterGroup(groupKeys().ParameterData.Key(), &m_parameterData, 0 );
  RegisterGroup(groupKeys().StateData.Key(), &m_stateData, 0);
}

ConstitutiveBase::~ConstitutiveBase()
{}



ConstitutiveBase::CatalogInterface::CatalogType& ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ConstitutiveBase::resize( localIndex newsize )
{
  ManagedGroup::resize(newsize);
  GetParameterData()->resize(newsize);
  GetStateData()->resize(newsize);
}

void ConstitutiveBase::SetVariableParameters()
{
  for( auto & viewBase : GetParameterData()->wrappers() )
  {
    viewBase.second->setSizedFromParent(1);
  }
}


}
} /* namespace geosx */
