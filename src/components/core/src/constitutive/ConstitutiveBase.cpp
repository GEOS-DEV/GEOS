// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/*
 * ConstitutiveBase.cpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#include "ConstitutiveBase.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ConstitutiveBase::ConstitutiveBase( std::string const & name,
                                    ManagedGroup * const parent ):
  ManagedGroup(name,parent)
{
  this->RegisterGroup<ManagedGroup>(groupKeys().ParameterData);
  this->RegisterGroup<ManagedGroup>(groupKeys().StateData);
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
