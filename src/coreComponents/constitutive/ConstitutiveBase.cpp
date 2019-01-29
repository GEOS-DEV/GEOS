/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
  ManagedGroup( name, parent )
{}

ConstitutiveBase::~ConstitutiveBase()
{}



ConstitutiveBase::CatalogInterface::CatalogType& ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ConstitutiveBase::AllocateConstitutiveData( dataRepository::ManagedGroup * const parent,
                                                 localIndex const )
{
  m_constitutiveDataGroup = parent;

  for( auto & group : this->GetSubGroups() )
  {
    for( auto & wrapper : group.second->wrappers() )
    {
      string const wrapperName = wrapper.first;
      std::unique_ptr<ViewWrapperBase> newWrapper = wrapper.second->clone( wrapperName, parent );
      parent->RegisterViewWrapper( makeFieldName(this->getName(), wrapperName), newWrapper.release() );
    }
  }

  for( auto & wrapper : this->wrappers() )
  {
    string const wrapperName = wrapper.first;
    std::unique_ptr<ViewWrapperBase> newWrapper = wrapper.second->clone( wrapperName, parent );
    parent->RegisterViewWrapper( makeFieldName(this->getName(), wrapperName), newWrapper.release() );
  }

}

void ConstitutiveBase::resize( localIndex newsize )
{
  ManagedGroup::resize( newsize );
}



}
} /* namespace geosx */
