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

/*
 * ConstitutiveManager.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{

template< typename VIEWTYPE >
using ViewAccessor = array < ReferenceWrapper< VIEWTYPE > > ;


class ConstitutiveManager : public dataRepository::ManagedGroup
{
public:
  ConstitutiveManager() = delete;

  ConstitutiveManager( std::string const & name,
                       ManagedGroup * const parent );

  void FillDocumentationNode() override final;
  virtual void CreateChild( string const & childKey, string const & childName ) override final;

  ~ConstitutiveManager() override;

//  using constitutiveMaps = std::pair< array<ManagedGroup const *> ,
// map<string,integer> > ;
//  constitutiveMaps & GetMaps( integer const reinit ) const;


  template< typename T >
  ViewAccessor< T > GetParameterData( string const & name );

  template< typename T >
  ViewAccessor< T const > GetParameterData( string const & name ) const;

  template< typename T >
  ViewAccessor< T > GetStateData( string const & name );

  template< typename T >
  ViewAccessor< T const > GetStateData( string const & name ) const;


};



template< typename T >
ViewAccessor< T > ConstitutiveManager::GetParameterData( string const & name )
{
  ViewAccessor< T > rval;

  rval.resize( this->GetSubGroups().size() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase * const material = GetGroup<ConstitutiveBase>(a);
    if( material->GetParameterData()->hasView(name) )
    {
      rval[a].set(material->GetParameterData()->getReference<T>(name));
    }
    else
    {
      rval[a].set(nullptr);
    }
  }
  return rval;
}

template< typename T >
ViewAccessor< T const > ConstitutiveManager::GetParameterData( string const & name ) const
{
  ViewAccessor< T const > rval;

  rval.resize( this->GetSubGroups().size() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase const * const material = GetGroup<ConstitutiveBase>(a);
    if( material->GetParameterData()->hasView(name) )
    {
      rval[a].set(material->GetParameterData()->getReference<T>(name));
    }
    else
    {
      rval[a].set(nullptr);
    }
  }
  return rval;
}


template< typename T >
ViewAccessor< T >
ConstitutiveManager::GetStateData( string const & name )
{
  ViewAccessor< T > rval;
  rval.resize( this->GetSubGroups().size() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase * const material = GetGroup<ConstitutiveBase>(a);
    if( material->GetStateData()->hasView(name) )
    {
      rval[a].set(material->GetStateData()->getReference<T>(name));
    }
    else
    {
      rval[a].set(nullptr);
    }
  }
  return rval;
}

template< typename T >
ViewAccessor< T const >
ConstitutiveManager::GetStateData( string const & name ) const
{
  ViewAccessor< T const > rval;
  rval.resize( this->GetSubGroups().size() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase const * const material = GetGroup<ConstitutiveBase>(a);
    if( material->GetStateData()->hasView(name) )
    {
      rval[a].set(material->GetStateData()->getReference<T>(name));
    }
    else
    {
      rval[a].set(nullptr);
    }
  }
  return rval;
}

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
