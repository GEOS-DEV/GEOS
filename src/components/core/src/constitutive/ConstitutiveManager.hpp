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
 * ConstitutiveManager.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#include "dataRepository/ManagedGroup.hpp"
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{


//template< typename T, bool HASPOINTERTYPE = std::is_pointer<T>::value >
//struct ConstitutiveWrapper
//{
//  T m_object;
//};
//
//template< typename T >
//struct ConstitutiveWrapper<T,false>
//{
//  T & m_object;
//};



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
  array< T * > GetParameterData( string const & name );

  template< typename T >
  array< T * > GetData( string const & name );

  template< typename T >
  array< T const * > GetData( string const & name ) const;


};



template< typename T >
array< T * > ConstitutiveManager::GetParameterData( string const & name )
{
  array< T * > rval;
//  string key = dataRepository::keys::parameterData;
//  this->forSubGroups( [this,&name, &rval, &key]( ManagedGroup * material ) ->
// void
//  {
//    dataRepository::view_rtype<T> temp0 =
// material->GetGroup(key)->getData<T>(name);
//    ConstitutiveWrapper< dataRepository::view_rtype<T> > temp( temp0 );
//    rval.push_back( std::move(temp) );
//  });

  for( auto& subGroupIter : this->GetSubGroups() )
  {
    T * temp{subGroupIter.second->GetGroup(dataRepository::keys::parameterData)->getData<T>(name)};
    rval.push_back( std::move(temp) );
  }
  return rval;
}

template< typename T >
array< T * >
ConstitutiveManager::GetData( string const & name )
{
  array< T * > rval;
  rval.resize( this->GetSubGroups().size() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase * const material = GetGroup<ConstitutiveBase>(a);
    if( material->GetStateData()->hasView(name) )
    {
      rval[a] = &(material->GetStateData()->getReference<T>(name));
    }
    else
    {
      rval[a] = nullptr;
    }
  }
  return rval;
}

template< typename T >
array< T const * >
ConstitutiveManager::GetData( string const & name ) const
{
  array< T const * > rval;
  rval.resize( this->GetSubGroups().size() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase const * const material = GetGroup<ConstitutiveBase>(a);
    if( material->GetStateData()->hasView(name) )
    {
      rval[a] = &(material->GetStateData()->getReference<T>(name));
    }
    else
    {
      rval[a] = nullptr;
    }
  }
  return rval;
}

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
