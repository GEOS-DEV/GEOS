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
 * @file ConstitutiveManager.hpp
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

/**
 * @typedef This is an alias for array of reference wrappers which are used to provide a
 * multidimensional array like interface via the operator[].
 */
template< typename VIEWTYPE >
using ViewAccessor = array1d < VIEWTYPE >;

/**
 * @class ConstitutiveManager
 * @brief Class to manage the allocation and access to constitutive relations
 */
class ConstitutiveManager : public dataRepository::ManagedGroup
{
public:
  ConstitutiveManager() = delete;

  ConstitutiveManager( string const & name,
                       ManagedGroup * const parent );

  virtual void CreateChild( string const & childKey, string const & childName ) override final;

  ConstitutiveBase *
  HangConstitutiveRelation( string const & constitutiveRelationInstanceName,
                            dataRepository::ManagedGroup * const parent,
                            localIndex const numConstitutivePointsPerParentIndex ) const;

  ~ConstitutiveManager() override;

  template< typename T = ConstitutiveBase >
  T const * GetConstitituveRelation( string const & constitutiveRelationInstanceName ) const
  {
    return this->GetGroup<T>( constitutiveRelationInstanceName );
  }

  template< typename T = ConstitutiveBase >
  T * GetConstitituveRelation( string const & constitutiveRelationInstanceName )
  {
    return this->GetGroup<T>( constitutiveRelationInstanceName );
  }

  template< typename T = ConstitutiveBase >
  T const * GetConstitituveRelation( localIndex const index ) const
  {
    return this->GetGroup<T>( index );
  }

  template< typename T = ConstitutiveBase >
  T * GetConstitituveRelation( localIndex const index )
  {
    return this->GetGroup<T>( index );
  }

  // template< typename T >
  // ViewAccessor< T >
  // GetConstitutiveData( string const & name,
  //                      dataRepository::ManagedGroup * const relationGroup );

  template< typename T >
  ViewAccessor< T >
  GetConstitutiveData( string const & name,
                       dataRepository::ManagedGroup const * const relationGroup ) const;


  struct groupKeyStruct
  {
    static constexpr auto constitutiveModelsString = "ConstitutiveModels";
  } m_ConstitutiveManagerGroupKeys;


};



template< typename T >
ViewAccessor< T >
ConstitutiveManager::GetConstitutiveData( string const & name,
                                          dataRepository::ManagedGroup const * const relationGroup ) const
{
  ViewAccessor< T const > rval( relationGroup->numSubGroups() );

  rval.resize( relationGroup->numSubGroups() );
  for( localIndex a=0 ; a<this->GetSubGroups().size() ; ++a )
  {
    ConstitutiveBase const * const material = relationGroup->GetGroup<ConstitutiveBase>( a );
    if( material->hasView( name ) )
    {
      rval[a] = material->getReference<T>( name );
    }
  }
  return rval;
}

// template< typename T >
// ViewAccessor< T >
// ConstitutiveManager::GetConstitutiveData( string const & name,
//                                           dataRepository::ManagedGroup * const relationGroup )
// {
//   return const_cast< ViewAccessor<T> >(const_cast<ConstitutiveManager const *>(this)->
//                                        GetConstitutiveData<T>( name, relationGroup ) );
// }

}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
