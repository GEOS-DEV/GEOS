/*
 * ConstitutiveManager.hpp
 *
 *  Created on: Aug 4, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_
#include "../dataRepository/ManagedGroup.hpp"
#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{


template< typename T , bool HASPOINTERTYPE = std::is_pointer<T>::value >
struct ConstitutiveWrapper
{
  T m_object;
};

template< typename T >
struct ConstitutiveWrapper<T,false>
{
  T & m_object;
};




class ConstitutiveManager : public dataRepository::ManagedGroup
{
public:
  ConstitutiveManager() = delete;

  ConstitutiveManager( std::string const & name,
                       ManagedGroup * const parent );

  void ReadXMLsub( pugi::xml_node const & targetNode );

  ~ConstitutiveManager();

  using constitutiveMaps = std::pair< array<ManagedGroup *> , map<string,int32> > ;
  constitutiveMaps & GetMaps( int32 const reinit ) const;


  template< typename T >
  array< ConstitutiveWrapper< dataRepository::view_rtype<T> > > GetParameterData( string const & name );

  template< typename T >
  array< ConstitutiveWrapper< dataRepository::view_rtype<T> > > GetData( string const & name );

  template< typename T >
  array< ConstitutiveWrapper< dataRepository::view_rtype_const<T> > > GetData( string const & name ) const;

};



template< typename T >
array< ConstitutiveWrapper< dataRepository::view_rtype<T> > > ConstitutiveManager::GetParameterData( string const & name )
{
  array< ConstitutiveWrapper< dataRepository::view_rtype<T> > > rval;
//  string key = dataRepository::keys::parameterData;
//  this->forSubGroups( [this,&name, &rval, &key]( ManagedGroup & material ) -> void
//  {
//    dataRepository::view_rtype<T> temp0 = material.GetGroup(key).getData<T>(name);
//    ConstitutiveWrapper< dataRepository::view_rtype<T> > temp( temp0 );
//    rval.push_back( std::move(temp) );
//  });

  for( auto& subGroupIter : this->GetSubGroups() )
  {
    ConstitutiveWrapper< dataRepository::view_rtype<T> > temp( {subGroupIter.second->GetGroup(dataRepository::keys::parameterData).getData<T>(name)} );
    rval.push_back( std::move(temp) );
  }
  return rval;
}

template< typename T >
array< ConstitutiveWrapper< dataRepository::view_rtype<T> > > ConstitutiveManager::GetData( string const & name )
{
  array< ConstitutiveWrapper< dataRepository::view_rtype<T> > > rval;
  this->forSubGroups( [this,&name, &rval]( ManagedGroup & material ) -> void
  {
    rval.push_back( material.getData<T>(name) );
  });
  return rval;
}


}
} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_CONSTITUTIVEMANAGER_HPP_ */
