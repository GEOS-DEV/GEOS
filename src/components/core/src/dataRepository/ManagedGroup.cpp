/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "ManagedGroup.hpp"

#include "dataRepository/SidreWrapper.hpp"

namespace geosx
{
namespace dataRepository
{

ManagedGroup::ManagedGroup( std::string const & name,
                                      ManagedGroup * const parent ) :
  m_keyLookup(),
  m_wrappers(),
  m_parent(parent),
  m_subObjectManagers(),
  m_sidreGroup(nullptr)
{
  asctoolkit::sidre::DataGroup * sidreParent = nullptr;
  if( m_parent==nullptr )
  {
    sidreParent = SidreWrapper::dataStore().getRoot();
  }
  else
  {
    sidreParent = parent->m_sidreGroup;
  }

  if( sidreParent->hasGroup(name) )
  {
    m_sidreGroup = sidreParent->getGroup(name);
  }
  else
  {
    m_sidreGroup = sidreParent->createGroup(name);
  }

  *(RegisterViewWrapper<localIndex>( "size" ).data()) = 0;
  RegisterViewWrapper<std::string>( "name" ).reference() = name;
  RegisterViewWrapper<std::string>( "path" );


}

ManagedGroup::~ManagedGroup()
{
  // TODO Auto-generated destructor stub
}

//DataObjectManager::DataObjectManager( DataObjectManager const & source ):
//    m_size( source.m_size ),
//    m_name( source.m_name ),
//    m_path( source.m_path ),
//    m_keyLookup( source.m_keyLookup ),
//    m_dataObjects( source.m_dataObjects ),
//    m_parent( source.m_parent )
//{}

ManagedGroup::ManagedGroup( ManagedGroup&& source ) :
  m_keyLookup( std::move(source.m_keyLookup) ),
  m_wrappers( std::move(source.m_wrappers) ),
  m_parent( std::move(source.m_parent) )
{}

ManagedGroup::CatalogInterface::CatalogType& ManagedGroup::GetCatalog()
{
  static ManagedGroup::CatalogInterface::CatalogType catalog;
  return catalog;
}

ViewWrapperBase& ManagedGroup::RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type )
{
  return *( rtTypes::ApplyTypeLambda( type,
                                      [this, &name]( auto a ) -> ViewWrapperBase*
      {
        return &( this->RegisterViewWrapper<decltype(a)>(name) );
      } ) );
}

void ManagedGroup::resize( std::size_t const newsize )
{
  for( auto&& i : this->m_wrappers )
  {
    i->resize(newsize);
  }
  *(this->getWrapper<std_size_t>( keys::Size ).data())=newsize;
}



}
} /* namespace ODS */
