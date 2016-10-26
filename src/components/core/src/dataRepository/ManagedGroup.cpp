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
  ManagedGroup( name, parent, nullptr )
{}


ManagedGroup::ManagedGroup( std::string const & name,
                            ManagedGroup * const parent,
                            cxx_utilities::DocumentationNode * docNode ) :
  m_docNode(docNode),
  m_keyLookup(),
  m_wrappers(),
  m_parent(parent),
  m_subGroups(),
  m_sidreGroup(nullptr)
{

  // SIDRE interaction
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



  // Setup DocumentationNode
  if( parent != nullptr )
  {
    if( parent->m_docNode != nullptr )
    {
      m_docNode = parent->m_docNode->AllocateChildNode( name,
                                                        name,
                                                        0,
                                                        "ManagedGroup",
                                                        "",
                                                        "ManagedGroup",
                                                        "ManagedGroup",
                                                        "",
                                                        parent->getName(),
                                                        0,
                                                        0 ) ;

    }
  }


  m_docNode->AllocateChildNode( "size",
                                "size",
                                -1,
                                "localIndex",
                                "",
                                "size of group",
                                "Number of entries in this group.",
                                "0",
                                "",
                                0,
                                0 );

  m_docNode->AllocateChildNode( "name",
                                "name",
                                -1,
                                "string",
                                "string",
                                "name of group",
                                "name of group.",
                                name,
                                "",
                                0,
                                0 );

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

void ManagedGroup::resize( localIndex const newsize )
{
  for( auto&& i : this->m_wrappers )
  {
    i->resize(newsize);
  }
  *(this->getWrapper<localIndex>( keys::Size ).data())=newsize;
}



void ManagedGroup::RegisterDocumentationNodes()
{
  for( auto&& subNode : m_docNode->getChildNodes() )
  {
    if( subNode.second.getDataType() != "DocumentationNode" )
    {
      RegisterViewWrapper( subNode.second.getStringKey(),
                           rtTypes::typeID(subNode.second.getDataType() ) );
    }
  }
}

void ManagedGroup::BuildDataStructure( dataRepository::ManagedGroup * const rootGroup )
{
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->BuildDataStructure( rootGroup );
  }
}

void ManagedGroup::FillDocumentationNode( dataRepository::ManagedGroup * const  )
{

}


void ManagedGroup::SetDocumentationNodes( dataRepository::ManagedGroup * const group )
{
  FillDocumentationNode(group);
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->SetDocumentationNodes(group);
  }
}


}
} /* namespace ODS */
