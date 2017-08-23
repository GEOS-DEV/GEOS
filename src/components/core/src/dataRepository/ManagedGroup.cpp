/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "ManagedGroup.hpp"

#include "dataRepository/SidreWrapper.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "spio/IOManager.hpp"
#include <mpi.h>



namespace geosx
{
namespace dataRepository
{

axom::sidre::Group * ManagedGroup::setSidreGroup( string const& name,
                                                            ManagedGroup * const parent )
{
  axom::sidre::Group * sidreParent = nullptr;
  axom::sidre::Group * sidreGroup  = nullptr;

  if( parent==nullptr )
  {
    sidreParent = SidreWrapper::dataStore().getRoot();
  }
  else
  {
    sidreParent = parent->m_sidreGroup;
  }

  if( sidreParent->hasGroup(name) )
  {
    sidreGroup = sidreParent->getGroup(name);
  }
  else
  {
    sidreGroup = sidreParent->createGroup(name);
  }
  return sidreGroup;
}

ManagedGroup::ManagedGroup( std::string const & name,
                            ManagedGroup * const parent ) :
  m_docNode(nullptr),
  m_keyLookup(),
  m_wrappers(),
  m_parent(parent),
  m_subGroups(),
  m_sidreGroup(ManagedGroup::setSidreGroup(name,parent)),
  m_size(0)
{

  // Setup DocumentationNode
  if( parent != nullptr )
  {
    if( parent->m_docNode != nullptr  )
    {
      if( this->m_docNode == nullptr )
      {
        m_docNode = parent->m_docNode->AllocateChildNode( name,
                                                          name,
                                                          0,
                                                          "",
                                                          "Node",
                                                          "",
                                                          "",
                                                          "",
                                                          parent->getName(),
                                                          0,
                                                          0,
                                                          0 ) ;
      }
    }
    else
    {
      m_docNode = new cxx_utilities::DocumentationNode( name,
                                                        name,
                                                        -1,
                                                        "Node",
                                                        "Node",
                                                        "The Root DocumentationNode for " + name,
                                                        "",
                                                        "",
                                                        "",
                                                        0,
                                                        0,
                                                        0,
                                                        0,
                                                        nullptr );
    }
  }
  else
  {
    m_docNode = new cxx_utilities::DocumentationNode( name,
                                                      name,
                                                      -1,
                                                      "Node",
                                                      "Node",
                                                      "The Root DocumentationNode for " + name,
                                                      "",
                                                      "",
                                                      "",
                                                      0,
                                                      0,
                                                      0,
                                                      0,
                                                      nullptr );
  }


  m_docNode->AllocateChildNode( "size",
                                "size",
                                -1,
                                "int32",
                                "int32",
                                "size of group",
                                "Number of entries in this group.",
                                "0",
                                "",
                                0,
                                0,
                                2 );

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
                                0,
                                0 );

  m_docNode->AllocateChildNode( "path",
                                "path",
                                -1,
                                "string",
                                "string",
                                "path in hierarchy",
                                "path in hierarchy",
                                "",
                                "",
                                0,
                                0,
                                2 );

  RegisterDocumentationNodes();
}


//ManagedGroup::ManagedGroup( std::string const & name,
//                            ManagedGroup * const parent,
//                            cxx_utilities::DocumentationNode * docNode ) :
//  m_docNode(docNode),
//  m_keyLookup(),
//  m_wrappers(),
//  m_parent(parent),
//  m_subGroups(),
//  m_sidreGroup(nullptr)
//{
//
//  // SIDRE interaction
//  asctoolkit::sidre::Group * sidreParent = nullptr;
//  if( m_parent==nullptr )
//  {
//    sidreParent = SidreWrapper::dataStore().getRoot();
//  }
//  else
//  {
//    sidreParent = parent->m_sidreGroup;
//  }
//
//  if( sidreParent->hasGroup(name) )
//  {
//    m_sidreGroup = sidreParent->getGroup(name);
//  }
//  else
//  {
//    m_sidreGroup = sidreParent->createGroup(name);
//  }
//
//
//
//  // Setup DocumentationNode
//  if( parent != nullptr )
//  {
//    if( parent->m_docNode != nullptr && this->m_docNode != nullptr )
//    {
//      m_docNode = parent->m_docNode->AllocateChildNode( name,
//                                                        name,
//                                                        0,
//                                                        "ManagedGroup",
//                                                        "Node",
//                                                        "ManagedGroup",
//                                                        "ManagedGroup",
//                                                        "",
//                                                        parent->getName(),
//                                                        0,
//                                                        0 ) ;
//
//    }
//    else
//    {
//
//    }
//  }
//
//
//  m_docNode->AllocateChildNode( "size",
//                                "size",
//                                -1,
//                                "int32",
//                                "int32",
//                                "size of group",
//                                "Number of entries in this group.",
//                                "0",
//                                "",
//                                0,
//                                0 );
//
//  m_docNode->AllocateChildNode( "name",
//                                "name",
//                                -1,
//                                "string",
//                                "string",
//                                "name of group",
//                                "name of group.",
//                                name,
//                                "",
//                                0,
//                                0 );
//
//  *(RegisterViewWrapper<int32>( "size" ).data()) = 0;
//  RegisterViewWrapper<std::string>( "name" ).reference() = name;
//  RegisterViewWrapper<std::string>( "path" );
//
//}

ManagedGroup::~ManagedGroup()
{
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
  m_parent( std::move(source.m_parent) ),
  m_size( source.m_size )
{}



ManagedGroup::CatalogInterface::CatalogType& ManagedGroup::GetCatalog()
{
  static ManagedGroup::CatalogInterface::CatalogType catalog;
  return catalog;
}

ViewWrapperBase& ManagedGroup::RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type )
{
  return *( rtTypes::ApplyTypeLambda1( type,
                                       [this, &name]( auto a ) -> ViewWrapperBase*
      {
        return &( this->RegisterViewWrapper<decltype(a)>(name) );
      } ) );
}

void ManagedGroup::resize( int32 const newsize )
{
  for( auto&& i : this->m_wrappers )
  {
    if( i->sizedFromParent() == 1 )
    {
      i->resize(newsize);
    }
  }
  m_size = newsize;
}



void ManagedGroup::RegisterDocumentationNodes()
{
//  std::cout<<std::string(m_docNode->m_level*2, ' ')<<"Registering Documentation Nodes for Group "<<this->getName()<<std::endl;
  for( auto&& subNode : m_docNode->getChildNodes() )
  {
//    std::cout<<subNode.first<<", "<<subNode.second.getName()<<std::endl;
    if( ( subNode.second.getSchemaType() != "DocumentationNode" ) &&
        ( subNode.second.getSchemaType() != "Node" ) &&
        ( subNode.second.m_isRegistered == 0 ) )
    {
//      std::cout<<std::string(subNode.second.m_level*2, ' ')<<"Register "<<subNode.second.getStringKey()<<" of type "<<subNode.second.getDataType()<<std::endl;
      ViewWrapperBase & view = RegisterViewWrapper( subNode.second.getStringKey(),
                                                    rtTypes::typeID(subNode.second.getDataType() ) );
      view.setSizedFromParent( subNode.second.m_managedByParent);
      subNode.second.m_isRegistered = 1;
    }
  }

  for( auto& subGroupIter : m_subGroups )
  {
    subGroupIter.second->RegisterDocumentationNodes();
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
  RegisterDocumentationNodes();
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->SetDocumentationNodes(group);
  }
}


void ManagedGroup::ReadXML( xmlWrapper::xmlNode const & targetNode )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  
  ReadXMLsub( targetNode );

  for( auto const & subDocEntry : docNode->m_child )
  {
    cxx_utilities::DocumentationNode subDocNode = subDocEntry.second;
    
    if (subDocNode.getIsInput() == 1)
    {
      xmlWrapper::ReadAttributeAsType( *this, subDocNode, targetNode );
    }
  }

  ReadXML_PostProcess();
}

void ManagedGroup::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  this->forSubGroups( [this,&targetNode]( ManagedGroup & subGroup ) -> void
  {
    subGroup.ReadXML( targetNode );
  });
}


void ManagedGroup::PrintDataHierarchy()
{
  for( auto& view : this->m_wrappers )
  {
    std::cout<<view->getName()<<", "<<view->get_typeid().name()<<std::endl;
  }

  for( auto& group : this->m_subGroups )
  {
    std::cout<<group.first<<std::endl;
    group.second->PrintDataHierarchy();
  }
}

void ManagedGroup::InitializationOrder( string_array & order )
{
  for( auto & subGroup : this->m_subGroups )
  {
    order.push_back(subGroup.first);
  }
}

void ManagedGroup::Initialize( ManagedGroup * const group )
{
  static int indent = 0;
  std::cout<<string(indent*2, ' ')<<"Calling ManagedGroup::Initialize() on "<<this->getName()<<" of type "<<cxx_utilities::demangle(this->get_typeid().name())<<std::endl;

  InitializePreSubGroups(group);

  string_array initOrder;
  InitializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    ++indent;
    this->GetGroup(groupName).Initialize(group);
    --indent;
  }

//  forSubGroups( [&]( ManagedGroup & subGroup ) -> void
//  {
//    ++indent;
//    subGroup.Initialize(group);
//    --indent;
//  });
  InitializePostSubGroups(group);
}

/* Add pointers to ViewWrapper data to the sidre tree. */
void ManagedGroup::registerSubViews() 
{
  for (std::unique_ptr<ViewWrapperBase> & wrapper : m_wrappers) 
  {
    wrapper->registerDataPtr();
  }

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.registerSubViews();
  });
}

/* Remove pointers to ViewWrapper data from the sidre tree. */
void ManagedGroup::unregisterSubViews()
{
  for (std::unique_ptr<ViewWrapperBase> & wrapper : m_wrappers) 
  {
    wrapper->unregisterDataPtr();
  }

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.unregisterSubViews();
  });
}

/* Save m_size to sidre views. */
void ManagedGroup::createSizeViews()
{
  m_sidreGroup->createView("__size__")->setScalar(m_size);

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.createSizeViews();
  });
}

/* Load m_size data from sidre views. */
void ManagedGroup::loadSizeViews()
{
  m_size = m_sidreGroup->getView("__size__")->getScalar();
  m_sidreGroup->destroyView("__size__");

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.loadSizeViews();
  });
}

/* Resize views to hold data from sidre. */
void ManagedGroup::resizeSubViews() 
{
  for (std::unique_ptr<ViewWrapperBase> & wrapper : m_wrappers) 
  {
    wrapper->resizeFromSidre();
  }

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.resizeSubViews();
  });
}


void ManagedGroup::storeSizedFromParent()
{
  for (std::unique_ptr<ViewWrapperBase> & wrapper : m_wrappers) 
  {
    wrapper->storeSizedFromParent();
  }

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.storeSizedFromParent();
  });
}

void ManagedGroup::loadSizedFromParent()
{
  for (std::unique_ptr<ViewWrapperBase> & wrapper : m_wrappers) 
  {
    wrapper->loadSizedFromParent();
  }
  

  forSubGroups([](ManagedGroup & subGroup) -> void 
  {
    subGroup.loadSizedFromParent();
  });
}

/* Write out a restart file. */
void ManagedGroup::writeRestart(int num_files, const string & path, const string & protocol, MPI_Comm comm) 
{
  SidreWrapper::dataStore().createAttributeScalar("__sizedFromParent__", -1);
  storeSizedFromParent();
  registerSubViews();
  createSizeViews();
  axom::spio::IOManager ioManager(comm);
  ioManager.write(m_sidreGroup, num_files, path, protocol);
}

/* Read in a restart file and reconstruct the sidre tree. */
void ManagedGroup::reconstructSidreTree(const string & root_path, const string & protocol, MPI_Comm comm)
{
  axom::spio::IOManager ioManager(comm);
  ioManager.read(m_sidreGroup, root_path, protocol);
}

/* Load sidre external data. */
void ManagedGroup::loadSidreExternalData(const string & root_path, MPI_Comm comm)
{
  loadSizedFromParent();
  SidreWrapper::dataStore().destroyAttribute("__sizedFromParent__");
  loadSizeViews();
  resizeSubViews();
  registerSubViews();
  axom::spio::IOManager ioManager(comm);
  ioManager.loadExternalData(m_sidreGroup, root_path);
  unregisterSubViews();
}


}
} /* namespace ODS */
