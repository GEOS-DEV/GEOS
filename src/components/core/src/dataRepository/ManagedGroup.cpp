/*
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "ManagedGroup.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include <mpi.h>

#ifdef USE_ATK
#include "dataRepository/SidreWrapper.hpp"
#include "spio/IOManager.hpp"
#endif

namespace geosx
{
namespace dataRepository
{

axom::sidre::Group * ManagedGroup::setSidreGroup( string const& name,
                                                  ManagedGroup * const parent )
{

  axom::sidre::Group * sidreParent = nullptr;
  axom::sidre::Group * sidreGroup  = nullptr;
#ifdef USE_ATK
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
#else
  return nullptr;

#endif
}

ManagedGroup::ManagedGroup( std::string const & name,
                            ManagedGroup * const parent ):
  m_docNode(nullptr),
  m_parent(parent),
  m_wrappers(),
  m_subGroups(),
#ifdef USE_ATK
  m_sidreGroup(ManagedGroup::setSidreGroup(name,parent)),
#endif
  m_size(0),
  m_name(name)
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
                                                          0 );
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

  RegisterDocumentationNodes();
}

ManagedGroup::~ManagedGroup()
{}

ManagedGroup::ManagedGroup( ManagedGroup&& source ):
  m_parent( std::move(source.m_parent) ),
  m_wrappers( std::move(source.m_wrappers) ),
#ifdef USE_ATK
  m_sidreGroup( std::move(source.m_sidreGroup) ),
#endif
  m_size( source.m_size ),
  m_name( source.m_name )
{}



ManagedGroup::CatalogInterface::CatalogType& ManagedGroup::GetCatalog()
{
  static ManagedGroup::CatalogInterface::CatalogType catalog;
  return catalog;
}

ViewWrapperBase * ManagedGroup::RegisterViewWrapper( std::string const & name, rtTypes::TypeIDs const & type )
{
  return rtTypes::ApplyTypeLambda1( type,
                                    [this, &name]( auto a ) -> ViewWrapperBase*
      {
        return this->RegisterViewWrapper<decltype(a)>(name);
      } );
}

void ManagedGroup::resize( indexType const newsize )
{
  for( auto&& i : this->wrappers() )
  {
    if( i.second->sizedFromParent() == 1 )
    {
      i.second->resize(newsize);
    }
  }
  m_size = newsize;
}



void ManagedGroup::RegisterDocumentationNodes()
{
  for( auto&& subNode : m_docNode->getChildNodes() )
  {
    if( ( subNode.second.getSchemaType().find("Node") == std::string::npos ) &&
        ( subNode.second.m_isRegistered == 0 ) )
    {
      ViewWrapperBase * const view = RegisterViewWrapper( subNode.second.getStringKey(),
                                                          rtTypes::typeID(subNode.second.getDataType() ) );
      view->setSizedFromParent( subNode.second.m_managedByParent);
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

// These fill the documentation and initialize fields on this:
void ManagedGroup::FillDocumentationNode()
{}

void ManagedGroup::SetDocumentationNodes()
{
  FillDocumentationNode();
  RegisterDocumentationNodes();
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->SetDocumentationNodes();
  }
}

// These fill the documentation and initialize fields on other objects:
void ManagedGroup::SetOtherDocumentationNodes(dataRepository::ManagedGroup * const rootGroup)
{
  FillOtherDocumentationNodes(rootGroup);
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->SetOtherDocumentationNodes(rootGroup);
  }
}

void ManagedGroup::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const )
{
}



void ManagedGroup::AddChildren( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    // Get the child tag and name
    std::string childName = childNode.attribute("name").value();
    if (childName.empty())
    {
      childName = childNode.name();
    }

    // Create children
    CreateChild(childNode.name(), childName);

    // Add grandchildren
    ManagedGroup * newChild = this->GetGroup<ManagedGroup>(childName);
    if (newChild != nullptr)
    {
      newChild->AddChildren(childNode);
    }
    else
    {
      if( !this->hasView(childName) )
      {
        //GEOS_ERROR("group with name " + childName + " not found in " + this->getName());
      }
    }
  }
}


void ManagedGroup::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Child not recognized: " << childKey << ", " << childName << std::endl;
}


void ManagedGroup::ReadXML( xmlWrapper::xmlNode const & targetNode )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  for( auto const & subDocEntry : docNode->m_child )
  {
    cxx_utilities::DocumentationNode subDocNode = subDocEntry.second;

    if (subDocNode.getIsInput() == 1)
    {
      xmlWrapper::ReadAttributeAsType( *this, subDocNode, targetNode );
    }
  }

  ReadXMLsub(targetNode);
  ReadXML_PostProcess();
}



void ManagedGroup::ReadXMLsub( xmlWrapper::xmlNode const & targetNode )
{
  for (xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    // Get the child tag and name
    std::string childName = childNode.attribute("name").value();
    if (childName.empty())
    {
      childName = childNode.name();
    }

    // Read the xml on children
    ManagedGroup * child = this->GetGroup<ManagedGroup>(childName);
    if (child != nullptr)
    {
      child->ReadXML(childNode);
    }
  }
}

void ManagedGroup::PrintDataHierarchy()
{
  for( auto& view : this->wrappers() )
  {
    std::cout<<view.second->getName()<<", "<<view.second->get_typeid().name()<<std::endl;
  }

  for( auto& group : this->m_subGroups )
  {
    group.second->PrintDataHierarchy();
  }
}

void ManagedGroup::InitializationOrder( string_array & order )
{
  for( auto & subGroupIter : this->m_subGroups )
  {
    order.push_back(subGroupIter.first);
  }
}

void ManagedGroup::Initialize( ManagedGroup * const group )
{
  static int indent = 0;
//  std::cout<<string(indent*2, ' ')<<"Calling ManagedGroup::Initialize() on
// "<<this->getName()<<" of type
// "<<cxx_utilities::demangle(this->get_typeid().name())<<std::endl;

  InitializePreSubGroups(group);

  string_array initOrder;
  InitializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    ++indent;
    this->GetGroup(groupName)->Initialize(group);
    --indent;
  }

//  forSubGroups( [&]( ManagedGroup * subGroup ) -> void
//  {
//    ++indent;
//    subGroup->Initialize(group);
//    --indent;
//  });
  InitializePostSubGroups(group);
}


void ManagedGroup::prepareToWrite() const
{
#ifdef USE_ATK
  if (!SidreWrapper::dataStore().hasAttribute("__sizedFromParent__"))
  {
    SidreWrapper::dataStore().createAttributeScalar("__sizedFromParent__", -1);
  }

  for (auto & pair : m_wrappers)
  {
    pair.second->registerToWrite();
  }

  if ( m_sidreGroup->hasView("__size__") ) {
    m_sidreGroup->getView("__size__")->setScalar(m_size);
  } else {
    m_sidreGroup->createView("__size__")->setScalar(m_size);  
  }

  forSubGroups([](const ManagedGroup * subGroup) -> void 
  {
    subGroup->prepareToWrite();
  });
#endif
}


void ManagedGroup::finishWriting() const
{
#ifdef USE_ATK
  axom::sidre::View* temp = m_sidreGroup->getView("__size__");
  m_sidreGroup->destroyView("__size__");

  for (auto & pair : m_wrappers)
  {
    pair.second->finishWriting();
  }

  forSubGroups([](const ManagedGroup * subGroup) -> void 
  {
    subGroup->finishWriting();
  });
#endif
}


void ManagedGroup::prepareToRead()
{
#ifdef USE_ATK
  axom::sidre::View* temp = m_sidreGroup->getView("__size__");
  m_size = temp->getScalar();
  m_sidreGroup->destroyView("__size__");

  for (auto & pair : m_wrappers)
  {
    pair.second->registerToRead();
  }

  forSubGroups([](ManagedGroup * subGroup) -> void 
  {
    subGroup->prepareToRead();
  });  
#endif
}


void ManagedGroup::finishReading()
{
#ifdef USE_ATK
  for (auto & pair : m_wrappers)
  {
    pair.second->finishReading();
  }

  forSubGroups([](ManagedGroup * subGroup) -> void 
  {
    subGroup->finishReading();
  });
#endif
}


} /* end namespace dataRepository */
} /* end namespace geosx  */
