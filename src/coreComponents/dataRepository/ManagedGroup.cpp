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
 * DataObjectManager.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: rrsettgast
 */

#include "ManagedGroup.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/TimingMacros.hpp"
#include <mpi.h>

#ifdef GEOSX_USE_ATK
#include "dataRepository/SidreWrapper.hpp"
#include "axom/sidre/core/sidre.hpp"
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
#ifdef GEOSX_USE_ATK
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
#ifdef GEOSX_USE_ATK
  m_sidreGroup(ManagedGroup::setSidreGroup(name,parent)),
#endif
  m_size(0),
  m_capacity(0),
  m_restart_flags(RestartFlags::WRITE_AND_READ),
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

//  RegisterViewWrapper( "size", &(this->m_size), false );
//  RegisterViewWrapper( "name", &(this->m_name), false );


  RegisterDocumentationNodes();
}

ManagedGroup::~ManagedGroup()
{}

ManagedGroup::ManagedGroup( ManagedGroup&& source ):
  m_parent( std::move(source.m_parent) ),
  m_wrappers( std::move(source.m_wrappers) ),
#ifdef GEOSX_USE_ATK
  m_sidreGroup( std::move(source.m_sidreGroup) ),
#endif
  m_size( source.m_size ),
  m_capacity( source.m_capacity ),
  m_restart_flags( source.m_restart_flags ),
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

ViewWrapperBase * ManagedGroup::RegisterViewWrapper( string const & name,
                                                     ViewWrapperBase * const wrapper )
{
  return m_wrappers.insert( name,
                            wrapper,
                            true );
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
  if( m_size > m_capacity )
  {
    m_capacity = m_size;
  }
}

void ManagedGroup::reserve( indexType const newsize )
{
  for( auto&& i : this->wrappers() )
  {
    if( i.second->sizedFromParent() == 1 )
    {
      i.second->reserve(newsize);
    }
  }
  m_capacity = newsize;
}


void ManagedGroup::RegisterDocumentationNodes()
{
  for( auto&& subNode : m_docNode->getChildNodes() )
  {
    if( ( subNode.second.getSchemaType().find("Node") == std::string::npos ) &&
        ( subNode.second.m_isRegistered == 0 ) )
    {
      std::string childType = subNode.second.getDataType();
      rtTypes::TypeIDs const typeID = rtTypes::typeID(childType);

      ViewWrapperBase * const view = RegisterViewWrapper( subNode.second.getStringKey(),
                                                          typeID );

      view->setSizedFromParent( subNode.second.m_managedByParent);
      view->setRestartFlags( subNode.second.getRestartFlags() );
      view->setPlotLevel( subNode.second.getVerbosity() );
      subNode.second.m_isRegistered = 1;

      string defVal = subNode.second.getDefault();

      if (defVal == "ExplicitDynamic")
      {
        double TEMP = 10;
        TEMP *= 2;
      }

      if( subNode.second.getIsInput() && defVal != "NONE" )
      {
        rtTypes::ApplyTypeLambda2 ( typeID,
                                    [&]( auto a, auto b ) -> void
        {
          ViewWrapper<decltype(a)>& dataView = ViewWrapper<decltype(a)>::cast(*view);
          std::vector<decltype(b)> values;
          stringutilities::StringToType( values, defVal );
          localIndex const size = LvArray::integer_conversion<localIndex>(values.size());
          dataView.resize( size );
          decltype(a) & data = dataView.reference();
          cxx_utilities::equateStlVector(data, values);
        });
      }



//      else if( defVal != "NONE" && defVal != "" )
//      {
//        rtTypes::ApplyTypeLambda2( typeID,
//                                   [&]( auto a, auto b ) -> void
//        {
//
//          ViewWrapper<decltype(a)>& dataView = ViewWrapper<decltype(a)>::cast(*view);
//          std::vector<decltype(b)> values;
//          stringutilities::StringToType( values, defVal );
//          localIndex const size = LvArray::integer_conversion<localIndex>(values.size());
//          if( size != 1 )
//          {
//            GEOS_ERROR("Expect size of default value to be a scalar");
//          }
//          decltype(a) & data = dataView.reference();
//          for( localIndex c=0 ; c<size ; ++c )
//          {
//            data[c] = values[0];
//          }
//          //data = values[0];
//        });
//      }
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
  for( auto&& subGroup : m_subGroups )
  {
    subGroup.second->SetDocumentationNodes();
  }
  RegisterDocumentationNodes();

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
        // GEOS_ERROR("group with name " + childName + " not found in " + this->getName());
      }
    }
  }
}


void ManagedGroup::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK("Child not recognized: " << childKey << ", " << childName);
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

void ManagedGroup::PrintDataHierarchy(integer indent)
{
  for( auto& view : this->wrappers() )
  {
    GEOS_LOG(string(indent, '\t')<<view.second->getName()<<", "<<view.second->get_typeid().name());
  }

  for( auto& group : this->m_subGroups )
  {
    GEOS_LOG(string(indent, '\t')<<group.first<<':');
    group.second->PrintDataHierarchy(indent + 1);
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
  static localIndex indent = 0;

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

void ManagedGroup::IntermediateInitializationRecursive( ManagedGroup * const rootGroup)
{
  IntermediateInitializationPreSubGroups(rootGroup);

  string_array initOrder;
  InitializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    this->GetGroup(groupName)->IntermediateInitializationRecursive( rootGroup );
  }

  IntermediateInitializationPostSubGroups( rootGroup );
}

void ManagedGroup::FinalInitializationRecursive( ManagedGroup * const rootGroup)
{
  FinalInitializationPreSubGroups(rootGroup);

  string_array initOrder;
  InitializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    this->GetGroup(groupName)->FinalInitializationRecursive( rootGroup );
  }

  FinalInitializationPostSubGroups( rootGroup );
}


localIndex ManagedGroup::PackSize( string_array const & wrapperNames,
                                   arrayView1d<localIndex const> const & packList,
                                   integer const recursive ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::PackSize(this->getName());

  packedSize += bufferOps::PackSize( string("Wrappers"));
  if( wrapperNames.size()==0 )
  {
    packedSize += bufferOps::PackSize( static_cast<int>(m_wrappers.size()) );
    for( auto const & wrapperPair : this->m_wrappers )
    {
      packedSize += bufferOps::PackSize( wrapperPair.first );
      if( packList.empty() )
      {
        packedSize += wrapperPair.second->PackSize();
      }
      else
      {
        packedSize += wrapperPair.second->PackSize(packList);
      }
    }
  }
  else
  {
    packedSize += bufferOps::PackSize( static_cast<int>(wrapperNames.size()) );
    for( auto const & wrapperName : wrapperNames )
    {
      ViewWrapperBase const * const wrapper = this->getWrapperBase(wrapperName);
      packedSize += bufferOps::PackSize( wrapperName );
      if( packList.empty() )
      {
        packedSize += wrapper->PackSize();
      }
      else
      {
        packedSize += wrapper->PackSize(packList);
      }
    }
  }
  if( recursive > 0 )
  {
    packedSize += bufferOps::PackSize( string("SubGroups"));
    packedSize += bufferOps::PackSize( m_subGroups.size() );
    for( auto const & keyGroupPair : this->m_subGroups )
    {
      packedSize += bufferOps::PackSize( keyGroupPair.first );
      packedSize += keyGroupPair.second->PackSize( wrapperNames, packList, recursive );
    }
  }

  return packedSize;
}


localIndex ManagedGroup::PackSize( string_array const & wrapperNames,
                                   integer const recursive ) const
{
  arrayView1d<localIndex> nullArray;
  return PackSize(wrapperNames,nullArray,recursive);
}


localIndex ManagedGroup::Pack( buffer_unit_type * & buffer,
                               string_array const & wrapperNames,
                               arrayView1d<localIndex const> const & packList,
                               integer const recursive ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack<true>( buffer, this->getName() );

  packedSize += bufferOps::Pack<true>( buffer, string("Wrappers") );
  if( wrapperNames.size()==0 )
  {
    packedSize += bufferOps::Pack<true>( buffer, m_wrappers.size() );
    for( auto const & wrapperPair : this->m_wrappers )
    {
      packedSize += bufferOps::Pack<true>( buffer, wrapperPair.first );
      if( packList.empty() )
      {
        packedSize += wrapperPair.second->Pack( buffer );
      }
      else
      {
        packedSize += wrapperPair.second->Pack( buffer, packList );
      }
    }
  }
  else
  {
    packedSize += bufferOps::Pack<true>( buffer, wrapperNames.size() );
    for( auto const & wrapperName : wrapperNames )
    {
      ViewWrapperBase const * const wrapper = this->getWrapperBase(wrapperName);
      packedSize += bufferOps::Pack<true>( buffer, wrapperName );
      if( packList.empty() )
      {
        packedSize += wrapper->Pack( buffer );
      }
      else
      {
        packedSize += wrapper->Pack( buffer, packList );
      }
    }
  }


  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack<true>( buffer, string("SubGroups") );
    packedSize += bufferOps::Pack<true>( buffer, m_subGroups.size() );
    for( auto const & keyGroupPair : this->m_subGroups )
    {
      packedSize += bufferOps::Pack<true>( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->Pack( buffer, wrapperNames, packList, recursive );
    }
  }

  return packedSize;
}

localIndex ManagedGroup::Pack( buffer_unit_type * & buffer,
                            string_array const & wrapperNames,
                            integer const recursive ) const
{
  arrayView1d<localIndex> nullArray;
  return Pack( buffer, wrapperNames, nullArray, recursive );
}

localIndex ManagedGroup::Unpack( buffer_unit_type const *& buffer,
                          arrayView1d<localIndex> & packList,
                          integer const recursive )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOS_ERROR_IF( groupName != this->getName(), "ManagedGroup::Unpack(): group names do not match");

  string wrappersLabel;
  unpackedSize += bufferOps::Unpack( buffer, wrappersLabel);
  GEOS_ERROR_IF( wrappersLabel != "Wrappers", "ManagedGroup::Unpack(): wrapper label incorrect");

  localIndex numWrappers;
  unpackedSize += bufferOps::Unpack( buffer, numWrappers);
  for( localIndex a=0 ; a<numWrappers ; ++a )
  {
    string wrapperName;
    unpackedSize += bufferOps::Unpack( buffer, wrapperName );
    ViewWrapperBase * const wrapper = this->getWrapperBase(wrapperName);
    wrapper->Unpack(buffer,packList);
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOS_ERROR_IF( subGroups != "SubGroups", "ManagedGroup::Unpack(): group names do not match");

    decltype( m_subGroups.size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOS_ERROR_IF( numSubGroups != m_subGroups.size(), "ManagedGroup::Unpack(): incorrect number of subGroups");

    for( auto const & index : this->m_subGroups )
    {
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup(subGroupName)->Unpack(buffer,packList,recursive);
    }
  }

  return unpackedSize;
}


void ManagedGroup::prepareToWrite() const
{
#ifdef GEOSX_USE_ATK
  if (getRestartFlags() == RestartFlags::NO_WRITE)
  {
    return;
  }

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
#ifdef GEOSX_USE_ATK
  if (getRestartFlags() == RestartFlags::NO_WRITE)
  {
    return;
  }

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
#ifdef GEOSX_USE_ATK
  if (getRestartFlags() != RestartFlags::WRITE_AND_READ)
  {
    return;
  }

  axom::sidre::View* temp = m_sidreGroup->getView("__size__");
  if ( temp != nullptr )
  {
    m_size = temp->getScalar();
    m_sidreGroup->destroyView("__size__");
  }

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
#ifdef GEOSX_USE_ATK
  if (getRestartFlags() != RestartFlags::WRITE_AND_READ)
  {
    return;
  }

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
