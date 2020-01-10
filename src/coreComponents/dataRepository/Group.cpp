/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "Group.hpp"
#include "ConduitRestart.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "common/TimingMacros.hpp"


namespace geosx
{
namespace dataRepository
{

conduit::Node & conduitNodeFromParent( string const & name, Group * const parent )
{
  if( parent == nullptr )
  {
    return rootConduitNode[ name ];
  }
  else
  {
    return parent->getConduitNode()[ name ];
  }
}

Group::Group( std::string const & name,
              Group * const parent ):
  m_parent( parent ),
  m_wrappers(),
  m_subGroups(),
  m_size( 0 ),
  m_capacity( 0 ),
  m_name( name ),
  m_logLevel( 0 ),
  m_restart_flags( RestartFlags::WRITE_AND_READ ),
  m_input_flags( InputFlags::INVALID ),
  m_conduitNode( conduitNodeFromParent( name, parent ) )
{}

Group::~Group()
{}

Group::CatalogInterface::CatalogType & Group::GetCatalog()
{
  static Group::CatalogInterface::CatalogType catalog;
  return catalog;
}

WrapperBase * Group::registerWrapper( string const & name,
                                      WrapperBase * const wrapper )
{
  return m_wrappers.insert( name,
                            wrapper,
                            true );
}

void Group::deregisterWrapper( string const & name )
{
  GEOSX_ERROR_IF( !hasWrapper( name ), "Wrapper " << name << " doesn't exist." );
  m_wrappers.erase( name );
  m_conduitNode.remove( name );
}


void Group::resize( indexType const newsize )
{
  for( auto && i : this->wrappers() )
  {
    if( i.second->sizedFromParent() == 1 )
    {
      i.second->resize( newsize );
    }
  }
  m_size = newsize;
  if( m_size > m_capacity )
  {
    m_capacity = m_size;
  }
}

void Group::reserve( indexType const newsize )
{
  for( auto && i : this->wrappers() )
  {
    if( i.second->sizedFromParent() == 1 )
    {
      i.second->reserve( newsize );
    }
  }
  m_capacity = newsize;
}

void Group::ProcessInputFileRecursive( xmlWrapper::xmlNode & targetNode )
{
  xmlWrapper::addIncludedXML( targetNode );

  // loop over the child nodes of the targetNode
  for( xmlWrapper::xmlNode childNode=targetNode.first_child() ; childNode ; childNode=childNode.next_sibling())
  {
    // Get the child tag and name
    std::string childName = childNode.attribute( "name" ).value();
    if( childName.empty())
    {
      childName = childNode.name();
    }

    // Create children
    Group * newChild = CreateChild( childNode.name(), childName );
    if( newChild == nullptr )
    {
      newChild = GetGroup( childName );
    }
    if( newChild != nullptr )
    {
      newChild->ProcessInputFileRecursive( childNode );
    }
  }

  ProcessInputFile( targetNode );
//  ProcessInputFile_PostProcess();
}

void Group::ProcessInputFile( xmlWrapper::xmlNode const & targetNode )
{

  std::set< string > processedXmlNodes;
  for( auto wrapperPair : m_wrappers )
  {
    WrapperBase * const wrapper = wrapperPair.second;
    InputFlags const inputFlag = wrapper->getInputFlag();
    if( inputFlag >= InputFlags::OPTIONAL )
    {
      string const & wrapperName = wrapperPair.first;
      rtTypes::TypeIDs const wrapperTypeID = rtTypes::typeID( wrapper->get_typeid());

      rtTypes::ApplyIntrinsicTypeLambda2( wrapperTypeID,
                                          [&]( auto a, auto GEOSX_UNUSED_ARG( b ) )
          {
            using COMPOSITE_TYPE = decltype( a );

            Wrapper< COMPOSITE_TYPE > & typedWrapper = Wrapper< COMPOSITE_TYPE >::cast( *wrapper );
            COMPOSITE_TYPE & objectReference = typedWrapper.reference();
            processedXmlNodes.insert( wrapperName );

            if( inputFlag == InputFlags::REQUIRED || !(typedWrapper.getDefaultValueStruct().has_default_value) )
            {
              bool const readSuccess = xmlWrapper::ReadAttributeAsType( objectReference,
                                                                        wrapperName,
                                                                        targetNode,
                                                                        inputFlag == InputFlags::REQUIRED );
              GEOSX_ERROR_IF( !readSuccess,
                              "Input variable " + wrapperName + " is required in " + targetNode.path()
                              + ". Available options are: \n"+ dumpInputOptions()
                              + "\nFor more details, please refer to documentation at: \n"
                              + "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html \n" );


            }
            else
            {
              xmlWrapper::ReadAttributeAsType( objectReference, wrapperName, targetNode, typedWrapper.getDefaultValueStruct() );
            }
          } );
    }
  }

  for( xmlWrapper::xmlAttribute attribute=targetNode.first_attribute() ; attribute ; attribute = attribute.next_attribute() )
  {
    string const childName = attribute.name();
    if( childName != "name" && childName != "xmlns:xsi" && childName != "xsi:noNamespaceSchemaLocation" )
    {
      GEOSX_ERROR_IF( processedXmlNodes.count( childName )==0,
                      "XML Node ("<<targetNode.name()<<") with attribute name=("<<
                      targetNode.attribute( "name" ).value()<<") contains child node named ("<<
                      childName<<") that is not read. Valid options are: \n" << dumpInputOptions()
                      + "\nFor more details, please refer to documentation at: \n"
                      + "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html \n" );
    }
  }

}

void Group::PostProcessInputRecursive()
{
  for( auto const & subGroupIter : m_subGroups )
  {
    subGroupIter.second->PostProcessInputRecursive();
  }
  PostProcessInput();
}



void Group::RegisterDataOnMeshRecursive( Group * const meshBodies )
{
  RegisterDataOnMesh( meshBodies );
  for( auto && subGroup : m_subGroups )
  {
    subGroup.second->RegisterDataOnMeshRecursive( meshBodies );
  }
}


Group * Group::CreateChild( string const & childKey, string const & childName )
{
  GEOSX_ERROR_IF( !(CatalogInterface::hasKeyName( childKey )),
                  "KeyName ("<<childKey<<") not found in Group::Catalog" );
  GEOSX_LOG_RANK_0( "Adding Object " << childKey<<" named "<< childName<<" from Group::Catalog." );
  return RegisterGroup( childName,
                        CatalogInterface::Factory( childKey, childName, this ) );
}


void Group::PrintDataHierarchy( integer indent )
{
  for( auto & view : this->wrappers() )
  {
    GEOSX_LOG( string( indent, '\t' )<<view.second->getName()<<", "<<view.second->get_typeid().name());
  }

  for( auto & group : this->m_subGroups )
  {
    GEOSX_LOG( string( indent, '\t' )<<group.first<<':' );
    group.second->PrintDataHierarchy( indent + 1 );
  }
}

string Group::dumpInputOptions()
{
  string rval;
  char temp[1000] = {0};
  sprintf( temp, "  |         name         |  opt/req  | Description \n" );
  rval.append( temp );
  sprintf( temp, "  |----------------------|-----------|-----------------------------------------\n" );
  rval.append( temp );

  for( auto const & wrapper : m_wrappers )
  {
    WrapperBase const * const wb = wrapper.second;
    if( wb->getInputFlag() == InputFlags::OPTIONAL ||
        wb->getInputFlag() == InputFlags::REQUIRED )
    {
      sprintf( temp, "  | %20s | %9s | %s \n",
               wb->getName().c_str(),
               InputFlagToString( wb->getInputFlag()).c_str(),
               wb->getDescription().c_str() );
      rval.append( temp );
    }
  }
  return rval;
}


void Group::InitializationOrder( string_array & order )
{
  for( auto & subGroupIter : this->m_subGroups )
  {
    order.push_back( subGroupIter.first );
  }
}

void Group::Initialize( Group * const group )
{
  static localIndex indent = 0;

  InitializePreSubGroups( group );

  string_array initOrder;
  InitializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    ++indent;
    this->GetGroup( groupName )->Initialize( group );
    --indent;
  }

  InitializePostSubGroups( group );
}


void Group::InitializePostInitialConditions( Group * const rootGroup )
{
  InitializePostInitialConditions_PreSubGroups( rootGroup );

  string_array initOrder;
  InitializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    this->GetGroup( groupName )->InitializePostInitialConditions( rootGroup );
  }

  InitializePostInitialConditions_PostSubGroups( rootGroup );
}

localIndex Group::PackSize( string_array const & wrapperNames,
                            arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool on_device ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::PackSize( this->getName());

  packedSize += bufferOps::PackSize( string( "Wrappers" ));
  if( wrapperNames.size()==0 )
  {
    packedSize += bufferOps::PackSize( static_cast< int >(m_wrappers.size()) );
    for( auto const & wrapperPair : this->m_wrappers )
    {
      packedSize += bufferOps::PackSize( wrapperPair.first );
      if( packList.empty() )
      {
        packedSize += wrapperPair.second->PackSize( on_device );
      }
      else
      {
        packedSize += wrapperPair.second->PackByIndexSize( packList, on_device );
      }
    }
  }
  else
  {
    packedSize += bufferOps::PackSize( static_cast< int >(wrapperNames.size()) );
    for( auto const & wrapperName : wrapperNames )
    {
      WrapperBase const * const wrapper = this->getWrapperBase( wrapperName );
      packedSize += bufferOps::PackSize( wrapperName );
      if( packList.empty() )
      {
        packedSize += wrapper->PackSize( on_device );
      }
      else
      {
        packedSize += wrapper->PackByIndexSize( packList, on_device );
      }
    }
  }
  if( recursive > 0 )
  {
    packedSize += bufferOps::PackSize( string( "SubGroups" ));
    packedSize += bufferOps::PackSize( m_subGroups.size() );
    for( auto const & keyGroupPair : this->m_subGroups )
    {
      packedSize += bufferOps::PackSize( keyGroupPair.first );
      packedSize += keyGroupPair.second->PackSize( wrapperNames, packList, recursive, on_device );
    }
  }

  return packedSize;
}


localIndex Group::PackSize( string_array const & wrapperNames,
                            integer const recursive,
                            bool on_device ) const
{
  arrayView1d< localIndex > nullArray;
  return PackSize( wrapperNames, nullArray, recursive, on_device );
}


localIndex Group::Pack( buffer_unit_type * & buffer,
                        string_array const & wrapperNames,
                        arrayView1d< localIndex const > const & packList,
                        integer const recursive,
                        bool on_device ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< true >( buffer, this->getName() );

  packedSize += bufferOps::Pack< true >( buffer, string( "Wrappers" ) );
  if( wrapperNames.size()==0 )
  {
    packedSize += bufferOps::Pack< true >( buffer, m_wrappers.size() );
    for( auto const & wrapperPair : this->m_wrappers )
    {
      packedSize += bufferOps::Pack< true >( buffer, wrapperPair.first );
      if( packList.empty() )
      {
        // invoke wrapper pack kernel
        packedSize += wrapperPair.second->Pack( buffer, on_device );
      }
      else
      {
        packedSize += wrapperPair.second->PackByIndex( buffer, packList, on_device );
      }
    }
  }
  else
  {
    packedSize += bufferOps::Pack< true >( buffer, wrapperNames.size() );
    for( auto const & wrapperName : wrapperNames )
    {
      WrapperBase const * const wrapper = this->getWrapperBase( wrapperName );
      packedSize += bufferOps::Pack< true >( buffer, wrapperName );
      if( packList.empty() )
      {
        packedSize += wrapper->Pack( buffer, on_device );
      }
      else
      {
        packedSize += wrapper->PackByIndex( buffer, packList, on_device );
      }
    }
  }


  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< true >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< true >( buffer, m_subGroups.size() );
    for( auto const & keyGroupPair : this->m_subGroups )
    {
      packedSize += bufferOps::Pack< true >( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->Pack( buffer, wrapperNames, packList, recursive, on_device );
    }
  }

  return packedSize;
}

localIndex Group::Pack( buffer_unit_type * & buffer,
                        string_array const & wrapperNames,
                        integer const recursive,
                        bool on_device ) const
{
  arrayView1d< localIndex > nullArray;
  return Pack( buffer, wrapperNames, nullArray, recursive, on_device );
}

localIndex Group::Unpack( buffer_unit_type const * & buffer,
                          arrayView1d< localIndex > & packList,
                          integer const recursive,
                          bool on_device )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF( groupName != this->getName(), "Group::Unpack(): group names do not match" );

  string wrappersLabel;
  unpackedSize += bufferOps::Unpack( buffer, wrappersLabel );
  GEOSX_ERROR_IF( wrappersLabel != "Wrappers", "Group::Unpack(): wrapper label incorrect" );

  localIndex numWrappers;
  unpackedSize += bufferOps::Unpack( buffer, numWrappers );
  for( localIndex a=0 ; a<numWrappers ; ++a )
  {
    string wrapperName;
    unpackedSize += bufferOps::Unpack( buffer, wrapperName );
    WrapperBase * const wrapper = this->getWrapperBase( wrapperName );
    wrapper->UnpackByIndex( buffer, packList, on_device );
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOSX_ERROR_IF( subGroups != "SubGroups", "Group::Unpack(): group names do not match" );

    decltype( m_subGroups.size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF( numSubGroups != m_subGroups.size(), "Group::Unpack(): incorrect number of subGroups" );

    for( auto const & index : this->m_subGroups )
    {
      GEOSX_UNUSED_VAR( index );
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += this->GetGroup( subGroupName )->Unpack( buffer, packList, recursive, on_device );
    }
  }

  return unpackedSize;
}


void Group::prepareToWrite()
{
  if( getRestartFlags() == RestartFlags::NO_WRITE )
  {
    return;
  }

  for( auto & pair : m_wrappers )
  {
    pair.second->registerToWrite();
  }

  m_conduitNode[ "__size__" ].set( m_size );

  forSubGroups( []( Group * subGroup )
      {
        subGroup->prepareToWrite();
      } );
}


void Group::finishWriting()
{
  if( getRestartFlags() == RestartFlags::NO_WRITE )
  {
    return;
  }

  for( auto & pair : m_wrappers )
  {
    pair.second->finishWriting();
  }

  forSubGroups( []( Group * subGroup )
      {
        subGroup->finishWriting();
      } );
}


void Group::loadFromConduit()
{
  if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
  {
    return;
  }

  m_size = m_conduitNode.fetch_child( "__size__" ).value();

  for( auto & pair : m_wrappers )
  {
    pair.second->loadFromConduit();
  }

  forSubGroups( []( Group * subGroup )
      {
        subGroup->loadFromConduit();
      } );
}

void Group::postRestartInitializationRecursive( Group * const domain )
{
  forSubGroups( [&]( Group * const subGroup )
      {
        subGroup->postRestartInitializationRecursive( domain );
      } );

  this->postRestartInitialization( domain );
}

void Group::SetSchemaDeviations( xmlWrapper::xmlNode GEOSX_UNUSED_ARG( schemaRoot ),
                                 xmlWrapper::xmlNode GEOSX_UNUSED_ARG( schemaParent ),
                                 integer GEOSX_UNUSED_ARG( documentationType ) )
{}

void Group::RegisterDataOnMesh( Group * const GEOSX_UNUSED_ARG( MeshBody ) )
{}


void Group::enableLogLevelInput()
{
  string const logLevelString = "logLevel";

  registerWrapper( logLevelString, &m_logLevel, false )->
    setApplyDefaultValue( 0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Log level" );
}

} /* end namespace dataRepository */
} /* end namespace geosx  */
