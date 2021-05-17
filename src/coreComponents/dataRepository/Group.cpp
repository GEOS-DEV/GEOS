/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

Group::Group( string const & name,
              Group * const parent ):
  Group( name, parent->getConduitNode() )
{
  GEOSX_ERROR_IF( parent == nullptr, "Should not be null." );
  m_parent = parent;
}

Group::Group( string const & name,
              conduit::Node & rootNode ):
  m_parent( nullptr ),
  m_sizedFromParent( 0 ),
  m_wrappers(),
  m_subGroups(),
  m_size( 0 ),
  m_capacity( 0 ),
  m_name( name ),
  m_logLevel( 0 ),
  m_restart_flags( RestartFlags::WRITE_AND_READ ),
  m_input_flags( InputFlags::INVALID ),
  m_conduitNode( rootNode[ name ] )
{}

Group::~Group()
{
// TODO enable this and fix bugs this exposes.
//  m_conduitNode.parent()->remove( m_name );
}

Group::CatalogInterface::CatalogType & Group::getCatalog()
{
  static Group::CatalogInterface::CatalogType catalog;
  return catalog;
}

WrapperBase & Group::registerWrapper( string const & name,
                                      std::unique_ptr< WrapperBase > wrapper )
{
  return *m_wrappers.insert( name, wrapper.release(), true );
}

void Group::deregisterWrapper( string const & name )
{
  GEOSX_ERROR_IF( !hasWrapper( name ), "Wrapper " << name << " doesn't exist." );
  m_wrappers.erase( name );
  m_conduitNode.remove( name );
}


void Group::resize( indexType const newSize )
{
  forWrappers( [newSize] ( WrapperBase & wrapper )
  {
    if( wrapper.sizedFromParent() == 1 )
    {
      wrapper.resize( newSize );
    }
  } );

  forSubGroups( [newSize] ( Group & subGroup )
  {
    if( subGroup.sizedFromParent() == 1 )
    {
      subGroup.resize( newSize );
    }
  } );

  m_size = newSize;
  if( m_size > m_capacity )
  {
    m_capacity = m_size;
  }
}

void Group::reserve( indexType const newSize )
{
  forWrappers( [newSize] ( WrapperBase & wrapper )
  {
    if( wrapper.sizedFromParent() == 1 )
    {
      wrapper.reserve( newSize );
    }
  } );

  forSubGroups( [newSize] ( Group & subGroup )
  {
    if( subGroup.sizedFromParent() == 1 )
    {
      subGroup.resize( newSize );
    }
  } );

  m_capacity = newSize;
}

string Group::getPath() const
{
  // In the Conduit node heirarchy everything begins with 'Problem', we should change it so that
  // the ProblemManager actually uses the root Conduit Node but that will require a full rebaseline.
  string const noProblem = getConduitNode().path().substr( sizeof( "Problem" ) -1 );
  return noProblem == "" ? "/" : noProblem;
}

void Group::processInputFileRecursive( xmlWrapper::xmlNode & targetNode )
{
  xmlWrapper::addIncludedXML( targetNode );

  // loop over the child nodes of the targetNode
  for( xmlWrapper::xmlNode childNode=targetNode.first_child(); childNode; childNode=childNode.next_sibling())
  {
    // Get the child tag and name
    string childName = childNode.attribute( "name" ).value();
    if( childName.empty())
    {
      childName = childNode.name();
    }

    // Create children
    Group * newChild = createChild( childNode.name(), childName );
    if( newChild == nullptr )
    {
      newChild = getGroupPointer( childName );
    }
    if( newChild != nullptr )
    {
      newChild->processInputFileRecursive( childNode );
    }
  }

  processInputFile( targetNode );
//  ProcessInputFile_PostProcess();
}

void Group::processInputFile( xmlWrapper::xmlNode const & targetNode )
{

  std::set< string > processedXmlNodes;
  for( std::pair< string const, WrapperBase * > & pair : m_wrappers )
  {
    if( pair.second->processInputFile( targetNode ) )
    {
      processedXmlNodes.insert( pair.first );
    }
  }

  for( xmlWrapper::xmlAttribute attribute=targetNode.first_attribute(); attribute; attribute = attribute.next_attribute() )
  {
    string const childName = attribute.name();
    if( childName != "name" && childName != "xmlns:xsi" && childName != "xsi:noNamespaceSchemaLocation" )
    {
      GEOSX_THROW_IF( processedXmlNodes.count( childName )==0,
                      "XML Node ("<<targetNode.name()<<") with attribute name=("<<
                      targetNode.attribute( "name" ).value()<<") contains child node named ("<<
                      childName<<") that is not read. Valid options are: \n" << dumpInputOptions()
                      + "\nFor more details, please refer to documentation at: \n"
                      + "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html \n",
                      InputError );
    }
  }
}

void Group::postProcessInputRecursive()
{
  for( auto const & subGroupIter : m_subGroups )
  {
    subGroupIter.second->postProcessInputRecursive();
  }
  postProcessInput();
}



void Group::registerDataOnMeshRecursive( Group & meshBodies )
{
  registerDataOnMesh( meshBodies );
  for( auto && subGroup : m_subGroups )
  {
    subGroup.second->registerDataOnMeshRecursive( meshBodies );
  }
}


Group * Group::createChild( string const & childKey, string const & childName )
{
  GEOSX_ERROR_IF( !(CatalogInterface::hasKeyName( childKey )),
                  "KeyName ("<<childKey<<") not found in Group::Catalog" );
  GEOSX_LOG_RANK_0( "Adding Object " << childKey<<" named "<< childName<<" from Group::Catalog." );
  return &registerGroup( childName,
                         CatalogInterface::factory( childKey, childName, this ) );
}


void Group::printDataHierarchy( integer const indent )
{
  for( auto & view : wrappers() )
  {
    GEOSX_LOG( string( indent, '\t' ) << view.second->getName() << ", " << LvArray::system::demangleType( view.second ) );
  }

  for( auto & group : m_subGroups )
  {
    GEOSX_LOG( string( indent, '\t' ) << group.first << ':' );
    group.second->printDataHierarchy( indent + 1 );
  }
}

string Group::dumpInputOptions() const
{
  string rval;

  bool writeHeader = true;
  for( auto const & wrapper : m_wrappers )
  {
    rval.append( wrapper.second->dumpInputOptions( writeHeader ) );
    writeHeader = false;
  }

  return rval;
}

void Group::deregisterGroup( string const & name )
{
  GEOSX_ERROR_IF( !hasGroup( name ), "Group " << name << " doesn't exist." );
  m_subGroups.erase( name );
  m_conduitNode.remove( name );
}

void Group::initializationOrder( string_array & order )
{
  for( auto & subGroupIter : m_subGroups )
  {
    order.emplace_back( subGroupIter.first );
  }
}

void Group::initialize()
{
  initializePreSubGroups();

  string_array initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initialize();
  }

  initializePostSubGroups();
}


void Group::initializePostInitialConditions()
{
  initializePostInitialConditionsPreSubGroups();

  string_array initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initializePostInitialConditions();
  }

  initializePostInitialConditionsPostSubGroups();
}

localIndex Group::packSize( string_array const & wrapperNames,
                            arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::PackSize( getName());

  packedSize += bufferOps::PackSize( string( "Wrappers" ));
  if( wrapperNames.size()==0 )
  {
    packedSize += bufferOps::PackSize( static_cast< int >(m_wrappers.size()) );
    for( auto const & wrapperPair : m_wrappers )
    {
      packedSize += bufferOps::PackSize( wrapperPair.first );
      if( packList.empty() )
      {
        packedSize += wrapperPair.second->packSize( true, onDevice, events );
      }
      else
      {
        packedSize += wrapperPair.second->packByIndexSize( packList, true, onDevice, events );
      }
    }
  }
  else
  {
    packedSize += bufferOps::PackSize( static_cast< int >(wrapperNames.size()) );
    for( auto const & wrapperName : wrapperNames )
    {
      WrapperBase const & wrapper = getWrapperBase( wrapperName );
      packedSize += bufferOps::PackSize( wrapperName );
      if( packList.empty() )
      {
        packedSize += wrapper.packSize( true, onDevice, events );
      }
      else
      {
        packedSize += wrapper.packByIndexSize( packList, true, onDevice, events );
      }
    }
  }
  if( recursive > 0 )
  {
    packedSize += bufferOps::PackSize( string( "SubGroups" ));
    packedSize += bufferOps::PackSize( m_subGroups.size() );
    for( auto const & keyGroupPair : m_subGroups )
    {
      packedSize += bufferOps::PackSize( keyGroupPair.first );
      packedSize += keyGroupPair.second->packSize( wrapperNames, packList, recursive, onDevice, events );
    }
  }

  return packedSize;
}


localIndex Group::packSize( string_array const & wrapperNames,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  arrayView1d< localIndex const > nullArray;
  return packSize( wrapperNames, nullArray, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        string_array const & wrapperNames,
                        arrayView1d< localIndex const > const & packList,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< true >( buffer, getName() );

  packedSize += bufferOps::Pack< true >( buffer, string( "Wrappers" ) );
  if( wrapperNames.size()==0 )
  {
    packedSize += bufferOps::Pack< true >( buffer, m_wrappers.size() );
    for( auto const & wrapperPair : m_wrappers )
    {
      packedSize += bufferOps::Pack< true >( buffer, wrapperPair.first );
      if( packList.empty() )
      {
        // invoke wrapper pack kernel
        packedSize += wrapperPair.second->pack( buffer, true, onDevice, events );
      }
      else
      {
        packedSize += wrapperPair.second->packByIndex( buffer, packList, true, onDevice, events );
      }
    }
  }
  else
  {
    packedSize += bufferOps::Pack< true >( buffer, wrapperNames.size() );
    for( auto const & wrapperName : wrapperNames )
    {
      WrapperBase const & wrapper = getWrapperBase( wrapperName );
      packedSize += bufferOps::Pack< true >( buffer, wrapperName );
      if( packList.empty() )
      {
        packedSize += wrapper.pack( buffer, true, onDevice, events );
      }
      else
      {
        packedSize += wrapper.packByIndex( buffer, packList, true, onDevice, events );
      }
    }
  }

  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< true >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< true >( buffer, m_subGroups.size() );
    for( auto const & keyGroupPair : m_subGroups )
    {
      packedSize += bufferOps::Pack< true >( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->pack( buffer, wrapperNames, packList, recursive, onDevice, events );
    }
  }

  return packedSize;
}

localIndex Group::pack( buffer_unit_type * & buffer,
                        string_array const & wrapperNames,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  arrayView1d< localIndex const > nullArray;
  return pack( buffer, wrapperNames, nullArray, recursive, onDevice, events );
}

localIndex Group::unpack( buffer_unit_type const * & buffer,
                          arrayView1d< localIndex > & packList,
                          integer const recursive,
                          bool onDevice,
                          parallelDeviceEvents & events )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOSX_ERROR_IF( groupName != getName(), "Group::Unpack(): group names do not match" );

  string wrappersLabel;
  unpackedSize += bufferOps::Unpack( buffer, wrappersLabel );
  GEOSX_ERROR_IF( wrappersLabel != "Wrappers", "Group::Unpack(): wrapper label incorrect" );

  localIndex numWrappers;
  unpackedSize += bufferOps::Unpack( buffer, numWrappers );
  for( localIndex a=0; a<numWrappers; ++a )
  {
    string wrapperName;
    unpackedSize += bufferOps::Unpack( buffer, wrapperName );
    getWrapperBase( wrapperName ).unpackByIndex( buffer, packList, true, onDevice, events );
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOSX_ERROR_IF( subGroups != "SubGroups", "Group::Unpack(): group names do not match" );

    decltype( m_subGroups.size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF( numSubGroups != m_subGroups.size(), "Group::Unpack(): incorrect number of subGroups" );

    for( auto const & index : m_subGroups )
    {
      GEOSX_UNUSED_VAR( index );
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += getGroup( subGroupName ).unpack( buffer, packList, recursive, onDevice, events );
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

  forWrappers( [] ( WrapperBase & wrapper )
  {
    wrapper.registerToWrite();
  } );

  m_conduitNode[ "__size__" ].set( m_size );

  forSubGroups( []( Group & subGroup )
  {
    subGroup.prepareToWrite();
  } );
}


void Group::finishWriting()
{
  if( getRestartFlags() == RestartFlags::NO_WRITE )
  {
    return;
  }

  forWrappers( [] ( WrapperBase & wrapper )
  {
    wrapper.finishWriting();
  } );

  forSubGroups( []( Group & subGroup )
  {
    subGroup.finishWriting();
  } );
}


void Group::loadFromConduit()
{
  if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
  {
    return;
  }

  m_size = m_conduitNode.fetch_child( "__size__" ).value();
  localIndex const groupSize = m_size;

  forWrappers( [&]( WrapperBase & wrapper )
  {
    if( !( wrapper.loadFromConduit()) )
    {
      if( wrapper.sizedFromParent() == 1 )
      {
        wrapper.resize( groupSize );
      }
    }
  } );

  forSubGroups( []( Group & subGroup )
  {
    subGroup.loadFromConduit();
  } );
}

void Group::postRestartInitializationRecursive()
{
  forSubGroups( [&]( Group & subGroup )
  {
    subGroup.postRestartInitializationRecursive();
  } );

  postRestartInitialization();
}

void Group::enableLogLevelInput()
{
  string const logLevelString = "logLevel";

  registerWrapper( logLevelString, &m_logLevel ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Log level" );
}

Group const & Group::getBaseGroupByPath( string const & path ) const
{
  Group const * currentGroup = this;
  string::size_type previousPosition = 0;

  if( path[ 0 ] == '/' )
  {
    bool foundTarget = false;
    for( int i=0; i<1000; ++i )
    {
      if( currentGroup->m_parent != nullptr )
      {
        currentGroup = currentGroup->m_parent;
      }
      else
      {
        foundTarget = true;
        previousPosition = 1;
        break;
      }
    }
    GEOSX_ERROR_IF( !foundTarget,
                    "Could not find the specified path from the starting group." );
  }

  string::size_type currentPosition;
  do
  {
    currentPosition = path.find( '/', previousPosition );
    string const curGroupName = path.substr( previousPosition, currentPosition - previousPosition );

    previousPosition = currentPosition + 1;

    if( curGroupName == "" || curGroupName == "." || curGroupName==currentGroup->m_name )
    {
      continue;
    }
    else if( curGroupName == ".." )
    {
      currentGroup = &this->getParent();
    }
    else
    {
      currentGroup = &currentGroup->getGroup( curGroupName );
    }
  }
  while( currentPosition != string::npos );

  return *currentGroup;
}

} /* end namespace dataRepository */
} /* end namespace geosx  */
