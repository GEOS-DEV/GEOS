/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
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
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"

#include <numeric>
#include <unordered_set>

#if defined(GEOSX_USE_PYGEOSX)
#include "python/PyGroupType.hpp"
#endif

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

WrapperBase & Group::registerWrapper( std::unique_ptr< WrapperBase > wrapper )
{
  // Extract `wrapperName` first to prevent from UB call order in the `insert` call.
  string const wrapperName = wrapper->getName();
  return *m_wrappers.insert( wrapperName, wrapper.release(), true );
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
  string const noProblem = getConduitNode().path().substr( std::strlen( dataRepository::keys::ProblemManager ) - 1 );
  return noProblem.empty() ? "/" : noProblem;
}

void Group::processInputFileRecursive( xmlWrapper::xmlNode & targetNode )
{
  xmlWrapper::addIncludedXML( targetNode );

  // Handle the case where the node was imported from a different input file
  // Set the path prefix to make sure all relative Path variables are interpreted correctly
  string const oldPrefix = Path::pathPrefix();
  xmlWrapper::xmlAttribute filePath = targetNode.attribute( xmlWrapper::filePathString );
  if( filePath )
  {
    Path::pathPrefix() = splitPath( filePath.value() ).first;
    targetNode.remove_attribute( filePath );
  }

  // Loop over the child nodes of the targetNode
  array1d< string > childNames;
  for( xmlWrapper::xmlNode childNode : targetNode.children() )
  {
    // Get the child tag and name
    string childName = childNode.attribute( "name" ).value();
    if( childName.empty() )
    {
      childName = childNode.name();
    }
    else
    {
      // Make sure child names are not duplicated
      GEOSX_ERROR_IF( std::find( childNames.begin(), childNames.end(), childName ) != childNames.end(),
                      GEOSX_FMT( "Error: An XML block cannot contain children with duplicated names ({}/{}). ",
                                 getPath(), childName ) );
      childNames.emplace_back( childName );
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

  // Restore original prefix once the node is processed
  Path::pathPrefix() = oldPrefix;
}

void Group::processInputFile( xmlWrapper::xmlNode const & targetNode )
{
  std::set< string > processedAttributes;
  for( std::pair< string const, WrapperBase * > & pair : m_wrappers )
  {
    if( pair.second->processInputFile( targetNode ) )
    {
      processedAttributes.insert( pair.first );
    }
  }

  for( xmlWrapper::xmlAttribute attribute : targetNode.attributes() )
  {
    string const attributeName = attribute.name();
    if( attributeName != "name" && attributeName != "xmlns:xsi" && attributeName != "xsi:noNamespaceSchemaLocation" )
    {
      GEOSX_THROW_IF( processedAttributes.count( attributeName ) == 0,
                      GEOSX_FMT( "XML Node '{}' with name='{}' contains unused attribute '{}'.\n"
                                 "Valid attributes are:\n{}\nFor more details, please refer to documentation at:\n"
                                 "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html",
                                 targetNode.path(), targetNode.attribute( "name" ).value(), attributeName, dumpInputOptions() ),
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

void Group::printMemoryAllocation( integer const indent, real64 const threshold )
{
  // static flag to keep track of whether or not a subgroup is the last one, and thus
  // visually will terminate the tree.
  static bool terminateBranch[64]{};

  // keep track of the address of the groups that have been printed. Not the best
  // way to do this, but doens't require infrastructure changes.
  static std::unordered_set< void * > groupsPrinted;

  // initialize the static variables at the head of the tree.
  if( indent==0 )
  {
    terminateBranch[0] = true;
    for( int i=1; i<64; ++i )
    {
      terminateBranch[i] = false;
    }
    groupsPrinted.clear();
  }

  // store the local allocations for each wrapper
  std::vector< double > localAllocations;

  // use the first index for the summation of all Wrappers in a Group
  localAllocations.emplace_back( 0 );
  for( auto & view : wrappers() )
  {
    double const bytesAllocated = view.second->bytesAllocated();
    localAllocations.emplace_back( bytesAllocated );
    localAllocations[0] += bytesAllocated;
  }

  int const numRanks = MpiWrapper::commSize();
  int const numValues = localAllocations.size();

  // storage for the gathered values of localAllocation
  std::vector< double > globalAllocations;
  if( MpiWrapper::commRank()==0 )
  {
    globalAllocations.resize( numRanks * numValues );
  }

  MpiWrapper::gather( localAllocations.data(),
                      numValues,
                      globalAllocations.data(),
                      numValues,
                      0,
                      MPI_COMM_GEOSX );

  // reduce data across ranks (min, max, sum)
  if( MpiWrapper::commRank()==0 )
  {
    array2d< double > allocationReductions( numValues, 3 );
    for( int a=0; a<numValues; ++a )
    {
      allocationReductions( a, 0 ) = 1e99;
      allocationReductions( a, 1 ) = 0;
      allocationReductions( a, 1 ) = 0;
      for( int b=0; b<numRanks; ++b )
      {
        int const recvIndex = a + b * numValues;
        double const value = globalAllocations[recvIndex];
        allocationReductions( a, 0 ) = std::min( allocationReductions( a, 0 ), value );
        allocationReductions( a, 1 ) = std::max( allocationReductions( a, 1 ), value );
        allocationReductions( a, 2 ) += value;
      }
    }

    // scope for the output of groups
    {
      double const * const groupAllocations = allocationReductions[0];

      if( indent==0 )
      {
        //                          1         2         3         4         5         6         7         8         9        10        11
        //        12
        //                 123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        GEOSX_LOG_RANK_0( "************************************************************************************************************************" );
        GEOSX_LOG_RANK_0( " Data repository memory allocations " );
        GEOSX_LOG_RANK_0( "                                                                                       min          max          sum    " );
        GEOSX_LOG_RANK_0( "                                                                                   -----------  -----------  -----------" );
      }


      string outputLine;
      int indentChars = 0;
      for( int i=0; i<indent; ++i )
      {
        outputLine += terminateBranch[i]==true ? "   " : "|  ";
        indentChars += 3;
      }
      // put a indent between groups in the tree
      GEOSX_LOG_RANK_0( outputLine.c_str()<<"|" );
      indentChars += 1;
      // only allocation data if it is above the threshold...
      if( groupAllocations[0] >= threshold )
      {
        indentChars += 3;
        outputLine += "|--{:.<" + std::to_string( 83-indentChars ) + "} {:>9s}    {:>9s}    {:>9s}";
        GEOSX_LOG_RANK_0( GEOSX_FMT( outputLine.c_str(),
                                     "[" + getName() + "]",
                                     stringutilities::toMetricPrefixString( groupAllocations[0] ) + 'B',
                                     stringutilities::toMetricPrefixString( groupAllocations[1] ) + 'B',
                                     stringutilities::toMetricPrefixString( groupAllocations[2] ) + 'B' ) );
      }
      else  // ...but we still need to output the group name to have a valid tree.
      {
        outputLine += "|--[{:<}]";
        GEOSX_LOG_RANK_0( GEOSX_FMT( outputLine.c_str(),
                                     getName() ) );
      }
    }

    // output the views
    {
      localIndex viewCount = 1;
      for( auto & view : wrappers() )
      {
        if( allocationReductions( viewCount, 0 ) > threshold )
        {
          string outputLine;
          int indentChars = 0;
          for( int i=0; i<=indent; ++i )
          {
            outputLine += terminateBranch[i]==true ? "   " : "|  ";
            indentChars += 3;
          }
          indentChars += 5;
          outputLine += "| - {:.<" + std::to_string( 83-indentChars ) + "} {:>9s}    {:>9s}    {:>9s}";
          GEOSX_LOG_RANK_0( GEOSX_FMT( outputLine.c_str(),
                                       view.second->getName(),
                                       stringutilities::toMetricPrefixString( allocationReductions( viewCount, 0 ) ) + 'B',
                                       stringutilities::toMetricPrefixString( allocationReductions( viewCount, 1 ) ) + 'B',
                                       stringutilities::toMetricPrefixString( allocationReductions( viewCount, 2 ) ) + 'B' ) );
        }
        ++viewCount;
      }
    }
  }

  // process subgroups
  localIndex const numSubGroups = m_subGroups.size();
  localIndex groupCounter = 0;
  for( auto & group : m_subGroups )
  {
    // set flag to indicate whether or not this is the last subgroup so that
    // the ascii art tree may be constructed properly
    terminateBranch[indent+1] = ++groupCounter==numSubGroups ? true : false;

    // check to see that the group hasen't been printed yet.
    if( groupsPrinted.count( group.second ) == 0 )
    {
      // insert the address of the Group into the map that holds printed groups.
      groupsPrinted.insert( group.second );
      group.second->printMemoryAllocation( indent + 1, threshold );
    }
  }
  if( indent==0 )
  {
    //         1         2         3         4         5         6         7         8         9        10
    //1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    GEOSX_LOG_RANK_0( "****************************************************************************************************" );
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

void Group::initialize_postMeshGeneration()
{
  array1d< string > initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initialize_postMeshGeneration();
  }
}


void Group::initialize()
{
  initializePreSubGroups();

  array1d< string > initOrder;
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

  array1d< string > initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initializePostInitialConditions();
  }

  initializePostInitialConditionsPostSubGroups();
}

template< bool DO_PACKING >
localIndex Group::packImpl( buffer_unit_type * & buffer,
                            array1d< string > const & wrapperNames,
                            arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, getName() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "Wrappers" ) );

  // `wrappers` are considered for packing if they match the size of this Group instance.
  // A way to check this is to check the sufficient (but not necessary...) condition `wrapper.sizedFromParent()`.
  std::vector< WrapperBase const * > wrappers;
  for( string const & wrapperName: wrapperNames )
  {
    if( hasWrapper( wrapperName ) )
    {
      WrapperBase const & wrapper = getWrapperBase( wrapperName );

      if( wrapper.sizedFromParent() )
      {
        wrappers.push_back( &wrapper );
      }
    }
    else
    {
      GEOSX_ERROR( "Wrapper " << wrapperName << " not found in Group " << getName() << "." );
    }
  }

  // Now we pack the `wrappers`.
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, LvArray::integerConversion< localIndex >( wrappers.size() ) );
  for( WrapperBase const * wrapper: wrappers )
  {
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, wrapper->getName() );
    if( packList.empty() )
    {
      packedSize += wrapper->pack< DO_PACKING >( buffer, true, onDevice, events );
    }
    else
    {
      packedSize += wrapper->packByIndex< DO_PACKING >( buffer, packList, true, onDevice, events );
    }
  }

  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_subGroups.size() );
    for( auto const & keyGroupPair : m_subGroups )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->packImpl< DO_PACKING >( buffer, wrapperNames, packList, recursive, onDevice, events );
    }
  }

  return packedSize;
}

localIndex Group::packSize( array1d< string > const & wrapperNames,
                            arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  buffer_unit_type * dummy;
  return this->packImpl< false >( dummy, wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::packSize( arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  std::vector< string > const tmp = mapKeys( m_wrappers );
  array1d< string > wrapperNames;
  wrapperNames.insert( 0, tmp.begin(), tmp.end() );
  return this->packSize( wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::packSize( array1d< string > const & wrapperNames,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  arrayView1d< localIndex const > nullArray;
  return packSize( wrapperNames, nullArray, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        array1d< string > const & wrapperNames,
                        arrayView1d< localIndex const > const & packList,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  return this->packImpl< true >( buffer, wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        arrayView1d< localIndex const > const & packList,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  std::vector< string > const tmp = mapKeys( m_wrappers );
  array1d< string > wrapperNames;
  wrapperNames.insert( 0, tmp.begin(), tmp.end() );
  return this->pack( buffer, wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        array1d< string > const & wrapperNames,
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
  GEOSX_ERROR_IF( groupName != getName(), "Group::unpack(): group names do not match" );

  string wrappersLabel;
  unpackedSize += bufferOps::Unpack( buffer, wrappersLabel );
  GEOSX_ERROR_IF( wrappersLabel != "Wrappers", "Group::unpack(): wrapper label incorrect" );

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
    GEOSX_ERROR_IF( subGroups != "SubGroups", "Group::unpack(): group names do not match" );

    decltype( m_subGroups.size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOSX_ERROR_IF( numSubGroups != m_subGroups.size(), "Group::unpack(): incorrect number of subGroups" );

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

  m_size = m_conduitNode.child( "__size__" ).value();
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

localIndex Group::getSubGroupIndex( keyType const & key ) const
{
  return getSubGroups().getIndex( key );
}

#if defined(GEOSX_USE_PYGEOSX)
PyTypeObject * Group::getPythonType() const
{ return geosx::python::getPyGroupType(); }
#endif

} /* end namespace dataRepository */
} /* end namespace geosx  */
