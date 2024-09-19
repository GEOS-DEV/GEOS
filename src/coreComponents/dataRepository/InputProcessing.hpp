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

/**
 * @file InputProcessing.hpp
 */

#ifndef GEOS_DATAREPOSITORY_INPUTPROCESSING_HPP_
#define GEOS_DATAREPOSITORY_INPUTPROCESSING_HPP_

#include "codingUtilities/ConstexprConstructs.hpp"
#include "InputExtension.hpp"
#include "InputExtensionGroup.hpp"
#include "InputParsing.hpp"
#include "TerseInputRegistry.hpp"
#include "Group.hpp"

#include <vector>

namespace geos
{

namespace dataRepository
{

namespace inputProcessing
{

struct NameAttrAsIdentity
{
  template < typename NodeType >
  string operator()( NodeType & docNode )
  {
    if ( docNode.attribute("name").empty() )
    {
      return string(docNode.name());
    }
    else
    {
      return docNode.attribute("name").value();
    }
  }
};

template < typename DocumentType, typename DataNode >
class InputProcessingPhase
{
  static_assert( inputExtension::is_valid_document< DocumentType >::value, "Provided DocumentType doesn't expose required document interfaces!" );
public:
  using document_type = DocumentType;
  using data_node = DataNode;
  using document_node = typename document_type::node_type;
  using document_node_pos = typename document_type::node_pos_type;
  using document_attribute = typename document_type::attribute_type;

  InputProcessingPhase( document_type & document, std::set< string > & mergableNodes, std::set< string > & allowedDeviations )
  : m_document( document ),
    m_inputExtender( mergableNodes ),
    m_allowedDeviations( allowedDeviations )
  {}

  virtual string const name() const = 0;

  void execute( data_node & dataRoot, document_node & docRoot )
  {
    // Map to track corresponding data nodes
    std::map< document_node, data_node * > doc2data;
    // Map to track the names of known child nodes for each data_node
    std::map< string, std::set< string > > childNamesMap;

    doc2data[ docRoot ] = &dataRoot;

    DynamicTreeIteration< NameAttrAsIdentity >::processTree( docRoot, [ this, &dataRoot, &docRoot, &doc2data, &childNamesMap ]( document_node & docNode )
    {
      string const oldPrefix = this->updatePrefixPath( docNode );
      string const dataNodeName = this->getDataNodeName( docNode );

      data_node * dataParent = nullptr;
      data_node * dataNode = doc2data[ docNode ];

      if ( docNode != docRoot )
      {
        if ( m_allowedDeviations.find( docNode.name() ) != m_allowedDeviations.end() )
        {
          std::cout << "deviation node " << docNode.name() << " ecountered, skipping" << std::endl;
          return;
        }
        std::cout << "looking for known ancestor of " << docNode.name() << std::endl;

        document_node currentDocNode = docNode.parent();
        std::vector< string > pathToKnownAncestor;

        // Walk up the tree to find the closest mapped ancestor
        while ( doc2data.find( currentDocNode ) == doc2data.end() && currentDocNode != docRoot )
        {
          string name = currentDocNode.attribute( "name" ).value();
          if( name == "" )
          {
            name = currentDocNode.name();
          }
          std::cout << "  unmapped: " << name << std::endl;
          if ( m_allowedDeviations.find( name ) != m_allowedDeviations.end() )
          {
            std::cout << "  deviation found while searching for known ancestor, skipping " << docNode.name() << std::endl;
            return;
          }
          pathToKnownAncestor.push_back( name );
          currentDocNode = currentDocNode.parent();
        }

        std::cout << "  found: " << currentDocNode.name() << std::endl;

        std::cout << "walking back to parent of data node" << std::endl;
        // walk back down the tree to the parent of the current node
        data_node * currentDataNode = doc2data[ currentDocNode ];
        for( auto riter = pathToKnownAncestor.rbegin(); riter != pathToKnownAncestor.rend(); ++riter )
        {
          std::cout << "  walking to " << *riter << std::endl;
          currentDataNode = &dereference( currentDataNode->getGroupPointer( *riter ) );
        }
        dataParent = currentDataNode;
        std::cout << "  found parent: " << dataParent->getName() << std::endl;

        // process the current node
        validateUniqueChildren( *dataParent, docNode, dataNodeName, childNamesMap );
        dataNode = this->retrieveOrCreateDataNode( dataParent, dataNodeName, docNode.name() );
        if ( dataNode != dataParent ) // this is a bit hacky, don't really like it
        {
          doc2data[ docNode ] = dataNode;
          this->processNode( dataNode, docNode );
        }
      }
      else
      {
        std::cout << " root " << docNode.name() << std::endl;
        validateUniqueChildren( dataRoot, docNode, dataNodeName, childNamesMap );
        this->processNode( &dataRoot, docNode );
      }
      restorePathPrefix( oldPrefix );
    } );
  }

protected:
  virtual void processNode( data_node *, document_node & ) = 0;
  virtual data_node * retrieveOrCreateDataNode( data_node * parent, string const & dataNodeName, string const & docNodeName ) = 0;

  void validateUniqueChildren( data_node & dataParent, document_node & docNode, string const & dataNodeName, std::map< string, std::set< string > > & childNamesMap )
  {
    auto & knownChildNames = childNamesMap[ dataParent.getName() ];
    document_node_pos docNodePos = m_document.getNodePosition( docNode );
    GEOS_ERROR_IF( std::find( knownChildNames.begin(), knownChildNames.end(), dataNodeName ) != knownChildNames.end(),
                  GEOS_FMT( "Error: An Input block cannot contain children with duplicated names.\n"
                            "Error detected at node {} with name = {} ({}:l.{})",
                            docNode.path(), dataNodeName, m_document.getFilePath(), docNodePos.line ) );
    knownChildNames.insert( dataNodeName ); // store the valid child name
  }

  string getDataNodeName( document_node & docNode )
  {
    string dataNodeName;
    inputParsing::readAttributeAsType< DocumentType >( dataNodeName,
                                                       "name",
                                                       rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupName ),
                                                       docNode,
                                                       string( "" ) );
    if ( dataNodeName.empty() )
    {
      dataNodeName = docNode.name();
    }
    return dataNodeName;
  }

  string updatePrefixPath( document_node & docNode )
  {
    string oldPrefix = std::string( Path::getPathPrefix() );
    document_attribute filePath = docNode.attribute( inputParsing::filePathString );
    if ( filePath )
    {
      Path::setPathPrefix( getAbsolutePath( splitPath( filePath.value() ).first ) );
    }
    return oldPrefix;
  }

  void restorePathPrefix( const string & oldPrefix )
  {
    Path::setPathPrefix( oldPrefix );
  }

  document_type & m_document;
  inputExtension::InputExtender< document_node > m_inputExtender;
  std::set< string > const m_allowedDeviations;
};

template < typename DocumentType, typename DataNode, typename... InputPhase >
class InputProcessor
{
public:
  using document_type = DocumentType;
  using data_node = DataNode;
  using document_node = typename document_type::node_type;

  static_assert((std::is_base_of< InputProcessingPhase< document_type, data_node >, InputPhase >::value && ...),
                "All Phases must be subclasses of InputProcessingPhase< DocumentType, DataNode >");

  explicit InputProcessor( document_type & document, std::set< string > & mergableNodes, std::set< string > & allowedDeviations )
    : m_phases( makePhases( document, mergableNodes, allowedDeviations, std::index_sequence_for<InputPhase...>{} ) )
  { }

  void execute( data_node & dataRoot, document_node & documentRoot )
  {
    compileTime::static_for< 0, getNumberOfPhases() >( [this, &dataRoot, &documentRoot] ( auto ii )
    {
      auto constexpr idx = decltype(ii)::value;
      this->template executePhase< idx >( dataRoot, documentRoot );
    } );
  }

  // Execute specific phase by index
  template < std::size_t I >
  void executePhase( data_node & dataRoot, document_node & documentRoot )
  {
    auto & phase = std::get< I >( m_phases );
    GEOS_LOG_RANK_0( GEOS_FMT( "-- Entering input processing phase '{}'", phase.name() ) );
    phase.execute( dataRoot, documentRoot );
    GEOS_LOG_RANK_0( GEOS_FMT( "-- Finished input processing phase '{}'", phase.name() ) );
  }

  // Get the number of registered phases
  static constexpr std::size_t getNumberOfPhases()
  {
    return sizeof...(InputPhase);
  }

private:
  template<std::size_t... Is>
  static std::tuple< InputPhase... > makePhases( document_type & document, std::set< string > & mergableNodes, std::set< string > & allowedDeviations, std::index_sequence<Is...> )
  {
    return std::make_tuple( InputPhase( document, mergableNodes, allowedDeviations )... );
  }

  std::tuple< InputPhase... > m_phases;
};

// Declaration: handles the initial processing of the document, applying input extension rules from internal Groups as encountered (e.g. deprecation)
//              at the end of this process all Groups needed for the simulation should be created
template < typename DocumentType, typename DataNode >
class Declaration : public InputProcessingPhase< DocumentType, DataNode >
{
public:
  using super = InputProcessingPhase< DocumentType, DataNode >;
  using super::super;
  using typename super::document_type;
  using typename super::data_node;
  using typename super::document_node;

  virtual string const name() const override { return "Declaration"; };

protected:
  using super::m_inputExtender;

  virtual data_node * retrieveOrCreateDataNode( data_node * parent, string const & dataNodeName, string const & docNodeName ) override
  {
    decltype(auto) child = parent->getGroupPointer( dataNodeName );
    if( child == nullptr )
    {
      child = parent->createChild( docNodeName, dataNodeName );
    }
    return child;
  }

  // this does traverse the internal node structure as we create the nodes as we encounter them
  virtual void processNode( data_node * dataNode, document_node & docNode ) override
  {
    auto extensions = dataNode->getInputExtensionRules( m_inputExtender );
    inputExtension::StaticParsingContext< document_node > context { docNode };
    m_inputExtender.applyExtensions( context, extensions );
  }
};

template < typename DocumentType, typename DataNode >
class TerseSyntax : public InputProcessingPhase< DocumentType, DataNode >
{
public:
  using super = InputProcessingPhase< DocumentType, DataNode >;
  using typename super::document_type;
  using typename super::data_node;
  using typename super::document_node;

  explicit TerseSyntax( document_type & document, std::set< string > & mergableNodes, std::set< string > & allowedDeviations )
    : super( document, mergableNodes, allowedDeviations ),
      m_terseExtensionRules( )
  {
    m_terseExtensionRules = inputExtension::TerseInputRegistry< document_node >::getInputExtensionRules();
  }

  virtual string const name( ) const override { return  "TerseSyntax"; }

protected:
  using super::m_inputExtender;

  virtual data_node * retrieveOrCreateDataNode( data_node * parent, string const & GEOS_UNUSED_PARAM( dataNodeName ), string const & GEOS_UNUSED_PARAM( docNodeName ) ) override
  {
    return parent;
  }

  virtual void processNode( data_node * GEOS_UNUSED_PARAM( _ ), document_node & docNode ) override
  {
    inputExtension::StaticParsingContext< document_node > context { docNode };
    auto extensions = softMapLookup( m_terseExtensionRules, docNode.name(), std::vector< inputExtension::Rule< document_node > > {} );
    m_inputExtender.applyExtensions( context, extensions );
  }
private:
  typename inputExtension::TerseInputRegistry< document_node >::map_type m_terseExtensionRules;
};

// Definition: Handles the final processing of the document, populating internal Wrapper values for the simulation.
//             No input extension takes place in this phase.
template < typename DocumentType, typename DataNode, typename DataAttribute >
class Definition : public InputProcessingPhase< DocumentType, DataNode >
{
public:
  using super = InputProcessingPhase< DocumentType, DataNode >;
  using typename super::document_type;
  using typename super::data_node;
  using internal_attribute_type = DataAttribute;
  using typename super::document_node;
  using document_node_pos = typename document_type::node_pos_type;
  using document_attribute = typename document_type::attribute_type;
  using document_attribute_pos = typename document_type::attribute_pos_type;

  explicit Definition( document_type & document, std::set< string > & mergableNodes, std::set< string > & allowedDeviations )
    : super( document, mergableNodes, allowedDeviations ),
      m_document(document)
  {}

  virtual string const name() const override { return "Definition"; };

protected:

  virtual data_node * retrieveOrCreateDataNode( data_node * parent, string const & dataNodeName, string const & docNodeName ) override
  {
    // Doing this is a "hack" to allow the TerseSyntax phase taking place after the Declaration phase to not instantiate any groups
    // and just focus on input expansion / modification of the input document structure, there are other solutions, but this
    // is probably more performant.
    decltype(auto) child = parent->getGroupPointer( dataNodeName );
    if( child == nullptr )
    {
      child = parent->createChild( docNodeName, dataNodeName );
    }
    return child;
  }

  virtual void processNode( data_node * dataNode, document_node & docNode ) override
  {
    document_node_pos docNodePos = m_document.getNodePosition( docNode );
    dataNode->registerDataContext( std::make_unique< DataFileContext< document_node, document_node_pos, document_attribute, document_attribute_pos > >( docNode, docNodePos ) );
    processAttributes( *dataNode, docNode, docNodePos );
  }

  void processAttributes( data_node & dataNode, document_node const & docNode, document_node_pos const & docNodePos )
  {
    std::set< string > processedAttributes;
    dataNode.forWrappers( [&] ( auto & dataAttribute )
    {
      if ( dataAttribute.initFromInput( docNode, docNodePos ) )
      {
        processedAttributes.insert( dataAttribute.getName() );
      }
    });

    // Check for unused attributes and throw if any are found
    for ( auto attribute : docNode.attributes() )
    {
      string attributeName = attribute.name();
      if ( ! inputParsing::isDocMetadataAttribute< DocumentType >( attributeName ) && processedAttributes.count( attributeName ) == 0)
      {
        GEOS_THROW_IF( processedAttributes.count( attributeName ) == 0,
                       GEOS_FMT( "Document Node at '{}' with name={} contains unused attribute '{}'.\n"
                                 "Valid attributes are:\n{}\nFor more details, please refer to documentation at:\n"
                                 "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html",
                                 docNode.path(), dataNode.getDataContext().toString(), attributeName,
                                 dataNode.dumpInputOptions() ),
                       InputError );
      }
    }
  }

private:
  document_type const & m_document;
};

using PreprocessOnly = InputProcessor< xmlWrapper::xmlDocument,
                                       Group,
                                       TerseSyntax< xmlWrapper::xmlDocument, Group >,
                                       Declaration< xmlWrapper::xmlDocument, Group >,
                                       TerseSyntax< xmlWrapper::xmlDocument, Group > >;

using AllProcessingPhases = InputProcessor< xmlWrapper::xmlDocument,
                                            Group,
                                            TerseSyntax< xmlWrapper::xmlDocument, Group >,
                                            Declaration< xmlWrapper::xmlDocument, Group >,
                                            TerseSyntax< xmlWrapper::xmlDocument, Group >,
                                            Definition< xmlWrapper::xmlDocument, Group, WrapperBase > >;


} // namespace inputProcessing

} // namespace dataRepository

} // namespace geos

#endif