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
 * @file xmlInputExtension.hpp
 */

#ifndef GEOS_DATAREPOSITORY_XMLINPUTEXTENSION_HPP_
#define GEOS_DATAREPOSITORY_XMLINPUTEXTENSION_HPP_

namespace geos
{

namespace xmlWrapper
{

using namespace dataRepository;

template < typename DocumentType >
struct InputExtensionConverter
{
public:
  using document_type = DocumentType;
  using node_type = typename document_type::node_type;
  using attribute_type = typename document_type::attribute_type;

  // Function to load an XML document and recursively convert its contents to InputExtension nodes
  static std::vector< inputExtension::Node< node_type > > convert( node_type & root )
  {
    // This won't actually work since we have no way of knowing who the actual parent is at any time during the static traversal
    std::map< node_type, inputExtension::Node< node_type > > nodeMap;
    StaticTreeIteration<>::processTree( root, [ &nodeMap ]( node_type & node )
    {
      inputExtension::Node< node_type > extensionNode;
      extensionNode.identifier = node.name();
      for( const attribute_type & attribute : node.attributes() )
      {
        extensionNode.attributes.push_back( inputExtension::mapAttribute::newAttribute( attribute.name(), attribute.value() ) );
      }
      // Check if the node has a parent and map it appropriately
      if ( ! node.parent().empty() )
      {
        nodeMap[node.parent()].subNodes.push_back(extensionNode);
      }
      nodeMap[node] = extensionNode;
    } );
    if (nodeMap.find( root ) != nodeMap.end())
    {
      return { nodeMap[ root ] }; // Wrap single root node into vector
    }
    return { };
  }
};

}

}

#endif