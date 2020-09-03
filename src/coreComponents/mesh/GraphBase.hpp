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

/**
 * @file GraphBase.hpp
 */

#ifndef GEOSX_MESH_GRAPHBASE_HPP
#define GEOSX_MESH_GRAPHBASE_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/GraphEdge.hpp"
#include "mesh/GraphVertex.hpp"




namespace geosx
{

namespace dataRepository
{}

class NodeManager;
class DomainPartition;

/**
 *  @class GraphBase
 *  @brief The GraphBase class provides an abstract base class implementation for different graph types.
 *	   The GraphBase is the Group specialization for different type of graph handling.
 */
class GraphBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for Graph base class.
   * @param[in] name of the Graph object
   * @param[in] parent the parent Group pointer for the Graph object
   */
  explicit GraphBase( std::string const & name,
                              Group * const parent );

  /**
   * @brief Destructor for GraphBase
   */
  virtual ~GraphBase();

  /**
   * @brief Function called to construct the graphs
  */ 
  virtual void GenerateGraph() = 0;

  /**
  * @brief Function called to partition of the graph using the partition of given mesh
  * @param[in] meshLevel the mesh already partitioned to copy from
  */
  virtual void PartitionGraph(const MeshLevel & meshLevel) = 0;

  /**
  * @brief Getter of edge list from the graph
  * @return edge list from the graph
  */ 
  std::vector<GraphEdge*> getEdges() { return m_edges; }

  /**
  * @brief Constant getter of edge list from the graph
  * @return constant edge list from the graph
  */ 
  std::vector<GraphEdge*> getEdges() const { return m_edges; }

  /**
  * @brief Getter of vertex list from the graph
  * @return vertex list from the graph
  */
  std::vector<GraphVertex*> getVertices() { return m_vertices; }

  /**
  * @brief Constant getter of vertex list from the graph
  * @return constant vertex list from the graph
  */
  std::vector<GraphVertex*> getVertices() const { return m_vertices; }

  /**
  * @brief Getter of association map between vertex and edge from the graph
  * @return association map between vertex and edge from the graph
  */
  std::map<GraphVertex*, std::vector<GraphEdge*>> getVertexWithEdgesMap() { return m_vertexWithEdgesMap; }

  /**
  * @brief Constant getter of association map between vertex and edge from the graph
  * @return constant association map between vertex and edge from the graph
  */
  std::map<GraphVertex*, std::vector<GraphEdge*>> getVertexWithEdgesMap() const { return m_vertexWithEdgesMap; }

  /**
   * @brief Return the name of the Graph in object catalog.
   * @return string that contains the catalog name of the Graph
   */
  static string CatalogName() { return "GraphBase"; }


  /// using alias for templated Catalog graph type
  using CatalogInterface = dataRepository::CatalogInterface< GraphBase, std::string const &, Group * const >;

/**
 * @brief Accessor for the singleton Catalog object
 * @return a static reference to the Catalog object
 *
 */
  static CatalogInterface::CatalogType & GetCatalog();
protected:
  std::vector<GraphEdge*> m_edges;
  std::vector<GraphVertex*> m_vertices;
  std::map<GraphVertex*, std::vector<GraphEdge*>> m_vertexWithEdgesMap;
};
}

#endif /* GEOSX_MESH_GRAPHBASE_ */
