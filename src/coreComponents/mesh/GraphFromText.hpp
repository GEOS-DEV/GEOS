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
 * @file GraphFromText.hpp
 */

#ifndef GEOSX_MESH_GRAPHFROMTEXT_HPP_
#define GEOSX_MESH_GRAPHFROMTEXT_HPP_

#include "dataRepository/Group.hpp"
#include "common/Path.hpp"
#include "mesh/GraphBase.hpp"
#include "mesh/GraphEdge.hpp"
#include "mesh/GraphVertex.hpp"
#include "mesh/MeshLevel.hpp"


namespace geosx
{

/**
 * @class GraphFromText
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class GraphFromText : public GraphBase
{
public:

  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  GraphFromText( const std::string & name,
                 Group * const parent );

  /// Destructor
  virtual ~GraphFromText() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string CatalogName() { return "GraphFromText"; }
  
  /**
  * @brief Creation and registration of a new edge
  * @param[in] ind index of the edge
  * @param[in] v1 first vertex on the edge
  * @param[in] v2 second vertex on the edge
  * @param[in] transm transmissibility associated with the edge
  */ 
  void AddEdge(localIndex ind, std::shared_ptr<GraphVertex> v1, std::shared_ptr<GraphVertex> v2, real64 transm);

  /**
  * @brief Deletion of an edge using its index
  * @param[in] ind index of the dege to remove
  */ 
  void RemoveEdge(localIndex ind);

  /**
  * @brief Deletion of a vertex using its indexes
  * @param[in] er Region index
  * @param[in] esr Surregion index
  * @param[in] ei Global vertex index
  */ 
  void RemoveVertex(localIndex er, localIndex esr, globalIndex ei);

  /**
  * @brief Deletion of a vertex using itself as parameter
  * @param[in] vertex Vertex to delete
  */
  void RemoveVertex(std::shared_ptr<GraphVertex> vertex);
  
  /**
  * @brief Creation and registration of a new edge
  * @param[in] er Region index
  * @param[in] esr Surregion index
  * @param[in] ei Global vertex index
  */ 
  void AddVertex(localIndex er, localIndex esr, globalIndex ei);

  /**
  * @brief Recover a vertex by using its indexes
  * @param[in] er Region index
  * @param[in] esr Surregion index
  * @param[in] ei Global vertex index
  * @return the vertex with given indexes, provided it exists
  */ 
  std::shared_ptr<GraphVertex> getVertexWithGlobalIndex(localIndex er, localIndex esr, globalIndex ei);
  
  /**
   * @brief Function called to construct the graphs
  */
  virtual void GenerateGraph() override;

  /**
  * @brief Function called to partition of the graph using the partition of given mesh
  * @param[in] meshLevel the mesh already partitioned to copy from
  */
  virtual void PartitionGraph(const MeshLevel & meshLevel) override;

  virtual void RemapFace(const MeshLevel & mesh);
 
  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr auto fileString = "file";
    static constexpr auto meshString = "mesh";

     } viewKeys;
  /// @endcond

  Path getFile() const { return m_file; }

  array1d<GraphEdge*> getBoundaryEdges() const { return m_boundaryEdges; }
  

private:
  Path m_file;
  string m_meshString;
  array1d<GraphEdge*> m_boundaryEdges;
  

};

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHFROMTEXT_HPP_ */
