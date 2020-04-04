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
 * @file MeshLevel.hpp
 */

#ifndef GEOSX_MESH_MESHLEVEL_HPP_
#define GEOSX_MESH_MESHLEVEL_HPP_

#include "NodeManager.hpp"
#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
class ElementRegionManager;

class MeshLevel : public dataRepository::Group
{
public:
  MeshLevel( string const & name,
             Group * const parent );
  virtual ~MeshLevel() override;

  void GenerateAdjacencyLists( arrayView1d< localIndex const > const & seedNodeList,
                               localIndex_array & nodeAdjacencyList,
                               localIndex_array & edgeAdjacencyList,
                               localIndex_array & faceAdjacencyList,
                               ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > & elementAdjacencyList,
                               integer const depth );


  virtual void InitializePostInitialConditions_PostSubGroups( Group * const ) override;


  struct viewStructKeys
  {
    dataRepository::ViewKey meshLevel                = { "meshLevel" };
  } viewKeys;

  struct groupStructKeys
  {
    dataRepository::GroupKey vertexManager  = { "vertexManager" };
    dataRepository::GroupKey cellManager    = { "cellManager" };

    static constexpr auto nodeManagerString = "nodeManager";
    static constexpr auto edgeManagerString = "edgeManager";
    static constexpr auto faceManagerString = "faceManager";

    // This key is defined in problem manager:
    static constexpr auto elemManagerString = "ElementRegions";

    dataRepository::GroupKey nodeManager = {nodeManagerString};
    dataRepository::GroupKey edgeManager = {edgeManagerString};
    dataRepository::GroupKey faceManager = {faceManagerString};
    dataRepository::GroupKey elemManager = {elemManagerString};
  } groupKeys;

  NodeManager const * getNodeManager() const { return &m_nodeManager; }
  NodeManager * getNodeManager()             { return &m_nodeManager; }

  EdgeManager const * getEdgeManager() const { return &m_edgeManager; }
  EdgeManager * getEdgeManager()             { return &m_edgeManager; }

  FaceManager const * getFaceManager() const { return &m_faceManager; }
  FaceManager * getFaceManager()             { return &m_faceManager; }

  ElementRegionManager const * getElemManager() const { return &m_elementManager; }
  ElementRegionManager * getElemManager()             { return &m_elementManager; }

private:

  NodeManager m_nodeManager;
  EdgeManager m_edgeManager;
  FaceManager m_faceManager;
  ElementRegionManager m_elementManager;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHLEVEL_HPP_ */
