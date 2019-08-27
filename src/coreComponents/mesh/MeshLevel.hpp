/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

/**
 * @file MeshLevel.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_

#include "NodeManager.hpp"
#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
class ElementRegionManager;

class MeshLevel : public dataRepository::ManagedGroup
{
public:
  MeshLevel( string const & name,
             ManagedGroup * const parent );
  virtual ~MeshLevel() override;

  void GenerateAdjacencyLists( localIndex_array & seedNodeList,
                               localIndex_array & nodeAdjacencyList,
                               localIndex_array & edgeAdjacencyList,
                               localIndex_array & faceAdjacencyList,
                               ElementRegionManager::ElementViewAccessor<ReferenceWrapper<localIndex_array>> & elementAdjacencyList,
                               integer const depth );


  virtual void InitializePostInitialConditions_PostSubGroups( ManagedGroup * const ) override;


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

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_MESHLEVEL_HPP_ */
