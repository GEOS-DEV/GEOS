/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_GENERATORS_MESHBLOCKMAPPINGS_HPP
#define GEOS_GENERATORS_MESHBLOCKMAPPINGS_HPP

#include "CellBlockManagerABC.hpp"

namespace geos
{

class MeshBlockMappings: public CellBlockManagerABC
{
public:
  MeshBlockMappings( string const & name, Group * const parent );

  virtual Group & getCellBlocks() override;

  virtual Group & getFaceBlocks() override;

  virtual LineBlockABC const & getLineBlock( string name ) const override;

  virtual Group const & getCellBlocks() const override;

  virtual Group const & getFaceBlocks() const override;

  virtual localIndex numNodes() const override;

  virtual localIndex numEdges() const override;

  virtual localIndex numFaces() const override;

  virtual array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodePositions() const override;

  virtual ArrayOfArrays< localIndex > getNodeToEdges() const override;

  virtual ArrayOfArrays< localIndex > getNodeToFaces() const override;

  virtual ToCellRelation< ArrayOfArrays< localIndex>> getNodeToElements() const override;

  virtual array2d< localIndex > getEdgeToNodes() const override;

  virtual ArrayOfArrays< localIndex > getEdgeToFaces() const override;

  virtual ArrayOfArrays< localIndex > getFaceToNodes() const override;

  virtual ArrayOfArrays< localIndex > getFaceToEdges() const override;

  virtual ToCellRelation< array2d< localIndex>> getFaceToElements() const override;

  virtual array1d< globalIndex > getNodeLocalToGlobal() const override;

  virtual std::map< string, SortedArray< localIndex>> const & getNodeSets() const override;

  virtual real64 getGlobalLength() const override;

  virtual void generateHighOrderMaps( localIndex const order,
                                      globalIndex const maxVertexGlobalID,
                                      globalIndex const maxEdgeGlobalID,
                                      globalIndex const maxFaceGlobalID,
                                      arrayView1d< globalIndex const > const edgeLocalToGlobal,
                                      arrayView1d< globalIndex const > const faceLocalToGlobal ) override;

private:

  struct viewKeyStruct
  {
    /// Cell blocks key
    static constexpr char const * cellBlocks()
    { return "cellBlocks"; }

    /// Face blocks key
    static constexpr char const * faceBlocks()
    { return "faceBlocks"; }

    /// Line blocks key
    static constexpr char const * lineBlocks()
    { return "lineBlocks"; }
  };

  std::map< string, SortedArray< localIndex > > m_nodeSets;

  real64 m_globalLength;

  localIndex m_numNodes;
  localIndex m_numFaces;
  localIndex m_numEdges;
};

} // geos

#endif //GEOS_GENERATORS_MESHBLOCKMAPPINGS_HPP
