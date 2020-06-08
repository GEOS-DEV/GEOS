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
 * @file FaceManager.hpp
 */

#ifndef GEOSX_MESH_FACEMANAGER_HPP_
#define GEOSX_MESH_FACEMANAGER_HPP_

#include "ToElementRelation.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class NodeManager;
class ElementRegionManager;
class CellElementSubRegion;

class FaceManager : public ObjectManagerBase
{
public:

  using NodeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  using ElemMapType = FixedToManyElementRelation;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static const string CatalogName()
  { return "FaceManager"; }

  virtual const string getCatalogName() const override
  { return FaceManager::CatalogName(); }

  static localIndex nodeMapExtraSpacePerFace()
  { return 4; }

  static localIndex edgeMapExtraSpacePerFace()
  { return 4; }

  ///@}
  ///
  ///
  ///
  ///
  FaceManager( string const &, Group * const parent );
  virtual ~FaceManager() override;

  virtual void resize( localIndex const newsize ) override;

  void BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elemManager );

  void computeGeometry( NodeManager const * const nodeManager );

  localIndex getMaxFaceNodes() const;

  void SortAllFaceNodes( NodeManager const * const nodeManager,
                         ElementRegionManager const * const elemManager );

  void SortFaceNodes( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                      arraySlice1d< real64 const > const elementCenter,
                      localIndex * const faceNodes,
                      localIndex const numFaceNodes );

  void SetDomainBoundaryObjects( NodeManager * const nodeManager );

  void SetIsExternal();

  virtual void ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  void FixUpDownMaps( bool const clearIfUnmapped );

  void compressRelationMaps();

  virtual void enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices ) override;

  void depopulateUpMaps( std::set< localIndex > const & receivedFaces,
                         ElementRegionManager const & elemRegionManager );

  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                   std::vector< std::vector< globalIndex > > & faceToNodes ) override;

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nodeListString              = "nodeList";
    static constexpr auto edgeListString              = "edgeList";
    static constexpr auto elementRegionListString     = "elemRegionList";
    static constexpr auto elementSubRegionListString  = "elemSubRegionList";
    static constexpr auto elementListString           = "elemList";
    constexpr static auto faceAreaString = "faceArea";
    constexpr static auto faceCenterString = "faceCenter";
    constexpr static auto faceNormalString = "faceNormal";
    constexpr static auto faceRotationMatrixString = "faceRotationMatrix";

    dataRepository::ViewKey nodeList              = { nodeListString };
    dataRepository::ViewKey edgeList              = { edgeListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };
  } viewKeys;

  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;

  constexpr int maxFacesPerNode() const { return 100; }

  arrayView1d< real64 const > const & faceArea() const { return m_faceArea; }

  arrayView2d< real64 const > const & faceCenter() const { return m_faceCenter; }

  arrayView2d< real64 > const & faceNormal() { return m_faceNormal; }
  arrayView2d< real64 const > const & faceNormal() const { return m_faceNormal; }

  arrayView3d< real64 const > const & faceRotationMatrix() const { return m_faceRotationMatrix; }

  NodeMapType & nodeList()                    { return m_nodeList; }
  NodeMapType const & nodeList() const { return m_nodeList; }

  EdgeMapType & edgeList()       { return m_edgeList; }
  EdgeMapType const & edgeList() const { return m_edgeList; }

  arrayView2d< localIndex > const & elementRegionList() { return m_toElements.m_toElementRegion; }
  arrayView2d< localIndex const > const & elementRegionList() const { return m_toElements.m_toElementRegion; }

  arrayView2d< localIndex > const & elementSubRegionList() { return m_toElements.m_toElementSubRegion; }
  arrayView2d< localIndex const > const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }

  arrayView2d< localIndex > const & elementList() { return m_toElements.m_toElementIndex; }
  arrayView2d< localIndex const > const & elementList() const { return m_toElements.m_toElementIndex; }

  ElemMapType & toElementRelation()       { return m_toElements; }
  ElemMapType const & toElementRelation() const { return m_toElements; }

private:

  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;


  NodeMapType m_nodeList;
  EdgeMapType m_edgeList;
  ElemMapType m_toElements;

  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  array1d< real64 > m_faceArea;
  array2d< real64 > m_faceCenter;
  array2d< real64 > m_faceNormal;
  array3d< real64 > m_faceRotationMatrix;

  constexpr static int MAX_FACE_NODES = 9;

  FaceManager() = delete;
  FaceManager( FaceManager const & ) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
