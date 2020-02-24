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

  virtual const string getCatalogName() const override final
  { return FaceManager::CatalogName(); }


  ///@}
  ///
  ///
  ///
  ///
  FaceManager( string const &, Group * const parent );
  virtual ~FaceManager() override final;


  void BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elemManager );

  void computeGeometry( NodeManager const * const nodeManager );

  localIndex getMaxFaceNodes() const;

  void SortAllFaceNodes( NodeManager const * const nodeManager,
                         ElementRegionManager const * const elemManager);

  void SortFaceNodes( arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X,
                      R1Tensor const & elemCenter,
                      localIndex * const faceNodes,
                      localIndex const numFaceNodes );

  void SetDomainBoundaryObjects( NodeManager * const nodeManager );

  void SetIsExternal();

  virtual void ViewPackingExclusionList( SortedArray<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  void FixUpDownMaps( bool const clearIfUnmapped );

  virtual void enforceStateFieldConsistencyPostTopologyChange( std::set<localIndex> const & targetIndices ) override;

  void depopulateUpMaps( std::set<localIndex> const & receivedFaces,
                         ElementRegionManager const & elemRegionManager );

  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const  nodeManager,
                                                   std::vector< std::vector< globalIndex > >& faceToNodes ) override final;

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

    dataRepository::ViewKey nodeList              = { nodeListString };
    dataRepository::ViewKey edgeList              = { edgeListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };
  } viewKeys;

  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;

  constexpr int maxFacesPerNode() const { return 100; }

  array1d<real64> &       faceArea()       { return m_faceArea; }
  array1d<real64> const & faceArea() const { return m_faceArea; }

  array1d<R1Tensor> &       faceCenter()       { return m_faceCenter; }
  array1d<R1Tensor> const & faceCenter() const { return m_faceCenter; }

  array1d<R1Tensor> &       faceNormal()       { return m_faceNormal; }
  array1d<R1Tensor> const & faceNormal() const { return m_faceNormal; }


  NodeMapType & nodeList()                    { return m_nodeList; }
  NodeMapType const & nodeList() const        { return m_nodeList; }

  EdgeMapType       & edgeList()       { return m_edgeList; }
  EdgeMapType const & edgeList() const { return m_edgeList; }

  array2d<localIndex>       & elementRegionList()       { return m_toElements.m_toElementRegion; }
  array2d<localIndex> const & elementRegionList() const { return m_toElements.m_toElementRegion; }

  array2d<localIndex>       & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
  array2d<localIndex> const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }

  array2d<localIndex>       & elementList()       { return m_toElements.m_toElementIndex; }
  array2d<localIndex> const & elementList() const { return m_toElements.m_toElementIndex; }

  ElemMapType       & toElementRelation()       { return m_toElements; }
  ElemMapType const & toElementRelation() const { return m_toElements; }

private:

  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;


  NodeMapType m_nodeList;
  EdgeMapType m_edgeList;
  ElemMapType m_toElements;

  map< localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInToNodes;
  map< localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInToEdges;

  array1d< real64 > m_faceArea;
  array1d< R1Tensor > m_faceCenter;
  array1d< R1Tensor > m_faceNormal;

  constexpr static int MAX_FACE_NODES = 9;

  FaceManager() = delete;
  FaceManager( FaceManager const &) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
