/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file FaceManager.h
 * @author settgast1
 */

#ifndef FACEMANAGER_H_
#define FACEMANAGER_H_

#include "ToElementRelation.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class NodeManager;
class ElementRegionManager;
class CellBlockSubRegion;

class FaceManager : public ObjectManagerBase
{
public:

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
  FaceManager( string const &, ManagedGroup * const parent );
  virtual ~FaceManager() override final;

//  void Initialize(  ){}

  virtual void FillDocumentationNode() override final;


  void BuildFaces( NodeManager * const nodeManager, ElementRegionManager * const elemManager );

  localIndex getMaxFaceNodes() const;

  void SortAllFaceNodes( NodeManager const * const nodeManager,
                         ElementRegionManager const * const elemManager);

  void SortFaceNodes( array1d<R1Tensor> const & X,
                      R1Tensor const & elemCenter,
                      arrayView1d<localIndex> & faceNodes,
                      localIndex const numFaceNodes );

  void SetDomainBoundaryObjects( NodeManager * const nodeManager );

  void SetIsExternal();

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex> const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                              arrayView1d<localIndex> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                arrayView1d<localIndex> const & packList ) override;


  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const  nodeManager,
                                                   array1d<globalIndex_array>& faceToNodes ) override final;
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

  OrderedVariableOneToManyRelation & nodeList()                    { return m_nodeList; }
  OrderedVariableOneToManyRelation const & nodeList() const        { return m_nodeList; }

  OrderedVariableOneToManyRelation       & edgeList()       { return m_edgeList; }
  OrderedVariableOneToManyRelation const & edgeList() const { return m_edgeList; }

  array2d<localIndex>       & elementRegionList()       { return m_toElements.m_toElementRegion; }
  array2d<localIndex> const & elementRegionList() const { return m_toElements.m_toElementRegion; }

  array2d<localIndex>       & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
  array2d<localIndex> const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }

  array2d<localIndex>       & elementList()       { return m_toElements.m_toElementIndex; }
  array2d<localIndex> const & elementList() const { return m_toElements.m_toElementIndex; }


private:

  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                             arrayView1d<localIndex> const & packList ) const;


  OrderedVariableOneToManyRelation m_nodeList;
  OrderedVariableOneToManyRelation m_edgeList;
  FixedToManyElementRelation m_toElements;

  array1d< R1Tensor > m_faceCenter;

  constexpr static int MAX_FACE_NODES = 9;

  FaceManager() = delete;
  FaceManager( FaceManager const &) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
