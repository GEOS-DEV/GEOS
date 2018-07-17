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

  void  AddNewFace( localIndex const & kReg,
                    localIndex const & kSubReg,
                    localIndex const & ke,
                    localIndex const & kelf,
                    localIndex & numFaces,
                    array<localIndex_array>& facesByLowestNode,
                    localIndex_array& tempNodeList,
                    array<localIndex_array>& tempFaceToNodeMap,
                    CellBlockSubRegion & elementRegion );



  void SortAllFaceNodes( NodeManager const & nodeManager,
                         ElementRegionManager const & elemManager);

  void SortFaceNodes( NodeManager const & nodeManager,
                      R1Tensor const & elementCenter,
                      const localIndex faceIndex );

  void SetDomainBoundaryObjects( NodeManager * const nodeManager );

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( localIndex_array const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                              localIndex_array const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                localIndex_array const & packList ) override;


  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const & nodeManager,
                                                   array<globalIndex_array>& faceToNodes ) override final;
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nodeListString              = "nodeList";
    static constexpr auto edgeListString              = "edgeList";
    static constexpr auto elementRegionListString     = "elemRegionList";
    static constexpr auto elementSubRegionListString  = "elemSubRegionList";
    static constexpr auto elementListString           = "elemList";

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

  Array2dT<localIndex>       & elementRegionList()       { return m_toElements.m_toElementRegion; }
  Array2dT<localIndex> const & elementRegionList() const { return m_toElements.m_toElementRegion; }

  Array2dT<localIndex>       & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
  Array2dT<localIndex> const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }

  Array2dT<localIndex>       & elementList()       { return m_toElements.m_toElementIndex; }
  Array2dT<localIndex> const & elementList() const { return m_toElements.m_toElementIndex; }


private:

  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const;


  OrderedVariableOneToManyRelation m_nodeList;
  OrderedVariableOneToManyRelation m_edgeList;
  FixedToManyElementRelation m_toElements;

  FaceManager() = delete;
  FaceManager( FaceManager const &) = delete;
  FaceManager( FaceManager && ) = delete;
};

}
#endif /* FACEMANAGERT_H_ */
