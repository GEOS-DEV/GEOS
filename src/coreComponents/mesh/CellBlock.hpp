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
 * @file CellBlock.hpp
 */

#ifndef ELEMENTOBJECTT_H_
#define ELEMENTOBJECTT_H_

#include "ElementSubRegionBase.hpp"
#include "FaceManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"


class StableTimeStep;

namespace geosx
{

/**
 * Class to manage the data stored at the element level.
 */
class CellBlock : public ElementSubRegionBase
{
public:

  using NodeMapType=FixedOneToManyRelation;
  using EdgeMapType=FixedOneToManyRelation;
  using FaceMapType=FixedOneToManyRelation;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   *
   * @return the name of this type in the catalog
   */
  static const string CatalogName()
  { return "CellBlock"; }

  /**
   *
   * @return the name of this type in the catalog
   */
  virtual const string getCatalogName() const override final
  { return CellBlock::CatalogName(); }


  ///@}


  /// deleted default constructor
  CellBlock() = delete;

  /**
   * @brief constructor
   * @param name the name of the object in the data repository
   * @param parent the parent object of this object in the data repository
   */
  CellBlock( string const & name, ManagedGroup * const parent );

  /**
   * @brief copy constructor
   * @param init the source to copy
   */
  CellBlock(const CellBlock& init);


  virtual ~CellBlock() override;

  /**
   * @brief function to return the localIndices of the nodes in a face of the element
   * @param elementIndex The localIndex of the target element
   * @param localFaceIndex the element local localIndex of the target face (this will be 0-numFacesInElement
   * @param nodeIndicies the node indices of the face
   */
  void GetFaceNodes( const localIndex elementIndex,
                     const localIndex localFaceIndex,
                     localIndex_array& nodeIndicies) const;

  /**
   * @brief function to return element center. this should be depricated.
   * @param k
   * @param nodeManager
   * @param useReferencePos
   * @return
   */
  R1Tensor const & calculateElementCenter( localIndex k,
                                           const NodeManager& nodeManager,
                                           const bool useReferencePos = true) const override;

  virtual void CalculateCellVolumes( array1d<localIndex> const & indices,
                                     array1d<R1Tensor> const & X ) override
  {
    ElementSubRegionBase::CalculateCellVolumes<CellBlock>( *this,
                                               indices,
                                               X );
  }

  inline void CalculateCellVolumesKernel( localIndex const k,
                                          array1d<R1Tensor> const & X )
  {
    R1Tensor & center = m_elementCenter[k];
    center = 0.0;

    R1Tensor Xlocal[10];

    for (localIndex a = 0; a < m_numNodesPerElement; ++a)
    {
      Xlocal[a] = X[m_toNodesRelation[k][a]];
      center += Xlocal[a];
    }
    center /= m_numNodesPerElement;

    if( m_numNodesPerElement == 8 )
    {
      m_elementVolume[k] = computationalGeometry::HexVolume(Xlocal);
    }
    else if( m_numNodesPerElement == 4)
    {
      m_elementVolume[k] = computationalGeometry::TetVolume(Xlocal);
    }
    else if( m_numNodesPerElement == 6)
    {
      m_elementVolume[k] = computationalGeometry::WedgeVolume(Xlocal);
    }
    else if ( m_numNodesPerElement == 5)
    {
      m_elementVolume[k] = computationalGeometry::PyramidVolume(Xlocal);
    }
    else
    {
        GEOS_ERROR("GEOX does not support cells with " << m_numNodesPerElement << " nodes");
    }
  }


  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;


  virtual arraySlice1dRval<localIndex const> nodeList( localIndex const k ) const override
  {
    return m_toNodesRelation[k];
  }

  virtual arraySlice1dRval<localIndex> nodeList( localIndex const k ) override
  {
    return m_toNodesRelation[k];
  }


  /**
   * @return the element to node map
   */
  FixedOneToManyRelation & nodeList()                    { return m_toNodesRelation; }

  /**
   * @return the element to node map
   */
  FixedOneToManyRelation const & nodeList() const        { return m_toNodesRelation; }

  /**
   * @return the element to node map
   */
  localIndex & nodeList( localIndex const k, localIndex a ) { return m_toNodesRelation[k][a]; }

  /**
   * @return the element to node map
   */
  localIndex const & nodeList( localIndex const k, localIndex a ) const { return m_toNodesRelation[k][a]; }

  /**
   * @return the element to edge map
   */
  FixedOneToManyRelation       & edgeList()       { return m_toEdgesRelation; }

  /**
   * @return the element to edge map
   */
  FixedOneToManyRelation const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @return the element to face map
   */
  FixedOneToManyRelation       & faceList()       { return m_toFacesRelation; }

  /**
   * @return the element to face map
   */
  FixedOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  string GetElementType() const { return m_elementType; }

  void SetElementType( string const & elementType);

protected:


  /// The elements to nodes relation
  NodeMapType  m_toNodesRelation;

  /// The elements to edges relation
  EdgeMapType  m_toEdgesRelation;

  /// The elements to faces relation
  FaceMapType  m_toFacesRelation;


//  CellBlock& operator=(const CellBlock& rhs);
  string m_elementType;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
