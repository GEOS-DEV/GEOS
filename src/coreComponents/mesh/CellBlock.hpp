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
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTOBJECTT_H_
#define ELEMENTOBJECTT_H_

#include "managers/ObjectManagerBase.hpp"
#include "FaceManager.hpp"


class StableTimeStep;

namespace geosx
{

/**
 * Class to manage the data stored at the element level.
 */
class CellBlock : public ObjectManagerBase
{
public:

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
  R1Tensor GetElementCenter(localIndex k, const NodeManager& nodeManager, const bool useReferencePos = true) const;

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto numNodesPerElementString     = "numNodesPerElement";
    static constexpr auto nodeListString               = "nodeList";
    static constexpr auto numEdgesPerElementString     = "numEdgesPerElement";
    static constexpr auto edgeListString               = "edgeList";
    static constexpr auto numFacesPerElementString     = "numFacesPerElement";
    static constexpr auto faceListString               = "faceList";
    static constexpr auto elementCenterString          = "elementCenter";
    static constexpr auto elementVolumeString          = "elementVolume";

    dataRepository::ViewKey numNodesPerElement = { numNodesPerElementString };
    dataRepository::ViewKey nodeList           = { nodeListString };
    dataRepository::ViewKey numEdgesPerElement = { numEdgesPerElementString };
    dataRepository::ViewKey edgeList           = { edgeListString };
    dataRepository::ViewKey numFacesPerElement = { numFacesPerElementString };
    dataRepository::ViewKey faceList           = { faceListString };
    dataRepository::ViewKey elementCenter      = { elementCenterString };
  } m_CellBlockViewKeys;


  /**
   *
   * @return reference to the viewKeyStruct member
   */
  virtual viewKeyStruct & viewKeys() override
  { return m_CellBlockViewKeys; }

  /**
   *
   * @return reference to const pointing to the viewKeyStruct member
   */
  virtual viewKeyStruct const & viewKeys() const override
  { return m_CellBlockViewKeys; }




  /**
   * @return number of nodes per element
   */
  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }

  /**
   * @return number of nodes per element
   */
  localIndex       & numNodesPerElement()       { return m_numNodesPerElement; }

  /**
   * @return number of edges per element
   */
  localIndex const & numEdgesPerElement() const { return m_numEdgesPerElement; }

  /**
   * @return number of edges per element
   */
  localIndex       & numEdgesPerElement()       { return m_numEdgesPerElement; }

  /**
   * @return number of faces per element
   */
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }

  /**
   * @return number of faces per element
   */
  localIndex       & numFacesPerElement()       { return m_numFacesPerElement; }

  /**
   * @return the element to node map
   */
  FixedOneToManyRelation & nodeList()                    { return m_toNodesRelation; }

  /**
   * @return the element to node map
   */
  FixedOneToManyRelation const & nodeList() const        { return m_toNodesRelation; }

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

private:

  /// The number of nodes per element in this cell block
  localIndex m_numNodesPerElement;

  /// The number of edges per element in this cell block
  localIndex m_numEdgesPerElement;

  /// The number of faces per element in this cell block
  localIndex m_numFacesPerElement;

  /// The elements to nodes relation
  FixedOneToManyRelation  m_toNodesRelation;

  /// The elements to edges relation
  FixedOneToManyRelation  m_toEdgesRelation;

  /// The elements to faces relation
  FixedOneToManyRelation  m_toFacesRelation;

  /// The member level field for the element center
  array1d< R1Tensor > m_elementCenter;

  /// The member level field for the element volume
  array1d< real64 > m_elementVolume;

//  CellBlock& operator=(const CellBlock& rhs);
  string m_elementType;

};



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
