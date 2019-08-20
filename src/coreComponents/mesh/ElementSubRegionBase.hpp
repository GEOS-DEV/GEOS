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
 * @file CellBase.hpp
 */

#ifndef SRC_CORECOMPONENTS_MESH_CELLBASE_HPP_
#define SRC_CORECOMPONENTS_MESH_CELLBASE_HPP_

#include "managers/ObjectManagerBase.hpp"
#include "finiteElement/ElementLibrary/FiniteElementBase.h"
namespace geosx
{

class NodeManager;
class FaceManager;
class MeshLevel;
class DomainPartition;

class ElementSubRegionBase : public ObjectManagerBase
{
public:
  ElementSubRegionBase( string const & name, dataRepository::ManagedGroup * const parent );
  ~ElementSubRegionBase();

  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                                   const NodeManager& nodeManager,
                                                   const bool useReferencePos = true) const = 0;

  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) = 0;

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) = 0;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) {}

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
  };

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto constitutiveModelsString = "ConstitutiveModels";
  };


  virtual arraySlice1dRval<localIndex const> nodeList( localIndex const k ) const = 0;
  virtual arraySlice1dRval<localIndex> nodeList( localIndex const k ) = 0;


  /**
   * @return number of nodes per element
   */
  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }

  /**
   * @return number of nodes per element
   */
  localIndex       & numNodesPerElement()       { return m_numNodesPerElement; }

  virtual localIndex numNodesPerElement( localIndex const ) const { return m_numNodesPerElement; }

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

  array1d< R1Tensor > const & getElementCenter() const
  {
    return m_elementCenter;
  }

  array1d< real64 > const & getElementVolume() const
  {
    return m_elementVolume;
  }

  dataRepository::ManagedGroup const * GetConstitutiveModels() const
  { return &m_constitutiveModels; }

  dataRepository::ManagedGroup * GetConstitutiveModels()
  { return &m_constitutiveModels; }

  virtual string GetElementTypeString() const { return m_elementTypeString; }

  FiniteElementBase::ElementType GetElementType() const { return m_elementType; }

  virtual void SetElementType( string const & elementType );

private:
  dataRepository::ManagedGroup m_constitutiveModels;

protected:
  /// The number of nodes per element in this cell block
  localIndex m_numNodesPerElement;

  /// The number of edges per element in this cell block
  localIndex m_numEdgesPerElement;

  /// The number of faces per element in this cell block
  localIndex m_numFacesPerElement;

  /// The member level field for the element center
  array1d< R1Tensor > m_elementCenter;

  /// The member level field for the element volume
  array1d< real64 > m_elementVolume;

  string m_elementTypeString;

  FiniteElementBase::ElementType m_elementType;
//  template< LAMBDA lambda >
//  void numNodesPerElemSwitchyard() const;
};


} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_MESH_CELLBASE_HPP_ */
