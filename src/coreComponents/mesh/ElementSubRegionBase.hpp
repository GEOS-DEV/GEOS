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
 * @file ElementSubRegionBase.hpp
 */

#ifndef GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_
#define GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_

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

  ElementSubRegionBase( string const & name, dataRepository::Group * const parent );
  ~ElementSubRegionBase();

  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) = 0;

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) = 0;

  virtual void FixUpDownMaps( bool const GEOSX_UNUSED_PARAM( clearIfUnmapped ) ) {}

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


  /**
   * @return number of nodes per element
   */
  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }

  /**
   * Set the number of independent nodes per element.
   * @param numNodes The number of independent nodes per element.
   */
  void setNumIndependentNodesPerElement( localIndex const numNodes )
  {
    m_numIndependentNodesPerElement = numNodes;
  }

  localIndex const & numIndependentNodesPerElement() const { return m_numIndependentNodesPerElement; }

  /**
   * Sets the number of nodes per element
   * @param numNodes The number of nodes per element
   */
  void setNumNodesPerElement( localIndex numNodes )
  {
    m_numNodesPerElement = numNodes;
  }

  virtual localIndex numNodesPerElement( localIndex const ) const { return m_numNodesPerElement; }

  /**
   * @return number of edges per element
   */
  localIndex const & numEdgesPerElement() const { return m_numEdgesPerElement; }

  /**
   * Sets the number of edges per element
   * @param numEdges The number of edges per element
   */
  void setNumEdgesPerElement( localIndex const numEdges )
  {
    m_numEdgesPerElement = numEdges;
  }

  /**
   * @return number of faces per element
   */
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }

  /**
   * Sets the number of faces per element
   * @param numFaces The number of faces per element
   */
  void setNumFacesPerElement( localIndex const numFaces )
  { m_numFacesPerElement = numFaces; }

  arrayView1d< R1Tensor const > const & getElementCenter() const
  { return m_elementCenter; }

  arrayView1d< R1Tensor > const & getElementCenter()
  { return m_elementCenter; }

  arrayView1d< real64 const > const & getElementVolume() const
  { return m_elementVolume; }

  dataRepository::Group const * GetConstitutiveModels() const
  { return &m_constitutiveModels; }

  dataRepository::Group * GetConstitutiveModels()
  { return &m_constitutiveModels; }

  virtual string GetElementTypeString() const { return m_elementTypeString; }

  virtual void SetElementType( string const & elementType );

  std::vector< int > getVTKNodeOrdering() const;

private:
  dataRepository::Group m_constitutiveModels;

protected:
  /// The number of nodes per element in this cell block
  localIndex m_numNodesPerElement;

  /// The number of independent nodes per element in this cell block
  localIndex m_numIndependentNodesPerElement;

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
};


} /* namespace geosx */

#endif /* GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_ */
