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

  dataRepository::Group const * GetConstitutiveModels() const
  { return &m_constitutiveModels; }

  dataRepository::Group * GetConstitutiveModels()
  { return &m_constitutiveModels; }

  virtual string GetElementTypeString() const { return m_elementTypeString; }

  virtual void SetElementType( string const & elementType );

private:
  dataRepository::Group m_constitutiveModels;

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

#endif /* GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_ */
