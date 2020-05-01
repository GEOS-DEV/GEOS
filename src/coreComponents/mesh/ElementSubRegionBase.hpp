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

/**
 * @class ElementSubRegionBase
 * Abstract class for a collection of mesh elements that
 * will be derived and specialized for cell elements,
 * face (fracture) elements, embedded surface (fracture)
 * elements, well elements, etc
 */
class ElementSubRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ElementSubRegionBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Destructor.
   */
  ~ElementSubRegionBase();

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  /**
   * @brief Calculate the geometric quantities for each element in the subregion.
   * @param[in] nodeManager the nodeManager (for geometrical info and connectivity involving nodes)
   * @param[in] faceManager the faceManager (for geometrical info and connectivity involving faces)
   *
   */
  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) = 0;

  /**
   * @brief Link the connectivity maps of the subregion to the managers storing the mesh information.
   * @param[in] mesh the meshLevel object (single level only)
   *
   * In the derived classes, this function is used to passe a pointer to the nodeManager,
   * faceManager, and (if needed) edgeManager to, respectively, the node list, face list, and
   * edge list of the subregion.
   */
  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) = 0;

  /**
   * @brief Call ObjectManagerBase::FixUpDownMaps for the connectivity maps needed by
   *        the derived class (i.e., element-to-node map, element-to-face map, etc)
   * @param[in] clearIfUnmapped clearIfUnmapped
   */
  virtual void FixUpDownMaps( bool const clearIfUnmapped ) { GEOSX_UNUSED_VAR( clearIfUnmapped ); }

  ///@}

  /**
   * @name Accessors / Setters
   */
  ///@{

  /**
   * @brief Const accessor for the number of nodes per element.
   * @return number of nodes per element
   */
  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }

  /**
   * @brief Const accessor for the number of nodes of a specific element.
   * @return the number of nodes in the element
   */
  virtual localIndex numNodesPerElement( localIndex const ) const { return m_numNodesPerElement; }

  /**
   * @brief Set the number of nodes per element.
   * @param[in] numNodes the number of nodes per element
   */
  void setNumNodesPerElement( localIndex numNodes )
  {
    m_numNodesPerElement = numNodes;
  }

  /**
   * @brief Const accessor for the number of independent nodes per element.
   * @return numNodes the number of independent nodes per element
   *
   * Currently, the number of independent nodes per element is always
   * equal to the number of nodes per element, except for the case
   * of a triangular prism
   */
  localIndex const & numIndependentNodesPerElement() const { return m_numIndependentNodesPerElement; }

  /**
   * @brief Set the number of independent nodes per element.
   * @param[in] numNodes the number of independent nodes per element
   */
  void setNumIndependentNodesPerElement( localIndex const numNodes )
  {
    m_numIndependentNodesPerElement = numNodes;
  }

  /**
   * @brief Const accessor for the number of edges per element.
   * @return the number of edges per element
   */
  localIndex const & numEdgesPerElement() const { return m_numEdgesPerElement; }

  /**
   * @brief Set the number of edges per element.
   * @param[in] numEdges the number of edges per element
   */
  void setNumEdgesPerElement( localIndex const numEdges )
  {
    m_numEdgesPerElement = numEdges;
  }

  /**
   * @brief Const accessor for the number of faces per element.
   * @return number of faces per element
   */
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }

  /**
   * @brief Set the number of faces per element.
   * @param[in] numFaces the number of faces per element
   */
  void setNumFacesPerElement( localIndex const numFaces )
  { m_numFacesPerElement = numFaces; }

  /**
   * @brief Const accessor for the center of each element in this subregion.
   * @return an arrayView1d of const element centers
   */
  arrayView1d< R1Tensor const > const & getElementCenter() const
  { return m_elementCenter; }

  /**
   * @brief Accessor for the center of each element in this subregion.
   * @return an arrayView1d of element centers
   */
  arrayView1d< R1Tensor > const & getElementCenter()
  { return m_elementCenter; }

  /**
   * @brief Accessor for the volume of each element in this subregion.
   * @return an arrayView1d of const element volumes
   */
  arrayView1d< real64 const > const & getElementVolume() const
  { return m_elementVolume; }

  /**
   * @brief Const accessor for the group in which the constitutive models of this subregion are registered.
   * @return a pointer to the const group in which the constitutive models are registered
   */
  dataRepository::Group const * GetConstitutiveModels() const
  { return &m_constitutiveModels; }

  /**
   * @brief Accessor for the group in which the constitutive models of this subregion are registered.
   * @return a pointer to the group in which the constitutive models are registered
   */
  dataRepository::Group * GetConstitutiveModels()
  { return &m_constitutiveModels; }

  /**
   * @brief Accessor for the type of element in this subregion.
   * @return a string specifying the type of element in this subregion
   *
   * See class FiniteElementBase for possible element type.
   */
  virtual string GetElementTypeString() const { return m_elementTypeString; }

  /**
   * @brief Set the type of element in this subregion.
   * @param[in] elementType a string specifying the element type
   */
  virtual void SetElementType( string const & elementType );

  /**
   * @brief Accessor for the VTK ordering for this subregion.
   * @return the VTK node ordering
   */
  std::vector< int > getVTKNodeOrdering() const;

  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {

    /// The string key for the number of nodes per element in this subregion.
    static constexpr auto numNodesPerElementString = "numNodesPerElement";
    /// The string key for the element-to-node relation
    static constexpr auto nodeListString           = "nodeList";
    /// The string key for the number of edges per element in this subregion.
    static constexpr auto numEdgesPerElementString = "numEdgesPerElement";
    /// The string key for the element-to-edge relation
    static constexpr auto edgeListString           = "edgeList";
    /// The string key for the number of faces per element in this subregion.
    static constexpr auto numFacesPerElementString = "numFacesPerElement";
    /// The string key for the element-to-face relation
    static constexpr auto faceListString           = "faceList";
    /// The string key for the member level field for the element center.
    static constexpr auto elementCenterString      = "elementCenter";
    /// The string key for the member level field for the element volume.
    static constexpr auto elementVolumeString      = "elementVolume";
  };

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// The key for the group in which the constitutive models of this subregion are registered.
    static constexpr auto constitutiveModelsString = "ConstitutiveModels";
  };

private:
  /// The group in which the constitutive models of this subregion are registered
  dataRepository::Group m_constitutiveModels;

protected:
  /// The number of nodes per element in this subregion.
  localIndex m_numNodesPerElement;

  /// The number of independent nodes per element in this subregion.
  localIndex m_numIndependentNodesPerElement;

  /// The number of edges per element in this subregion.
  localIndex m_numEdgesPerElement;

  /// The number of faces per element in this subregion.
  localIndex m_numFacesPerElement;

  /// The member level field for the element center.
  array1d< R1Tensor > m_elementCenter;

  /// The member level field for the element volume.
  array1d< real64 > m_elementVolume;

  /// The type of element in this subregion.
  string m_elementTypeString;

  /// The type of element in this subregion.
  FiniteElementBase::ElementType m_elementType;
};


} /* namespace geosx */

#endif /* GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_ */
