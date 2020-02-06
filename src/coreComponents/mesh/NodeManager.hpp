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
 * @file NodeManager.hpp
 */


#ifndef GEOSX_MESH_NODEMANAGER_HPP_
#define GEOSX_MESH_NODEMANAGER_HPP_

#include "managers/ObjectManagerBase.hpp"
#include <string.h>
#include "CellBlockManager.hpp"
#include "ToElementRelation.hpp"


// *********************************************************************************************************************
// *********************************************************************************************************************
class SiloFile;

namespace geosx
{

class CellBlock;
class FaceManager;
class EdgeManager;
class ElementRegionManager;


/**
 * @class NodeManager
 * @brief The NodeManager class provides an interface to ObjectManagerBase in order to manage node data.
 *
 * The NodeManagerT class manages the node data using the
 * ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry
 * corresponds to a node.
 */
class NodeManager : public ObjectManagerBase
{
public:

  //START_SPHINX_INCLUDE_01
  using EdgeMapType = InterObjectRelation< ArrayOfSets< localIndex > >;
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;
  using ElemMapType = OrderedVariableToManyElementRelation;
  //END_SPHINX_INCLUDE_01

  /**
   * @brief main constructor for NodeManager Objects
   * @param name the name of this instantiation of NodeManager in the repository
   * @param parent the parent group of this instantiation of NodeManager
   */
  NodeManager( std::string const & name,
               dataRepository::Group * const parent );

  /**
   *  @brief default destructor
   */
  ~NodeManager() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "NodeManager"; }

  /**
   * @brief virtual access to CatalogName()
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  const string getCatalogName() const override final
  { return NodeManager::CatalogName(); }

  void SetEdgeMaps( EdgeManager const * const edgeManager );

  void SetFaceMaps( FaceManager const * const faceManager );

  void SetElementMaps( ElementRegionManager const * const elementRegionManager );

//  void Initialize();

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  void FixUpDownMaps( bool const clearIfUnmapped );

  void depopulateUpMaps( std::set<localIndex> const & receivedNodes,
                         array2d< localIndex > const & edgesToNodes,
                         ArrayOfArraysView< localIndex const > const & facesToNodes,
                         ElementRegionManager const & elemRegionManager );

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto referencePositionString     = "ReferencePosition";
    static constexpr auto totalDisplacementString     = "TotalDisplacement";
    static constexpr auto edgeListString              = "edgeList";
    static constexpr auto faceListString              = "faceList";
    static constexpr auto elementRegionListString     = "elemRegionList";
    static constexpr auto elementSubRegionListString  = "elemSubRegionList";
    static constexpr auto elementListString           = "elemList";

    dataRepository::ViewKey referencePosition = { referencePositionString };
    dataRepository::ViewKey totalDisplacement = { totalDisplacementString };
    dataRepository::ViewKey edgeList           = { edgeListString };
    dataRepository::ViewKey faceList           = { faceListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };

  } viewKeys;


  /**
   * @struct
   */
  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;


  /**
   * \defgroup accessors for NodeManager fixed data
   * @{
   */



  /**
   * @brief const accessor to the node->edge relation
   * @return const reference to relation
   */
  EdgeMapType const & edgeList() const
  { return m_toEdgesRelation; }

  /**
   * @brief accessor to the node->edge relation
   * @return reference to relation
   */
  EdgeMapType & edgeList()
  { return m_toEdgesRelation; }

  FaceMapType       & faceList()       { return m_toFacesRelation; }
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  OrderedVariableToManyElementRelation & toElementRelation() {return m_toElements;}
  OrderedVariableToManyElementRelation const & toElementRelation() const {return m_toElements;}

  ArrayOfArrays< localIndex > & elementRegionList()       { return m_toElements.m_toElementRegion; }
  ArrayOfArraysView< localIndex const > const & elementRegionList() const
  { return m_toElements.m_toElementRegion.toViewCC(); }

  ArrayOfArrays< localIndex > & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
  ArrayOfArraysView< localIndex const > const & elementSubRegionList() const
  { return m_toElements.m_toElementSubRegion.toViewCC(); }

  ArrayOfArrays< localIndex > & elementList()       { return m_toElements.m_toElementIndex; }
  ArrayOfArraysView< localIndex const > const & elementList() const
  { return m_toElements.m_toElementIndex.toViewCC(); }

  /**
   * @brief const accessor to the reference position array
   * @return const reference to reference position
   */
  array1d<R1Tensor> const & referencePosition() const
  { return m_referencePosition; }

  /**
   * @brief accessor to the reference position array
   * @return reference to reference position
   */
  array1d<R1Tensor> & referencePosition()
  { return m_referencePosition; }

protected:

private:
  /**
   * @brief function to pack the upward and downward pointing maps.
   * @tparam DOPACK template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                packed and the function returns the size of the packing that would have occured if set to TRUE.
   * @param buffer the buffer to pack data into
   * @param packList the indices of nodes that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;

  /// reference position of the nodes
  array1d<R1Tensor> m_referencePosition;

  /// nodeToEdge relation
  EdgeMapType m_toEdgesRelation;

  /// nodeToFace relation
  FaceMapType m_toFacesRelation;

  /// nodeToElement relation
  ElemMapType m_toElements;

  map< localIndex, set<globalIndex> > m_unmappedGlobalIndicesInToEdges;
  map< localIndex, set<globalIndex> > m_unmappedGlobalIndicesInToFaces;
  map< localIndex, array1d< array1d< set<globalIndex> > > > m_unmappedGlobalIndicesInToElems;



  /// deleted constructor
  NodeManager() = delete;

  /// deleted copy constructor
  NodeManager( const NodeManager& init ) = delete;

  /// deleted assignement operator
  NodeManager& operator=( const NodeManager&) = delete;

};
}


#endif // MESH_NODEMANAGER_HPP_
