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
 * @file NodeManager.hpp
 */


#ifndef MESH_NODEMANAGER_HPP_
#define MESH_NODEMANAGER_HPP_

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

  /**
   * @brief main constructor for NodeManager Objects
   * @param name the name of this instantiation of NodeManager in the repository
   * @param parent the parent group of this instantiation of NodeManager
   */
  NodeManager( std::string const & name,
               dataRepository::ManagedGroup * const parent );

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


  void FillDocumentationNode() override final;


  void SetEdgeMaps( EdgeManager const * const edgeManager );

  void SetFaceMaps( FaceManager const * const faceManager );

  void SetElementMaps( ElementRegionManager const * const elementRegionManager );

//  void Initialize();

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex> const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                              arrayView1d<localIndex> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                arrayView1d<localIndex> const & packList ) override;

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
  UnorderedVariableOneToManyRelation const & edgeList() const
  { return m_toEdgesRelation; }

  /**
   * @brief accessor to the node->edge relation
   * @return reference to relation
   */
  UnorderedVariableOneToManyRelation & edgeList()
  { return m_toEdgesRelation; }

  UnorderedVariableOneToManyRelation       & faceList()       { return m_toFacesRelation; }
  UnorderedVariableOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  OrderedVariableToManyElementRelation & toElementRelation() {return m_toElements;}
  OrderedVariableToManyElementRelation const & toElementRelation() const {return m_toElements;}

  array1d<localIndex_array>       & elementRegionList()       { return m_toElements.m_toElementRegion; }
  array1d<localIndex_array> const & elementRegionList() const { return m_toElements.m_toElementRegion; }

  array1d<localIndex_array>       & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
  array1d<localIndex_array> const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }

  array1d<localIndex_array>        & elementList()       { return m_toElements.m_toElementIndex; }
  array1d<localIndex_array>  const & elementList() const { return m_toElements.m_toElementIndex; }
//  array1d<set<localIndex>>       & elementRegionList()       { return m_toElements.m_toElementRegion; }
//  array1d<set<localIndex>> const & elementRegionList() const { return m_toElements.m_toElementRegion; }
//
//  array1d<set<localIndex>>       & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }
//  array1d<set<localIndex>> const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }
//
//  array1d<set<localIndex>>        & elementList()       { return m_toElements.m_toElementIndex; }
//  array1d<set<localIndex>>  const & elementList() const { return m_toElements.m_toElementIndex; }


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
                                    arrayView1d<localIndex> const & packList ) const;

   /// reference position of the nodes
  array1d<R1Tensor> m_referencePosition;

  /// nodeToEdge relation
  UnorderedVariableOneToManyRelation m_toEdgesRelation;

  /// nodeToFace relation
  UnorderedVariableOneToManyRelation m_toFacesRelation;

  /// nodeToElement relation
  OrderedVariableToManyElementRelation m_toElements;

  /// deleted constructor
  NodeManager() = delete;

  /// deleted copy constructor
  NodeManager( const NodeManager& init ) = delete;

  /// deleted assignement operator
  NodeManager& operator=( const NodeManager&) = delete;

};
}


#endif // MESH_NODEMANAGER_HPP_
