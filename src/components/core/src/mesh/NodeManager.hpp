// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
/**
 * @file NodeManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */


#ifndef NODEMANAGERT_H_
#define NODEMANAGERT_H_

#include "managers/ObjectManagerBase.hpp"
#include <string.h>
#include "CellBlockManager.hpp"


// *********************************************************************************************************************
// *********************************************************************************************************************
class SiloFile;

namespace geosx
{

class CellBlock;
class FaceManager;
class EdgeManager;
class ElementRegionManager;

namespace dataRepository
{
namespace keys
{
std::string const nodeManager    = "NodeManager";
std::string const elementRegionMap("elementRegionMap");
std::string const elementSubRegionMap("elementSubRegionMap");
std::string const elementMap("elementMap");

}
}

using namespace dataRepository;
/**
 * @author Randolph Settgast
 *
 * The NodeManagerT class manages the node data using the
 * ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry
 * corresponds to a node.
 */
class NodeManager : public ObjectManagerBase
{
public:


  /// default constructor
  NodeManager( std::string const & name,
               ManagedGroup * const parent );



  /// default destructor
  ~NodeManager() override;

  static string CatalogName() { return dataRepository::keys::nodeManager; }
  const string getCatalogName() const override final
  { return NodeManager::CatalogName(); }


  void FillDocumentationNode() override final;


  void SetElementMaps( ElementRegionManager const * const elementRegionManager );

//  void Initialize();

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( localIndex_array const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                              localIndex_array const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                localIndex_array const & packList ) override;

public:


  /** @name Maps
   * The Maps
   */
  ///@{


//  UnorderedVariableOneToManyRelation&  m_nodeToFaceMap;
//  UnorderedVariableOneToManyRelation&  m_nodeToEdgeMap;

//  UnorderedVariableOneToManyRelation& m_toCrackSurfacesRelation;

  ///@}


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


  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;

  view_rtype_const<r1_array> referencePosition() const { return this->getData<r1_array>(viewKeys.referencePosition); }
  view_rtype<r1_array>       referencePosition()       { return this->getData<r1_array>(viewKeys.referencePosition); }
//  view_rtype_const<r1_array> totalDisplacement() const { return this->getData<r1_array>(viewKeys.totalDisplacement); }
//  view_rtype<r1_array>       totalDisplacement()       { return this->getData<r1_array>(viewKeys.totalDisplacement); }

  UnorderedVariableOneToManyRelation       & edgeList()       { return m_toEdgesRelation; }
  UnorderedVariableOneToManyRelation const & edgeList() const { return m_toEdgesRelation; }

  UnorderedVariableOneToManyRelation       & faceList()       { return m_toFacesRelation; }
  UnorderedVariableOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  array<localIndex_array>       & elementRegionList()       { return m_toElementRegionList; }
  array<localIndex_array> const & elementRegionList() const { return m_toElementRegionList; }

  array<localIndex_array>       & elementSubRegionList()       { return m_toElementSubRegionList; }
  array<localIndex_array> const & elementSubRegionList() const { return m_toElementSubRegionList; }

  UnorderedVariableOneToManyRelation       & elementList()       { return m_toElementList; }
  UnorderedVariableOneToManyRelation const & elementList() const { return m_toElementList; }


protected:

private:
  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                             localIndex_array const & packList ) const;


  /// copy constructor
  NodeManager() = delete;
  NodeManager( const NodeManager& init ) = delete;
  NodeManager& operator=( const NodeManager&) = delete;

  r1_array m_referencePosition;
  UnorderedVariableOneToManyRelation m_toEdgesRelation;
  UnorderedVariableOneToManyRelation m_toFacesRelation;

  array<localIndex_array> m_toElementRegionList ;
  array<localIndex_array> m_toElementSubRegionList ;
  UnorderedVariableOneToManyRelation m_toElementList ;



};
// *********************************************************************************************************************
// *********************************************************************************************************************


// *********************************************************************************************************************


}


#endif /* NODEMANAGERT_H_ */
