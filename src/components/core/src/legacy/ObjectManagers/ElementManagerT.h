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
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTMANAGERT_H_
#define ELEMENTMANAGERT_H_

#include "../../dataRepository/ManagedGroup.hpp"
#include "Common/Common.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "DataStructures/VectorFields/ElementRegionT.h"


/**
 * Class to manage the data stored at the element level.
 */
class ElementManagerT : public ObjectDataStructureBaseT
{
public:
  ElementManagerT();
  virtual ~ElementManagerT();

  void Initialize(  ){}

  void SetIsAttachedToSendingGhostNode( const NodeManager& nodeManager )
  {
    for( std::map< RegKeyType, ElementRegionT >::iterator i=m_ElementRegions.begin() ; i!=m_ElementRegions.end() ; ++i )
    {
      ElementRegionT& elementRegion = i->second;
      elementRegion.SetIsAttachedToSendingGhostNode(nodeManager);
    }
  }


  using ObjectDataStructureBaseT::resize;
  globalIndex resize( const lvector& numElements,
                      const array<string>& elementRegionNames,
                      const array<string>& elementTypes );


  void HACKInitialConditions();

  void SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject );
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL)
  {
    (void)referenceObject;
  }

  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<gArray1d>& objectToCompositionObject  )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("ElementManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }

  //void ConstructListOfBoundaryObjects( gArray1d& objectList ) const ;


  int AddElementToRegion( const int regionNumber, const int* const elemsToNodes );

  void ConstructListOfIndexesFromMap( const array< set< std::pair<ElementRegionT*,localIndex> > >& toElementMap,
                                      const lArray1d& nodeList,
                                      std::map< std::string, lArray1d>& localIndexes,
                                      const int depth );

  void ModifyToElementMapsFromSplit( const std::map< std::string, lSet>& modifiedElements,
                                     NodeManager& nodeManager,
                                     FaceManager& faceManager );

  void UpdateExternalityFromSplit( const std::map< std::string, lSet>& modifiedElements,
                                   NodeManager& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManager& faceManager );

  void InitializeFlowFaceRegion();
  void GenerateFlowFaceRegion(FaceManager& faceManager);

  using ObjectDataStructureBaseT::WriteSilo;
  using ObjectDataStructureBaseT::ReadSilo;

  void WriteSilo( SiloFile& siloFile,
                  const std::string& meshname,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const std::string& regionName = "none",
                  const lArray1d& mask = lArray1d() );

  void ReadSilo( const SiloFile& siloFile,
                 const std::string& meshname,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart );


  void ResetGlobalToLocalMap();

  template< typename T_indices >
  unsigned int PackElements( bufvector& buffer,
                             lSet& sendnodes,
                             lSet& sendfaces,
                             const std::map<std::string,T_indices>& elementList,
                             const NodeManager& nodeManager,
                             const FaceManager& faceManager,
                             const bool packConnectivityToGlobal,
                             const bool packFields,
                             const bool packMaps,
                             const bool packSets ) const;

  unsigned int UnpackElements( const bufvector& buffer,
                               const NodeManager& nodeManager,
                               const FaceManager& faceManager,
                               std::map< std::string, lArray1d>& elementRegionReceiveLocalIndices,
                               const bool unpackConnectivityToLocal,
                               const bool unpackFields,
                               const bool unpackMaps,
                               const bool unpackSets );

  unsigned int UnpackElements( const char*& pbuffer,
                               const NodeManager& nodeManager,
                               const FaceManager& faceManager,
                               std::map< std::string, lArray1d>& elementRegionReceiveLocalIndices,
                               const bool unpackConnectivityToLocal,
                               const bool unpackFields,
                               const bool unpackMaps,
                               const bool unpackSets );

  void ConnectivityFromGlobalToLocal( const std::map< std::string, lSet>& allReceivedElements,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );

  localIndex* ElementToNodeMap( const localIndex elemIndex )
  {
    const RegKeyType regionKey = m_regionIndexToKey[m_ElementIndexToRegionIndex[elemIndex]];
    return m_ElementRegions[regionKey].m_toNodesRelation[m_ElementIndexToRegionLocalIndex[elemIndex]];
  }

  void ElementToNodeMap( const localIndex elemIndex, lArray1d& indices )
  {
    indices.clear();
    const RegKeyType regionKey = m_regionIndexToKey[m_ElementIndexToRegionIndex[elemIndex]];
    const localIndex* iptr = m_ElementRegions[regionKey].m_toNodesRelation[m_ElementIndexToRegionLocalIndex[elemIndex]];
    indices.resize(m_ElementRegions[regionKey].m_numNodesPerElem);
    localIndex a = 0;
    for(lArray1d::iterator it = indices.begin() ; it != indices.end() ; ++it, ++a)
      *it = iptr[a];
  }

  int m_numElems;

  void ResizeNumberOfElementsAfterSplit();

  array<lArray1d>& m_ElementToNodeMap;
  array<lArray1d>& m_ElementToFaceMap;
  array<lArray1d>& m_ElementToElementMap;

  array<integer>& m_ElementIndexToRegionIndex;
  lArray1d& m_ElementIndexToRegionLocalIndex;

  array<string> m_regionIndexToKey;


  typedef std::string RegKeyType;
  std::map< RegKeyType, ElementRegionT > m_ElementRegions;

private:
  ElementManagerT( const ElementManagerT& );
  ElementManagerT& operator=( const ElementManagerT&);


};

#endif /* ELEMENTMANAGERT_H_ */
