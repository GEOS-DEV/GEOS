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
 * @file PhysicalDomainT.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#include "PhysicalDomainT.h"
//#include "ElementLibrary/IntegrationRuleT.h"
#include "IO/BinStream.h"
#include "IO/silo/SiloFile.h"
#include "BoundaryConditions/BoundaryConditions.h"
#include "ElementLibrary/FiniteElement.h"


PhysicalDomainT::PhysicalDomainT():
  m_globalDomainNumber(),
  m_ellipsoidalDiscreteElementManager(),
  m_ellipsoidalContactManager(),
  m_discreteElementSurfaceFaces(),
  m_discreteElementSurfaceNodes(),
  m_discreteElementManager(&m_discreteElementSurfaceNodes, &m_discreteElementSurfaceFaces),
  m_contactManager(),
  m_feNodeManager(),
  m_feElementManager(),
  m_feFaceManager(),
  m_crackManager(NULL),
  m_crackSurfaceVertex(NULL),
  m_externalFaces(&m_feFaceManager, &m_discreteElementSurfaceFaces),
  m_cartesianGridManager()
#ifdef SRC_EXTERNAL
  ,
  m_faultPatchFaces(),
  m_faultPatchNodes(),
  m_faultPorePressure(),
  m_faultElementManager(&m_faultPatchNodes, &m_faultPatchFaces),
  m_wellboreManager()
#endif
#ifdef SRC_INTERNAL2
  ,
  m_microseismicElementManager(),
  m_jointSets2(),
  m_xfemManager(NULL)
#endif
{
  /*
     m_managedObjectDataStructures.resize(Last_ObjectDataStructureNames_Index +
        1);
     m_managedObjectDataStructures[EllipsoidalDiscreteElementManager] =
        &m_ellipsoidalDiscreteElementManager;
     m_managedObjectDataStructures[EllipsoidalContactManager] =
        &m_ellipsoidalContactManager;
     m_managedObjectDataStructures[DiscreteElementFaceManager] =
        &m_discreteElementSurfaceFaces;
     m_managedObjectDataStructures[DiscreteElementNodeManager] =
        &m_discreteElementSurfaceNodes;
     m_managedObjectDataStructures[DiscreteElementManager] =
        &m_discreteElementManager;
     m_managedObjectDataStructures[ContactManager] = &m_contactManager;
     m_managedObjectDataStructures[FiniteElementNodeManager] = &m_feNodeManager;
     m_managedObjectDataStructures[FiniteElementElementManager] =
        &m_feElementManager;
     m_managedObjectDataStructures[FiniteElementFaceManager] = &m_feFaceManager;
     m_managedObjectDataStructures[FiniteElementEdgeManager] = &m_feEdgeManager;
     m_managedObjectDataStructures[ExternalFaceManager] = &m_externalFaces;
     m_managedObjectDataStructures[VirtualFaceManager] = &m_virtualFaceManager;
     m_managedObjectDataStructures[VirtualEdgeManager] = &m_virtualEdgeManager;
   */

  m_feNodeManager.m_nodeToEdgeMap.SetRelatedObject( &m_feEdgeManager );


}

PhysicalDomainT::~PhysicalDomainT()
{}

int PhysicalDomainT::Initialize()
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator elementRegionIter =
         m_feElementManager.m_ElementRegions.begin() ; elementRegionIter
       != m_feElementManager.m_ElementRegions.end() ; ++elementRegionIter)
  {
    const std::string& RegionKey = elementRegionIter->first;
    ElementRegionT& elementRegion = elementRegionIter->second;

    std::cout << "Rank, Region Key, NumElems: " << rank << ", " << RegionKey << ", "
              << elementRegion.m_numElems << std::endl;

    //TODO: change this to be determined based on solver type.  We only use
    // these arrays for small def solvers.
//    if ( !elementRegion.m_elementGeometryID.compare(0, 3, "S4R") )
//    {
//      //Do nothing, since we don't use any mechanical solver for this shell
// element.
//    }
    /*
       else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "TRSH") )
       {
       //Do nothing, since we don't use any mechanical solver for this shell
          element.
       }
     */
//    else
//    {
    elementRegion.CalculateShapeFunctionDerivatives(m_feNodeManager);
//    }

    elementRegion.Initialize();

    // make sets from nodesets
    for (std::map<std::string, lSet>::const_iterator i = m_feNodeManager.m_Sets.begin() ; i
         != m_feNodeManager.m_Sets.end() ; ++i)
    {
      const std::string& setname = i->first;
      const lSet& set = i->second;
      //std::cout << "setname = " << setname << "set = " << set;
      elementRegion.ConstructSetFromSetAndMap(set, elementRegion.m_toNodesRelation, setname);
    }
  }

#ifdef SRC_EXTERNAL
  m_faultElementManager.resize(m_faultPatchFaces.DataLengths());
#endif

  SetInterObjectRelations();

  return 0;
}


void PhysicalDomainT::SetInterObjectRelations()
{

  m_ellipsoidalDiscreteElementManager.m_neighborList.SetRelatedObject( &m_ellipsoidalDiscreteElementManager );
  m_ellipsoidalDiscreteElementManager.m_neighborListInverse.SetRelatedObject( &m_ellipsoidalDiscreteElementManager );

//  m_ellipsoidalContactManager;

  m_discreteElementSurfaceFaces.m_toNodesRelation.SetRelatedObject( &m_discreteElementSurfaceNodes );
  m_discreteElementSurfaceNodes.m_nodeToFaceMap.SetRelatedObject( &m_discreteElementSurfaceFaces );
  m_discreteElementManager.m_discreteElementToExternalFacesMap.SetRelatedObject( &m_discreteElementSurfaceFaces );
  m_discreteElementManager.m_discreteElementToExternalNodesMap.SetRelatedObject( &m_discreteElementSurfaceNodes );

  m_contactManager.m_contactToIntersectionPolygonPointsMap.SetRelatedObject( &m_externalFaces );

#ifdef SRC_EXTERNAL
  m_faultPatchFaces.m_toNodesRelation.SetRelatedObject( &m_faultPatchNodes );
  m_faultPatchNodes.m_nodeToFaceMap.SetRelatedObject( &m_faultPatchFaces );
#endif

  m_feNodeManager.m_nodeToEdgeMap.SetRelatedObject( &m_feEdgeManager );
  m_feNodeManager.m_nodeToFaceMap.SetRelatedObject( &m_feFaceManager );

  m_feFaceManager.m_toEdgesRelation.SetRelatedObject( &m_feEdgeManager );
  m_feFaceManager.m_toNodesRelation.SetRelatedObject( &m_feNodeManager );

  m_feEdgeManager.m_toNodesRelation.SetRelatedObject( &m_feNodeManager );
  m_feEdgeManager.m_toFacesRelation.SetRelatedObject( &m_feFaceManager );

  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegion=m_feElementManager.m_ElementRegions.begin() ;
       elementRegion!=m_feElementManager.m_ElementRegions.end() ; ++elementRegion )
  {
    elementRegion->second.m_toNodesRelation.SetRelatedObject( &m_feNodeManager );
    //elementRegion->second.m_toEdgesRelation.SetRelatedObject( &m_feEdgeManager
    // );
    elementRegion->second.m_toFacesRelation.SetRelatedObject( &m_feFaceManager );
  }
}


/**
 * @brief Get the map of all enumerated and managed ObjectDataStructureBase
 * managers
 * @author Scott Johnson
 */
/*const array<ObjectDataStructureBaseT*>&
   PhysicalDomainT::GetObjectDataStructures()
   {
   return m_managedObjectDataStructures;
   }*/



template< typename T_indices >
void PhysicalDomainT::Pack( const ObjectDataStructureKeys key,
                            const T_indices& sendIndices,
                            bufvector& buffer,
                            const bool packConnectivityToGlobal,
                            const bool packFields,
                            const bool packMaps,
                            const bool packSets )
{
  switch (key)
  {
  case PhysicalDomainT::FiniteElementNodeManager:
    m_feNodeManager.PackNodes( sendIndices, m_feFaceManager, buffer, packConnectivityToGlobal, packFields, packMaps, packSets );
    break;
  case PhysicalDomainT::FiniteElementEdgeManager:
    m_feEdgeManager.PackEdges( sendIndices,
                               m_feNodeManager,
                               m_feFaceManager,
                               buffer,
                               packConnectivityToGlobal,
                               packFields, packMaps, packSets );
    break;
  case PhysicalDomainT::FiniteElementFaceManager:
    m_feFaceManager.PackFaces( sendIndices,
                               m_feNodeManager,
                               &m_feEdgeManager,
                               buffer,
                               packConnectivityToGlobal,
                               packFields, packMaps, packSets );
    break;
  case PhysicalDomainT::EllipsoidalDiscreteElementManager:
  //m_ellipsoidalDiscreteElementManager.Pack( sendIndices, buffer );
  //break;
  case PhysicalDomainT::DiscreteElementManager:
  //m_discreteElementManager.Pack( sendIndices, buffer );
  //break;
  case PhysicalDomainT::DiscreteElementFaceManager:
  //m_discreteElementSurfaceFaces.PackFaces( sendIndices,
  //                                         m_feNodeManager,
  //                                         &m_feEdgeManager,
  //                                         buffer,
  //                                         packConnectivityToGlobal,
  //                                         packFields, packMaps, packSets );
  //break;
  case PhysicalDomainT::DiscreteElementNodeManager:
  //m_discreteElementSurfaceNodes.PackNodes( sendIndices,
  //                                         m_discreteElementSurfaceFaces,
  //                                         buffer,
  //                                         packConnectivityToGlobal,
  //                                         packFields, packMaps, packSets );
  //break;
  case PhysicalDomainT::FaultPatchNodeManager:
  //m_faultPatchNodes.PackNodes( sendIndices, m_faultPatchFaces, buffer,
  // packConnectivityToGlobal, packFields, packMaps, packSets );
  //break;
  case PhysicalDomainT::FaultPatchFaceManager:
  //m_faultPatchFaces.PackFaces( sendIndices,
  //                             m_faultPatchNodes,
  //                             &m_feEdgeManager,
  //                             buffer,
  //                             packConnectivityToGlobal,
  //                             packFields );
  //break;
  case PhysicalDomainT::FaultPatchElementManager:
  case PhysicalDomainT::EllipsoidalContactManager:
  case PhysicalDomainT::ContactManager:
  case PhysicalDomainT::FiniteElementElementManager:
  case PhysicalDomainT::FiniteElementElementRegion:
  case PhysicalDomainT::ExternalFaceManager:
  case PhysicalDomainT::CartesianGridManager:
  case PhysicalDomainT::VirtualFaceManager:
  case PhysicalDomainT::VirtualEdgeManager:
  case PhysicalDomainT::numObjectDataStructureNames:
  default:
    throw GPException("PhysicalDomainT::GetObjectDataStructure: inappropriate type for PhysicalDomainT::Pack "+ toString<int>(key));
    break;

  }
}
template void PhysicalDomainT::Pack( const ObjectDataStructureKeys, const lArray1d&, bufvector&, const bool, const bool, const bool, const bool );
template void PhysicalDomainT::Pack( const ObjectDataStructureKeys, const lSet&, bufvector&, const bool, const bool,const bool, const bool );
template<> void PhysicalDomainT::Pack( const ObjectDataStructureKeys key,
                                       const std::map<std::string,lArray1d>& sendIndices,
                                       bufvector& buffer,
                                       const bool packConnectivityToGlobal,
                                       const bool packFields,
                                       const bool packMaps,
                                       const bool packSets )
{
  if( key == PhysicalDomainT::FiniteElementElementManager )
  {
    lSet junk1;
    lSet junk2;

    m_feElementManager.PackElements( buffer,
                                     junk1, junk2,
                                     sendIndices,
                                     m_feNodeManager,
                                     m_feFaceManager,
                                     packConnectivityToGlobal,
                                     packFields,
                                     packMaps,
                                     packSets );
  }
  else
  {
    throw GPException("PhysicalDomainT::GetObjectDataStructure: inappropriate type for PhysicalDomainT::Pack "+ toString<int>(key));
  }


}


/**
 * @brief Unpack the requested data structure from the buffer
 * @author Scott Johnson
 */

template< typename T_indices >
void PhysicalDomainT::Unpack(const ObjectDataStructureKeys name,
                             const char*& pbuffer,
                             T_indices& indices,
                             const bool unpackConnectivityToLocal,
                             const bool unpackFields,
                             const bool unpackMaps,
                             const bool unpackSets  )
{
  switch (name)
  {
  case PhysicalDomainT::EllipsoidalDiscreteElementManager:
    m_ellipsoidalDiscreteElementManager.Unpack( pbuffer, indices );
    break;
  case PhysicalDomainT::DiscreteElementManager:
    m_discreteElementManager.Unpack( pbuffer, indices );
    break;
  case PhysicalDomainT::DiscreteElementFaceManager:
    m_discreteElementSurfaceFaces.UnpackFaces( pbuffer,
                                               m_discreteElementSurfaceNodes,
                                               NULL,
                                               NULL,
                                               indices,
                                               unpackConnectivityToLocal,
                                               unpackFields, unpackMaps, unpackSets );
    break;
  case PhysicalDomainT::DiscreteElementNodeManager:
    m_discreteElementSurfaceNodes.UnpackNodes( pbuffer,
                                               m_discreteElementSurfaceFaces,
                                               indices,
                                               unpackConnectivityToLocal,
                                               unpackFields, unpackMaps, unpackSets );
    break;
  case PhysicalDomainT::FiniteElementNodeManager:
    m_feNodeManager.UnpackNodes( pbuffer,
                                 m_feFaceManager,
                                 indices,
                                 unpackConnectivityToLocal,
                                 unpackFields, unpackMaps, unpackSets );
    break;
  case PhysicalDomainT::FiniteElementFaceManager:
    m_feFaceManager.UnpackFaces( pbuffer,
                                 m_feNodeManager,
                                 &m_feEdgeManager,
                                 &m_externalFaces,
                                 indices,
                                 unpackConnectivityToLocal,
                                 unpackFields, unpackMaps, unpackSets );
    break;
  case PhysicalDomainT::FiniteElementEdgeManager:
    m_feEdgeManager.UnpackEdges( pbuffer,
                                 m_feNodeManager,
                                 m_feFaceManager,
                                 indices,
                                 unpackConnectivityToLocal,
                                 unpackFields, unpackMaps, unpackSets );
    break;
#ifdef SRC_EXTERNAL
  case PhysicalDomainT::FaultPatchFaceManager:
    m_faultPatchFaces.UnpackFaces( pbuffer,
                                   m_faultPatchNodes,
                                   NULL,
                                   NULL,
                                   indices,
                                   unpackConnectivityToLocal,
                                   unpackFields, unpackMaps, unpackSets );
    break;
  case PhysicalDomainT::FaultPatchNodeManager:
    m_faultPatchNodes.UnpackNodes( pbuffer,
                                   m_faultPatchFaces,
                                   indices,
                                   unpackConnectivityToLocal,
                                   unpackFields, unpackMaps, unpackSets );
    break;
  case PhysicalDomainT::FaultPatchElementManager:
#endif
  case PhysicalDomainT::FiniteElementElementManager:
  case PhysicalDomainT::EllipsoidalContactManager:
  case PhysicalDomainT::ContactManager:
  case PhysicalDomainT::FiniteElementElementRegion:
  case PhysicalDomainT::ExternalFaceManager:
  case PhysicalDomainT::CartesianGridManager:
  case PhysicalDomainT::VirtualFaceManager:
  case PhysicalDomainT::VirtualEdgeManager:
  case PhysicalDomainT::numObjectDataStructureNames:
  default:
    throw GPException("PhysicalDomainT::GetObjectDataStructure: Unrecognized object type "
                      + toString<int>(name));
  }
}
template void PhysicalDomainT::Unpack( const ObjectDataStructureKeys, const char*&, lArray1d&, const bool, const bool, const bool, const bool );
template<> void PhysicalDomainT::Unpack( const ObjectDataStructureKeys name,
                                         const char*& pbuffer,
                                         std::map< std::string, lArray1d>& indices,
                                         const bool unpackConnectivityToLocal,
                                         const bool unpackFields,
                                         const bool unpackMaps,
                                         const bool unpackSets  )
{
  if( name == PhysicalDomainT::FiniteElementElementManager )
  {
    std::map< std::string, lArray1d> elementRegionReceiveLocalIndices;
    m_feElementManager.UnpackElements( pbuffer,
                                       m_feNodeManager,
                                       m_feFaceManager,
                                       elementRegionReceiveLocalIndices,
                                       unpackConnectivityToLocal,
                                       unpackFields, unpackMaps, unpackSets );

    for( std::map< std::string, lArray1d>::iterator i =elementRegionReceiveLocalIndices.begin() ;
         i!=elementRegionReceiveLocalIndices.end() ; ++i )
    {
      indices[i->first].insert( indices[i->first].end(), i->second.begin(), i->second.end() );
    }
  }
  else
  {
    throw GPException("PhysicalDomainT::GetObjectDataStructure: inappropriate type for PhysicalDomainT::Unpack "
                      + toString<int>(name));

  }

}


/**
 * @brief Get the associated data structure for a particular enumerated name
 * @author Scott Johnson
 * @param[in] name Enumerated name
 * @return ObjectDataStructureBase manager associated with the given enumerated
 * name
 */
ObjectDataStructureBaseT& PhysicalDomainT::GetObjectDataStructure( PhysicalDomainT::ObjectDataStructureKeys key,
                                                                   const std::string& regionName )
{
  switch (key)
  {
  case EllipsoidalDiscreteElementManager:
    return m_ellipsoidalDiscreteElementManager;
  case EllipsoidalContactManager:
    return m_ellipsoidalContactManager;
  case DiscreteElementNodeManager:
    return m_discreteElementSurfaceNodes;
  case DiscreteElementFaceManager:
    return m_discreteElementSurfaceFaces;
  case DiscreteElementManager:
    return m_discreteElementManager;
  case ContactManager:
    return m_contactManager;
  case FiniteElementNodeManager:
    return m_feNodeManager;
  case FiniteElementEdgeManager:
    return m_feEdgeManager;
  case FiniteElementFaceManager:
    return m_feFaceManager;
  case FiniteElementElementManager:
    return m_feElementManager;
  case FiniteElementElementRegion:
    return m_feElementManager.m_ElementRegions[regionName];
  case ExternalFaceManager:
    return m_externalFaces;
  case CartesianGridManager:
    return m_cartesianGridManager;
#ifdef SRC_EXTERNAL
  case FaultPatchNodeManager:
    return m_faultPatchNodes;
  case FaultPatchFaceManager:
    return m_faultPatchFaces;
  case FaultPatchElementManager:
    return m_faultElementManager;
#endif
  case VirtualEdgeManager:
  case VirtualFaceManager:
  case numObjectDataStructureNames:
  case WellboreElementManager:
    return m_wellboreManager;
  default:
    throw GPException("PhysicalDomainT::GetObjectManager: Unrecognized object type "+ toString<int> (key));
  }

}

/**
 * @brief Get the associated numerator for a particular managed data structure
 * @author Scott Johnson
 * @param[in] name Name of the managed data structure
 * @return Enumerated name
 */
PhysicalDomainT::ObjectDataStructureKeys PhysicalDomainT::GetObjectDataStructureConditionKey(const std::string name)
{
  if(name.compare( PhysicalDomainT::EllipsoidalConditionStr() ) == 0)
    return PhysicalDomainT::EllipsoidalDiscreteElementManager;
  else if(name.compare( PhysicalDomainT::DiscreteElementConditionStr() ) == 0)
    return PhysicalDomainT::DiscreteElementManager;
  else if(name.compare( PhysicalDomainT::FiniteElementNodeConditionStr() ) == 0)
    return PhysicalDomainT::FiniteElementNodeManager;
  else if(name.compare( PhysicalDomainT::FiniteElementElementConditionStr() ) == 0)
    return PhysicalDomainT::FiniteElementElementRegion;
  else if(name.compare( PhysicalDomainT::FiniteElementFaceConditionStr() ) == 0)
    return PhysicalDomainT::FiniteElementFaceManager;
  else if(name.compare( PhysicalDomainT::FiniteElementEdgeConditionStr() ) == 0)
    return PhysicalDomainT::FiniteElementEdgeManager;
  else if(name.compare( PhysicalDomainT::ExternalFaceConditionStr() ) == 0)
    return PhysicalDomainT::ExternalFaceManager;
  else if(name.compare( PhysicalDomainT::FaultElementConditionStr() ) == 0)
    return PhysicalDomainT::FaultPatchElementManager;
  else if(name.compare( PhysicalDomainT::FaultFaceManagerStr() ) == 0)
    return PhysicalDomainT::FaultPatchFaceManager;
  else if(name.compare( PhysicalDomainT::CartesianGridManagerStr() ) == 0)
    return PhysicalDomainT::CartesianGridManager;
  else if(name.compare( PhysicalDomainT::WellboreElementManagerStr() ) == 0)
    return PhysicalDomainT::WellboreElementManager;
  else
    throw GPException("PhysicalDomainT::GetObjectDataStructureConditionKey: Unrecognized object type "
                      + name);
}


/**
 * @brief Get the associated numerator for a particular managed data structure
 * @author Scott Johnson
 * @param[in] name Name of the managed data structure
 * @return Enumerated name
 */
PhysicalDomainT::ObjectDataStructureKeys PhysicalDomainT::GetObjectDataStructureKey(const std::string name)
{
  if(name.compare( PhysicalDomainT::EllipsoidalDiscreteElementManagerStr() ) == 0)
    return PhysicalDomainT::EllipsoidalDiscreteElementManager;
  else if(name.compare( PhysicalDomainT::EllipsoidalContactManagerStr() ) == 0)
    return PhysicalDomainT::EllipsoidalContactManager;
  else if(name.compare( PhysicalDomainT::DiscreteElementFaceManagerStr() ) == 0)
    return PhysicalDomainT::DiscreteElementFaceManager;
  else if(name.compare( PhysicalDomainT::DiscreteElementNodeManagerStr() ) == 0)
    return PhysicalDomainT::DiscreteElementNodeManager;
  else if(name.compare( PhysicalDomainT::DiscreteElementManagerStr() ) == 0)
    return PhysicalDomainT::DiscreteElementManager;
  else if(name.compare( PhysicalDomainT::ContactManagerStr() ) == 0)
    return PhysicalDomainT::ContactManager;
  else if(name.compare( PhysicalDomainT::FiniteElementNodeManagerStr() ) == 0)
    return PhysicalDomainT::FiniteElementNodeManager;
  else if(name.compare( PhysicalDomainT::FiniteElementElementManagerStr() ) == 0)
    return PhysicalDomainT::FiniteElementElementManager;
  else if(name.compare( PhysicalDomainT::FiniteElementElementRegionStr() ) == 0)
    return PhysicalDomainT::FiniteElementElementRegion;
  else if(name.compare( PhysicalDomainT::FiniteElementFaceManagerStr() ) == 0)
    return PhysicalDomainT::FiniteElementFaceManager;
  else if(name.compare( PhysicalDomainT::FiniteElementEdgeManagerStr() ) == 0)
    return PhysicalDomainT::FiniteElementEdgeManager;
  else if(name.compare( PhysicalDomainT::ExternalFaceManagerStr() ) == 0)
    return PhysicalDomainT::ExternalFaceManager;
  else if(name.compare( PhysicalDomainT::CartesianGridManagerStr() ) == 0)
    return PhysicalDomainT::CartesianGridManager;
  else if(name.compare( PhysicalDomainT::FaultElementManagerStr() ) == 0)
    return PhysicalDomainT::FaultPatchElementManager;
  else if(name.compare( PhysicalDomainT::FaultFaceManagerStr() ) == 0)
    return PhysicalDomainT::FaultPatchFaceManager;
  else if(name.compare( PhysicalDomainT::FaultNodeManagerStr() ) == 0)
    return PhysicalDomainT::FaultPatchNodeManager;
  else if(name.compare( PhysicalDomainT::VirtualEdgeManagerStr() ) == 0)
    return PhysicalDomainT::VirtualEdgeManager;
  else if(name.compare( PhysicalDomainT::VirtualFaceManagerStr() ) == 0)
    return PhysicalDomainT::VirtualFaceManager;
  else if(name.compare( PhysicalDomainT::WellboreElementManagerStr() ) == 0)
    return PhysicalDomainT::WellboreElementManager;
  else
    throw GPException("PhysicalDomainT::GetObjectDataStructureKey: Unrecognized object type "
                      + name);
}

std::string PhysicalDomainT::GetObjectDataStructureName(const PhysicalDomainT::ObjectDataStructureKeys key)
{
  std::string names[] = { EllipsoidalDiscreteElementManagerStr(),
                          EllipsoidalContactManagerStr(),
                          DiscreteElementNodeManagerStr(),
                          DiscreteElementFaceManagerStr(),
                          DiscreteElementManagerStr(),
                          ContactManagerStr(),
                          FiniteElementNodeManagerStr(),
                          FiniteElementEdgeManagerStr(),
                          FiniteElementFaceManagerStr(),
                          FiniteElementElementManagerStr(),
                          FiniteElementElementRegionStr(),
                          ExternalFaceManagerStr(),
                          CartesianGridManagerStr(),
                          FaultNodeManagerStr(),
                          FaultFaceManagerStr(),
                          FaultElementManagerStr(),
                          VirtualEdgeManagerStr(),
                          VirtualFaceManagerStr(),
                          WellboreElementManagerStr() };
  return names[key];
}
/*
 * @author S. Walsh
 *
 * Get an unknown (at compile time) object manager at runtime.
 * If the desired object manager is known ahead of time the object itself should
 * be used instead.
 *
 * Returns the object data structure corresponding to the requested type (and
 * region for element regions).
 * Should only be used when object type is not known at compile time (eg.
 * initilization).
 *
 * @param objectType - type of manager to return
 * @param regionName - used to get named element region - otherwise ignored
 *
 *
   ObjectDataStructureBaseT& PhysicalDomainT::GetObjectManager(const
 * ObjectManagerType objectType,
                                                            const std::string&
 * regionName)
   {
   //TODO: this should all be switched out for the logic above;
 * ObjectManagerType is ambiguous for this purpose!
   switch (objectType)
   {
    case ObjectDataStructureBaseT::FiniteElementElementManager:
      return m_feElementManager;
    case ObjectDataStructureBaseT::FiniteElementElementRegion:
      return m_feElementManager.m_ElementRegions[regionName];
    case ObjectDataStructureBaseT::FiniteElementNodeManager:
      return m_feNodeManager;
    case ObjectDataStructureBaseT::FiniteElementFaceManager:
      return m_feFaceManager;
    case ObjectDataStructureBaseT::FiniteElementEdgeManager:
      return m_feEdgeManager;
    case ObjectDataStructureBaseT::ExternalFaceManager:
      return m_externalFaces;
    case ObjectDataStructureBaseT::DiscreteElementManager:
      return m_discreteElementManager;
    case ObjectDataStructureBaseT::EllipsoidalDiscreteElementManager:
      return m_ellipsoidalDiscreteElementManager;
    case ObjectDataStructureBaseT::EllipsoidalContactManager:
      return m_ellipsoidalContactManager;
    default:
      throw GPException("PhysicalDomainT::GetObjectManager: Unrecognized object
 * type "+ toString<int> (objectType));

   }

   return m_feElementManager;
   }
 */
/**
 * @brief Write to SILO
 * @author R Settgast
 * Write all objects to SILO format
 * @param[out] siloFile SILO file object
 * @param[in] cycleNum Timestep index
 * @param[in] problemTime Current simulation time
 * @param[in] writeFEMMesh Flag whether to write out the finite element mesh
 * @param[in] writeFEMFaces Flag whether to write out the finite element faces
 * @param[in] writeFEMEdges Flag whether to write out the finite element edges
 * @param[in] writeDE Flag whether to write out the discrete elements
 * @param[in] writeCP Flag whether to write out the common plane contacts
 * @param[in] writeCG Flag whether to write out the cartesian grid
 */
void PhysicalDomainT::WriteSilo(SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart,
                                const bool WriteXFEM,
                                const bool writeFEMMesh,
                                const bool writeFEMFaces,
                                const bool writeFEMEdges,
                                const bool writeDE,
                                const bool writeCP,
                                const bool writeCG,
                                const bool writeFaultElements,
                                const bool writeMicroseismicSources)
{
  //--------------initialize-----------------
  if(WriteXFEM)
  {
    WriteXFEMElements( siloFile, cycleNum, problemTime, isRestart, true );
  }

  WriteFiniteElementMesh( siloFile, cycleNum, problemTime, isRestart, writeFEMMesh, writeFEMFaces, writeFEMEdges);

  WriteDiscreteElements( siloFile, cycleNum, problemTime, isRestart, writeDE );

  WriteEllipsoidalDiscreteElements( siloFile, cycleNum, problemTime, isRestart, writeDE );

#ifdef SRC_EXTERNAL
  // WriteMicroseismicSources( siloFile, cycleNum, problemTime, isRestart,
  // writeMicroseismicSources );
  WriteFaultElementMesh( siloFile, cycleNum, problemTime, isRestart, writeFaultElements);
#endif
#ifdef SRC_INTERNAL2
  m_microseismicElementManager.WriteMicroseismicSilo(siloFile, cycleNum, problemTime, isRestart, writeMicroseismicSources);
#endif
  WriteCommonPlanes( siloFile, cycleNum, problemTime, isRestart, writeCP );

  WriteCartesianGrid( siloFile, cycleNum, problemTime, isRestart, writeCG );
#ifdef SRC_EXTERNAL
  m_wellboreManager.WriteWellboreSilo( siloFile, cycleNum, problemTime, isRestart );
#endif


  //--------------clean up-----------------


  if( isRestart )
  {
    siloFile.DBWriteWrapper("m_globalDomainNumber",m_globalDomainNumber);
  }

}

void PhysicalDomainT::ReadSilo( const SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart,
                                const bool writeFEMMesh,
                                const bool writeFEMFaces,
                                const bool writeFEMEdges,
                                const bool writeDE,
                                const bool writeCP,
                                const bool writeCG,
                                const bool writeFaultElements,
                                const bool writeMicroseismicSources )
{

  ReadFiniteElementMesh( siloFile, cycleNum, problemTime, isRestart );

  ReadDiscreteElements( siloFile, cycleNum, problemTime, isRestart );

  ReadEllipsoidalDiscreteElements( siloFile, cycleNum, problemTime, isRestart );

#ifdef SRC_EXTERNAL
  ReadFaultElementMesh( siloFile, cycleNum, problemTime, isRestart );
#endif

#ifdef SRC_INTERNAL2
  m_microseismicElementManager.m_elemManagerHACK = &m_feElementManager;
  m_microseismicElementManager.ReadMicroseismicSilo(siloFile, cycleNum, problemTime, isRestart);
#endif

  ReadCommonPlanes( siloFile, cycleNum, problemTime, isRestart );

  ReadCartesianGrid( siloFile, cycleNum, problemTime, isRestart );

#ifdef SRC_EXTERNAL
  m_wellboreManager.ReadSilo( siloFile, "WellboreFields", "wellbore_mesh",
                              DB_NODECENT, cycleNum, problemTime, isRestart );
#endif
}

#ifdef SRC_EXTERNAL
void PhysicalDomainT::WriteFaultElementMesh(SiloFile& siloFile,
                                            const int cycleNum,
                                            const realT problemTime,
                                            const bool isRestart,
                                            const bool writeFaultElements)
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE FAULT PATCH DATA-----------------
  if (writeFaultElements && m_faultPatchFaces.DataLengths() > 0)
  {
    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(m_faultPatchNodes.m_numNodes);
    dvector ycoords(m_faultPatchNodes.m_numNodes);
    dvector zcoords(m_faultPatchNodes.m_numNodes);

    const array<R1Tensor>& refpos = m_faultPatchNodes.GetFieldData<FieldInfo::referencePosition>();
    const array<R1Tensor>& disp = m_faultPatchNodes.GetFieldData<FieldInfo::displacement>();

    for (localIndex a = 0 ; a < m_faultPatchNodes.m_numNodes ; ++a)
    {
      R1Tensor nodePosition;
      nodePosition = refpos[a];
      nodePosition += disp[a];

      xcoords[a] = nodePosition(0);
      ycoords[a] = nodePosition(1);
      zcoords[a] = nodePosition(2);
    }

    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();

    // TODO this only works for quad faces!!!
    // face mesh

    const int numFaceTypes = 1;
    array<localIndex*> faceConnectivity(numFaceTypes);
    array<globalIndex*> globalFaceNumbers(numFaceTypes);
    ivector fshapecnt(numFaceTypes);
    ivector fshapetype(numFaceTypes);
    ivector fshapesize(numFaceTypes);

    array<lArray2d> faceToNodeMap(numFaceTypes);

    for (int faceType = 0 ; faceType < numFaceTypes ; ++faceType)
    {
      faceToNodeMap[faceType].resize2(m_faultPatchFaces.m_numFaces, 4);

      for (localIndex k = 0 ; k < m_faultPatchFaces.m_numFaces ; ++k)
      {
        for (int a = 0 ; a < 4 ; ++a)
        {
          faceToNodeMap[faceType][k][a] = m_faultPatchFaces.m_toNodesRelation[k][a];
        }
      }

      faceConnectivity[faceType] = faceToNodeMap[faceType].data();

      globalFaceNumbers[faceType] = m_faultPatchFaces.m_localToGlobalMap.data();
      fshapecnt[faceType] = m_faultPatchFaces.m_numFaces;
      fshapetype[faceType] = DB_ZONETYPE_QUAD;
      fshapesize[faceType] = 4;
    }


    const std::string meshName("fault_mesh");
    siloFile.WriteMeshObject(meshName, m_faultPatchNodes.m_numNodes, coords,
                             m_faultPatchNodes.m_localToGlobalMap.data(), numFaceTypes,
                             fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                             NULL, fshapetype.data(), fshapesize.data(), cycleNum, problemTime);

    //-----------------------------
    //Write node data
    //-----------------------------

    m_faultPatchNodes.WriteSilo( siloFile, "FaultNodeFields", meshName, DB_NODECENT, cycleNum, problemTime, isRestart );


    //-----------------------------
    //Write face data
    //-----------------------------
    m_faultPatchFaces.WriteSilo( siloFile, "FaultFaceFields", meshName, DB_ZONECENT, cycleNum, problemTime, isRestart );

    //-----------------------------
    //Write element data
    //-----------------------------
    m_faultElementManager.WriteSiloFaultElements( siloFile, "FaultElementFields", meshName, DB_ZONECENT, cycleNum, problemTime, isRestart );
  }

  //Write pore pressure file
  if (writeFaultElements && m_faultElementManager.HasPorePressure())
  {
    //get the grid(s) of pore pressures for the current time step
    const Table4D* ppptr = m_faultElementManager.PorePressure();
    std::vector<realT> x0[3];
    localIndex nXYZ = 1;
    lArray1d nX(3);
    for(localIndex i = 0 ; i < 3 ; i++)
    {
      x0[i] = ppptr->AxisValues(i);
      nX[i] = x0[i].size();
      nXYZ *= nX[i];
    }

    unsigned int i0 = 0, i1 = 0;
    const realT tfct0 = ppptr->LookupIndices(3, problemTime, i0, i1);
    const realT tfct1 = 1.0 - tfct0;
    i0 *= nXYZ;
    i1 *= nXYZ;
    {
      const std::string meshName("pore_pressure_grid");
      std::map< std::string, array<real64> > realData;

      //-----------------------------
      //Write grid mesh
      //-----------------------------
      {
        dvector xcoords(nXYZ);
        dvector ycoords(nXYZ);
        dvector zcoords(nXYZ);

        realData["Pressure"].resize(nXYZ);
        array<real64>& pp = realData["Pressure"];
        const std::vector<realT>& ppAll = ppptr->Values();

        localIndex a = 0;
        for(localIndex k = 0 ; k < nX[2] ; ++k)
        {
          const realT z = x0[2][k];
          for(localIndex j = 0 ; j < nX[1] ; ++j)
          {
            const realT y = x0[1][j];
            for(localIndex i = 0 ; i < nX[0] ; ++i, ++a)
            {
              const realT x = x0[0][i];
              xcoords[a] = x;
              ycoords[a] = y;
              zcoords[a] = z;
              pp[a] = tfct0 * ppAll[i0 + a] + tfct1 * ppAll[i1 + a];
            }
          }
        }

        realT *coords[3];
        coords[0] = xcoords.data();
        coords[1] = ycoords.data();
        coords[2] = zcoords.data();

        const int ret = siloFile.WriteQuadMeshObject(meshName, nX[0], nX[1], nX[2], coords, cycleNum, problemTime);
        if(ret == -1)
          throw GPException("PhysicalDomainT::WriteFaultElementMesh: Unable to write SILO pore pressure mesh");
      }

      //-----------------------------
      //Write grid data
      //-----------------------------
      {
        std::string subDirectory = "PorePressureGrid";
        std::string rootDirectory = "/" + subDirectory;
        siloFile.MakeSubDirectory( subDirectory, rootDirectory );
        int ret = DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());
        if(ret == -1)
          throw GPException("PhysicalDomainT::WriteFaultElementMesh: Unable to change to PorePressureGrid SILO directory");

        siloFile.WriteFieldMapToSilo<realT>( meshName, realData, DB_NODECENT, cycleNum,
                                             problemTime, false, rootDirectory,
                                             "none", lArray1d() );

        ret = DBSetDir(siloFile.m_dbFilePtr, "..");
        if(ret == -1)
          throw GPException("PhysicalDomainT::WriteFaultElementMesh: Unable to change to parent directory of PorePressureGrid SILO directory");
      }
    }
  }
}

void PhysicalDomainT::ReadFaultElementMesh( const SiloFile& siloFile,
                                            const int cycleNum,
                                            const realT problemTime,
                                            const bool isRestart )
{
  const std::string meshName("fault_mesh");

  const int err = m_faultPatchNodes.ReadSilo( siloFile, "FaultNodeFields", meshName,
                                              DB_NODECENT, cycleNum, problemTime, isRestart );
  if(err)
    return;

  m_faultPatchFaces.ReadSilo( siloFile, "FaultFaceFields", meshName, DB_ZONECENT, cycleNum, problemTime, isRestart );

  m_faultElementManager.ReadSiloFaultElements(siloFile, "FaultElementFields", meshName, DB_ZONECENT, cycleNum, problemTime, isRestart);
}
#endif

void PhysicalDomainT::WriteFiniteElementMesh( SiloFile& siloFile,
                                              const int cycleNum,
                                              const realT problemTime,
                                              const bool isRestart,
                                              const bool writeFEMMesh,
                                              const bool writeFEMFaces,
                                              const bool writeFEMEdges)
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE FE DATA-----------------
  if (writeFEMMesh && m_feElementManager.m_numElems > 0)
  {

    bool writeArbitraryPolygon(false);
    const std::string meshName("volume_mesh");
    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(m_feNodeManager.m_numNodes);
    dvector ycoords(m_feNodeManager.m_numNodes);
    dvector zcoords(m_feNodeManager.m_numNodes);
    for (localIndex a = 0 ; a < m_feNodeManager.m_numNodes ; ++a)
    {
      R1Tensor nodePosition;
      nodePosition = (*m_feNodeManager.m_refposition)[a];
      nodePosition += (*m_feNodeManager.m_displacement)[a];

      xcoords[a] = nodePosition(0);
      ycoords[a] = nodePosition(1);
      zcoords[a] = nodePosition(2);
    }

    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();

    const int numElementRegions = m_feElementManager.m_ElementRegions.size();
    array<localIndex*> meshConnectivity(numElementRegions);
    array<int*> isGhostElement(numElementRegions);
    array<globalIndex*> globalElementNumbers(numElementRegions);
    ivector shapecnt(numElementRegions);
    ivector shapetype(numElementRegions);
    ivector shapesize(numElementRegions);

    array<FixedOneToManyRelation> elementToNodeMap(m_feElementManager.m_ElementRegions.size());

    int count = 0;
    for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator elementRegionIter =
           m_feElementManager.m_ElementRegions.begin() ; elementRegionIter
         != m_feElementManager.m_ElementRegions.end() ; ++elementRegionIter)
    {
      ElementRegionT& elementRegion = elementRegionIter->second;

      // The following line seems to be redundant. It's actual function is to
      // size this temp array.(pfu)
      elementToNodeMap[count] = elementRegion.m_toNodesRelation;

      for (localIndex k = 0 ; k < elementRegion.m_numElems ; ++k)
      {
        const localIndex* const elemToNodeMap =
          elementRegion.m_toNodesRelation[k];

        const array<integer> nodeOrdering = elementRegion.SiloNodeOrdering();
        for (localIndex a = 0 ; a < elementRegion.m_numNodesPerElem ; ++a)
          elementToNodeMap[count](k, a) = elemToNodeMap[nodeOrdering[a]];

      }

      //      meshConnectivity[count] = elementRegion.m_ElementToNodeMap.data();
      meshConnectivity[count] = elementToNodeMap[count].data();

      isGhostElement[count] = (elementRegion.GetFieldData<FieldInfo::ghostRank>()).data();

      globalElementNumbers[count] = elementRegion.m_localToGlobalMap.data();
      shapecnt[count] = elementRegion.m_numElems;

      if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D8") )
      {
        shapetype[count] = DB_ZONETYPE_HEX;
      }
      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D6") )
      {
        shapetype[count] = DB_ZONETYPE_HEX;
        writeArbitraryPolygon = true;
      }
      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D4") )
      {
        shapetype[count] = DB_ZONETYPE_TET;
      }
      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "CPE4") || !elementRegion.m_elementGeometryID.compare(0, 3, "S4R") )
      {
        shapetype[count] = DB_ZONETYPE_QUAD;
      }
      else if ( !elementRegion.m_elementGeometryID.compare(0, 4,
                                                           "STRI") ||
                !elementRegion.m_elementGeometryID.compare(0, 4, "TRSH") || !elementRegion.m_elementGeometryID.compare(0, 4, "CPE3"))
      {
        shapetype[count] = DB_ZONETYPE_TRIANGLE;
      }
      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "CPE2") )
      {
        shapetype[count] = DB_ZONETYPE_TRIANGLE;
      }
      else
      {
        throw GPException("PhysicalDomainT::WriteFiniteElementMesh: Do not recognize geometry type " + elementRegion.m_elementGeometryID + " \n");
      }

      shapesize[count] = elementRegion.m_toNodesRelation.Dimension(1);
      count++;
    }

    siloFile.WriteMeshObject(meshName, m_feNodeManager.m_numNodes, coords,
                             m_feNodeManager.m_localToGlobalMap.data(), numElementRegions,
                             shapecnt.data(), meshConnectivity.data(), globalElementNumbers.data(),
                             isGhostElement.data(), shapetype.data(), shapesize.data(), cycleNum, problemTime);


    // write node fields in silo mesh, and all restart data as unassociated
    // variables.


    m_feNodeManager.WriteSilo(siloFile, "NodalFields", meshName, DB_NODECENT, cycleNum, problemTime, isRestart);



    m_feElementManager.WriteSilo( siloFile, meshName, cycleNum, problemTime, isRestart );



    if ( (isRestart || (writeFEMFaces && m_feFaceManager.DataLengths() > 0)) )
    {

      // face mesh
      const std::string facemeshName("face_mesh");

      if (writeArbitraryPolygon)
      {
        const int numFaceTypes = 1;
        int dbZoneType = DB_ZONETYPE_POLYGON;
        // See a discussion of silo's arbitrary polygon implementation at
        // https://visitbugs.ornl.gov/projects/7/wiki/Arbitrary_Polygons_and_Polyhedra_in_Silo
        // It is not documented in silo manual.
        array<localIndex*> faceConnectivity(numFaceTypes);
        array<globalIndex*> globalFaceNumbers(numFaceTypes);
        ivector fshapecnt(numFaceTypes);
        ivector fshapetype(numFaceTypes);
        ivector fshapesize(numFaceTypes);

        array<lArray1d> faceToNodeMap(numFaceTypes);
        {
          for (localIndex k = 0 ; k < m_feFaceManager.m_numFaces ; ++k)
          {
            faceToNodeMap[0].push_back(m_feFaceManager.m_toNodesRelation[k].size());
            for (localIndex a = 0 ; a < m_feFaceManager.m_toNodesRelation[k].size() ; ++a)
            {
              faceToNodeMap[0].push_back(m_feFaceManager.m_toNodesRelation[k][a]);
            }
          }

          faceConnectivity[0] = faceToNodeMap[0].data();

          globalFaceNumbers[0] = m_feFaceManager.m_localToGlobalMap.data();
          fshapecnt[0] = m_feFaceManager.m_numFaces;
          fshapetype[0] = dbZoneType;
          fshapesize[0] = 0;
        }
        int lnodelist = faceToNodeMap[0].size();

        siloFile.WritePolygonMeshObject(facemeshName, m_feNodeManager.m_numNodes, coords,
                                        m_feNodeManager.m_localToGlobalMap.data(), numFaceTypes,
                                        fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                                        NULL, fshapetype.data(), fshapesize.data(), cycleNum, problemTime, lnodelist);


      }
      else  //The old way
      {
        const int numFaceTypes = 1;
        int numNodesPerFace = m_feFaceManager.m_toNodesRelation[0].size(); // TODO
                                                                           // assumes
                                                                           // all
                                                                           // faces
                                                                           // have
                                                                           // same
                                                                           // number
                                                                           // of
                                                                           // nodes
        int dbZoneType = DB_ZONETYPE_POLYGON;
        if(numNodesPerFace == 3)
        {
          dbZoneType = DB_ZONETYPE_TRIANGLE;
        }
        else if(numNodesPerFace == 4)
        {
          dbZoneType = DB_ZONETYPE_QUAD;
        }
        else if(numNodesPerFace == 2)
        {
          dbZoneType = DB_ZONETYPE_BEAM;
        }

        array<localIndex*> faceConnectivity(numFaceTypes);
        array<globalIndex*> globalFaceNumbers(numFaceTypes);
        ivector fshapecnt(numFaceTypes);
        ivector fshapetype(numFaceTypes);
        ivector fshapesize(numFaceTypes);

        array<lArray2d> faceToNodeMap(numFaceTypes);


        for (int faceType = 0 ; faceType < numFaceTypes ; ++faceType)
        {
          faceToNodeMap[faceType].resize2(m_feFaceManager.m_numFaces, numNodesPerFace);

          for (localIndex k = 0 ; k < m_feFaceManager.m_numFaces ; ++k)
          {
            for (int a = 0 ; a < numNodesPerFace ; ++a)
            {
              faceToNodeMap[faceType][k][a] = m_feFaceManager.m_toNodesRelation[k][a];
            }
          }

          faceConnectivity[faceType] = faceToNodeMap[faceType].data();

          globalFaceNumbers[faceType] = m_feFaceManager.m_localToGlobalMap.data();
          fshapecnt[faceType] = m_feFaceManager.m_numFaces;
          fshapetype[faceType] = dbZoneType;
          fshapesize[faceType] = numNodesPerFace;
        }

        siloFile.WriteMeshObject(facemeshName, m_feNodeManager.m_numNodes, coords,
                                 m_feNodeManager.m_localToGlobalMap.data(), numFaceTypes,
                                 fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                                 NULL, fshapetype.data(), fshapesize.data(), cycleNum, problemTime);
      }


      m_feFaceManager.WriteSilo( siloFile, "FaceFields",facemeshName, DB_ZONECENT, cycleNum, problemTime, isRestart );

#if WRITE_FACENODE
      if (!isRestart)
        m_feNodeManager.WriteSilo(siloFile, "FaceNodeFields", facemeshName, DB_NODECENT, cycleNum, problemTime, isRestart);
#endif

    }

    if ( isRestart || (writeFEMFaces && m_externalFaces.DataLengths() > 0) )
    {

      // external face mesh
      const std::string facemeshName("external_face_mesh");

      if (writeArbitraryPolygon)
      {
        const int numFaceTypes = 1;
        int dbZoneType = DB_ZONETYPE_POLYGON;
        array<localIndex*> faceConnectivity(numFaceTypes);
        array<globalIndex*> globalFaceNumbers(numFaceTypes);
        ivector fshapecnt(numFaceTypes);
        ivector fshapetype(numFaceTypes);
        ivector fshapesize(numFaceTypes);

        array<lArray1d> faceToNodeMap(numFaceTypes);
        {
          for (localIndex k = 0 ; k < m_externalFaces.DataLengths() ; ++k)
          {
            bool is_fe;
            const localIndex faceIndex = m_externalFaces.FaceIndex(k, is_fe);

            faceToNodeMap[0].push_back(m_feFaceManager.m_toNodesRelation[faceIndex].size());
            for (localIndex a = 0 ; a < m_feFaceManager.m_toNodesRelation[faceIndex].size() ; ++a)
            {
              faceToNodeMap[0].push_back(m_feFaceManager.m_toNodesRelation[faceIndex][a]);
            }
          }

          faceConnectivity[0] = faceToNodeMap[0].data();

          globalFaceNumbers[0] = m_feFaceManager.m_localToGlobalMap.data();
          fshapecnt[0] = m_externalFaces.DataLengths();
          fshapetype[0] = dbZoneType;
          fshapesize[0] = 0;
        }
        int lnodelist = faceToNodeMap[0].size();

        siloFile.WritePolygonMeshObject(facemeshName, m_feNodeManager.m_numNodes, coords,
                                        m_feNodeManager.m_localToGlobalMap.data(), numFaceTypes,
                                        fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                                        NULL, fshapetype.data(), fshapesize.data(), cycleNum, problemTime, lnodelist);
      }
      else
      {
        const int numFaceTypes = 1;
        int numNodesPerFace = m_feFaceManager.m_toNodesRelation[0].size(); //TODO
                                                                           // assumes
                                                                           // all
                                                                           // faces
                                                                           // have
                                                                           // same
                                                                           // number
                                                                           // of
                                                                           // nodes
        int dbZoneType = DB_ZONETYPE_POLYGON;
        if(numNodesPerFace == 3)
        {
          dbZoneType = DB_ZONETYPE_TRIANGLE;
        }
        else if(numNodesPerFace == 4)
        {
          dbZoneType = DB_ZONETYPE_QUAD;
        }

        array<localIndex*> faceConnectivity(numFaceTypes);
        array<globalIndex*> globalFaceNumbers(numFaceTypes);
        ivector fshapecnt(numFaceTypes);
        ivector fshapetype(numFaceTypes);
        ivector fshapesize(numFaceTypes);

        array<lArray2d> faceToNodeMap(numFaceTypes);

        for (int faceType = 0 ; faceType < numFaceTypes ; ++faceType)
        {
          faceToNodeMap[faceType].resize2(m_externalFaces.DataLengths(), numNodesPerFace);

          for (localIndex k = 0 ; k < m_externalFaces.DataLengths() ; ++k)
          {
            bool is_fe;
            const localIndex faceIndex = m_externalFaces.FaceIndex(k, is_fe);
            for (int a = 0 ; a < numNodesPerFace ; ++a)
            {
              faceToNodeMap[faceType][k][a] = m_feFaceManager.m_toNodesRelation[faceIndex][a];
            }
          }

          faceConnectivity[faceType] = faceToNodeMap[faceType].data();

          globalFaceNumbers[faceType] = m_externalFaces.m_localToGlobalMap.data();
          fshapecnt[faceType] = m_externalFaces.DataLengths();
          fshapetype[faceType] = dbZoneType;
          fshapesize[faceType] = numNodesPerFace;
        }

        siloFile.WriteMeshObject(facemeshName, m_feNodeManager.m_numNodes, coords,
                                 m_feNodeManager.m_localToGlobalMap.data(), numFaceTypes,
                                 fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                                 NULL, fshapetype.data(), fshapesize.data(), cycleNum, problemTime);
      }

      m_externalFaces.WriteSiloExternalFaces(siloFile, "ExternalFaceFields", facemeshName, DB_ZONECENT, cycleNum, problemTime, isRestart );

    }

    if ( isRestart || (writeFEMEdges && m_feEdgeManager.DataLengths() > 0) )
    {
      // write edges

      const std::string edgeMeshName("edge_mesh");

      const int numEdgeTypes = 1;
      const int numNodesPerEdge = 2;
      int dbZoneType = DB_ZONETYPE_BEAM;

      array<localIndex*> edgeConnectivity(numEdgeTypes);
      array<globalIndex*> globalEdgeNumbers(numEdgeTypes);
      ivector eshapecnt(numEdgeTypes);
      ivector eshapetype(numEdgeTypes);
      ivector eshapesize(numEdgeTypes);

      array<lArray2d> edgeToNodeMap(numEdgeTypes);


      for (int edgeType = 0 ; edgeType < numEdgeTypes ; ++edgeType)
      {
        edgeToNodeMap[edgeType].resize2(m_feEdgeManager.DataLengths(), numNodesPerEdge);

        for (localIndex k = 0 ; k < m_feEdgeManager.DataLengths() ; ++k)
        {
          for (int a = 0 ; a < numNodesPerEdge ; ++a)
          {
            if ( m_feFaceManager.m_toNodesRelation[0].size() == 2 && a > 0)
            {
              edgeToNodeMap[edgeType][k][a] = m_feEdgeManager.m_toNodesRelation[k][0];
            }
            else
            {
              edgeToNodeMap[edgeType][k][a] = m_feEdgeManager.m_toNodesRelation[k][a];
            }
          }
        }

        edgeConnectivity[edgeType] = edgeToNodeMap[edgeType].data();

        globalEdgeNumbers[edgeType] = m_feEdgeManager.m_localToGlobalMap.data();
        eshapecnt[edgeType] = m_feEdgeManager.DataLengths();
        eshapetype[edgeType] = dbZoneType;
        eshapesize[edgeType] = numNodesPerEdge;
      }

      siloFile.WriteMeshObject(edgeMeshName, m_feNodeManager.m_numNodes, coords,
                               m_feNodeManager.m_localToGlobalMap.data(), numEdgeTypes,
                               eshapecnt.data(), edgeConnectivity.data(), globalEdgeNumbers.data(),
                               NULL, eshapetype.data(), eshapesize.data(), cycleNum, problemTime);


      m_feEdgeManager.WriteSilo( siloFile, "EdgeFields",edgeMeshName, DB_ZONECENT, cycleNum, problemTime, isRestart );



//
//
//
//
//      const std::string edgeMeshName("edge_mesh");
//      const localIndex numEdges = m_feEdgeManager.DataLengths();
//
//      array<lArray1d> edgeToNodes(1);
//
//      count = 0;
//
//      for (localIndex k = 0; k <
// m_feEdgeManager.m_toNodesRelation.Dimension(0); ++k)
//      {
//        if ( m_feEdgeManager.m_toNodesRelation.Dimension(1) == 1 )
//        {
//          // This is for 2D only.
//          // Silo polygon does not work if we only have one node for each
// element.
//          edgeToNodes[count].push_back(
// m_feEdgeManager.m_toNodesRelation.Dimension(1)+1 );
//        }
//        else
//        {
//          edgeToNodes[count].push_back(m_feEdgeManager.m_toNodesRelation.Dimension(1));
//        }
//        for (localIndex a = 0; a <
// m_feEdgeManager.m_toNodesRelation.Dimension(1); ++a)
//        {
//          edgeToNodes[count].push_back(m_feEdgeManager.m_toNodesRelation(k,
// a));
//          if ( m_feEdgeManager.m_toNodesRelation.Dimension(1) == 1 )
//          {
//            edgeToNodes[count].push_back(m_feEdgeManager.m_toNodesRelation(k,
// a));
//          }
//        }
//      }
//
//      array<localIndex*> edgeMeshConnectivity(1);
//      ivector eshapecnt(1);
//      ivector eshapetype(1);
//      ivector eshapesize(1);
//
//      edgeMeshConnectivity[count] = edgeToNodes[count].data();
//
//      eshapecnt[count] = numEdges;
//      eshapetype[count] = DB_ZONETYPE_POLYGON;
//      eshapesize[count] = -edgeToNodes[count].size();
//      ++count;
//
//      siloFile.WriteMeshObject(edgeMeshName, m_feNodeManager.m_numNodes,
// coords, NULL, 1,
//                               eshapecnt.data(), edgeMeshConnectivity.data(),
// NULL, NULL,
//                               eshapetype.data(), eshapesize.data(), cycleNum,
// problemTime);
//
//      this->m_feEdgeManager.WriteSilo(siloFile, "EdgeFields",
// edgeMeshName.c_str(), DB_ZONECENT, cycleNum, problemTime, isRestart );

    }

  }//end FE write
}

void PhysicalDomainT::ReadFiniteElementMesh( const SiloFile& siloFile,
                                             const int cycleNum,
                                             const realT problemTime,
                                             const bool isRestart )
{


  int err = m_feNodeManager.ReadSilo( siloFile, "NodalFields", "volume_mesh",
                                      DB_NODECENT, cycleNum, problemTime, isRestart );
//  err = m_feNodeManager.ReadSilo( siloFile, "NodalFieldsB", "face_mesh",
//                                      DB_NODECENT, cycleNum, problemTime,
// isRestart );
  if(err)
    return;

  m_feElementManager.ReadSilo( siloFile, "volume_mesh",
                               cycleNum, problemTime, isRestart );

  m_feNodeManager.ConstructNodeToElementMap( m_feElementManager );

  m_feFaceManager.m_elemManagerHACK = &m_feElementManager;
  err = m_feFaceManager.ReadSilo( siloFile, "FaceFields","face_mesh", DB_ZONECENT, cycleNum, problemTime, isRestart );
  if(err)
    return;

  m_externalFaces.ReadSiloExternalFaces(siloFile, "ExternalFaceFields", "external_face_mesh", DB_ZONECENT, cycleNum, problemTime, isRestart );

  // TODO This is a temporary fix to correct the size of edge to node map.
  //  Later we need to redefine edge to nodes map as
  // OrderedVariableOneToManyRelation
  if (m_feFaceManager.m_toNodesRelation.size() == 0)
  {
    throw GPException("PhysicalDomainT::ReadFiniteElementMesh: No m_feFaceManager.m_toNodesRelation members during FE mesh silo read");
  }
  if (m_feFaceManager.m_toNodesRelation[0].size() == 2)
  {
    m_feEdgeManager.m_toNodesRelation.resize2( 0, 1 );
  }
  m_feEdgeManager.ReadSilo(siloFile, "EdgeFields", "edge_mesh", DB_ZONECENT, cycleNum, problemTime, isRestart );
}



void PhysicalDomainT::WriteDiscreteElements( SiloFile& siloFile,
                                             const int cycleNum,
                                             const realT problemTime,
                                             const bool isRestart,
                                             const bool writeDE )
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE DISCRETE ELEMENT DATA-----------------
  //NOTE: CURRENTLY VISIT DELETES EDGES ON SINGLE ELEMENTS, BUT IF
  //THAT IS FIXED WE CAN SWITCH BACK TO ELEMENT RENDERING ... UNTIL THEN, USE
  // POLYGONS
  //ALSO CHANGE LAGRANGEEXPLICITDYNAMICSSOLVERT.CPP TO _NOT_ DUPLICATE DE FIELDS
  // IN EXTERNAL FACES
  if (writeDE && m_discreteElementManager.DataLengths() > 0)
  {
    const std::string meshName = "de_mesh";

    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(m_discreteElementSurfaceNodes.m_numNodes);
    dvector ycoords(m_discreteElementSurfaceNodes.m_numNodes);
    dvector zcoords(m_discreteElementSurfaceNodes.m_numNodes);
    array<R1Tensor>& cpos =
      m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::currentPosition>();
    for (localIndex a = 0 ; a < m_discreteElementSurfaceNodes.m_numNodes ; ++a)
    {
      xcoords[a] = cpos[a](0);
      ycoords[a] = cpos[a](1);
      zcoords[a] = cpos[a](2);
    }
    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();

    {
      //fill arrays
      //nfaces Total number of unique faces
      //nodecnts Number of nodes per face (size nfaces)
      //sumnodecnts Sum of nodecnts for all i
      //nodelist List of node ids associated with faces (size sumnodects)
      int sumnodecnts = 0, ntmp;
      std::vector<int> nodecnts, nodelist;
      int numfaces = this->m_discreteElementSurfaceFaces.m_toNodesRelation.size();
      for (int i = 0 ; i < numfaces ; i++)
      {
        ntmp = this->m_discreteElementSurfaceFaces.m_toNodesRelation[i].size();
        sumnodecnts += ntmp;
        nodecnts.push_back(ntmp);
        for (int j = 0 ; j < ntmp ; j++)
          nodelist.push_back(this->m_discreteElementSurfaceFaces.m_toNodesRelation[i][j]);
      }

      //facecnts Number of faces per DE (nDiscreteElements)
      //sumfacecnts Sum of facecnts for all i
      //facelist List of face ids associated with DE (size sumfacecnts)
      std::vector<int> facecnts, facelist;
      int sumfacecnts = 0;
      for (localIndex i = 0 ; i < this->m_discreteElementManager.DataLengths() ; i++)
      {
        size_t nfcs = this->m_discreteElementManager.m_discreteElementToExternalFacesMap[i].size();
        sumfacecnts += nfcs;
        facecnts.push_back(nfcs);
        for (localIndex j = 0 ; j < nfcs ; j++)
          facelist.push_back(
            this->m_discreteElementManager.m_discreteElementToExternalFacesMap[i][j]);
      }

      siloFile.WriteDiscreteElementMeshObject(
        meshName.c_str(), m_discreteElementSurfaceNodes.m_numNodes, coords,
        m_discreteElementSurfaceNodes.m_localToGlobalMap.data(),
        m_discreteElementManager.DataLengths(), numfaces, nodecnts.data(), sumnodecnts,
        nodelist.data(), facecnts.data(), sumfacecnts, facelist.data(), 0,   //(int*)this->m_discreteElementSurfaceNodes.m_localToGlobalMap.data(),
        NULL, cycleNum, problemTime);

      if(m_discreteElementManager.m_discreteElementToDiscreteElementContactsMap.size() > 0)
      {
        array<integer> tnodelist;
        std::map<std::string, array<real64> > realFields;
        std::map<std::string, array<R1Tensor> > R1Fields;
        int i = 0;
        for(std::map<localIndex, std::map<localIndex, DiscreteElementContact> >::const_iterator it =
              m_discreteElementManager.m_discreteElementToDiscreteElementContactsMap.begin() ;
            it != m_discreteElementManager.m_discreteElementToDiscreteElementContactsMap.end() ; ++it, ++i)
        {
          for(std::map<localIndex, DiscreteElementContact>::const_iterator it1 = it->second.begin() ; it1 != it->second.end() ; ++it1)
          {
            const int j = static_cast<int>(it1->first);
            tnodelist.push_back(i);
            tnodelist.push_back(j);
            it1->second.Serialize(realFields, R1Fields);
          }
        }
        if(nodelist.size() > 0)
        {
          const std::string contactmeshName = "de_contact_mesh";
          realT* tcoords[3];
          dvector txcoords(m_discreteElementManager.DataLengths());
          dvector tycoords(m_discreteElementManager.DataLengths());
          dvector tzcoords(m_discreteElementManager.DataLengths());
          array<R1Tensor>& tcpos =
            m_discreteElementManager.GetFieldData<FieldInfo::currentPosition>();
          for (localIndex a = 0 ; a < m_discreteElementManager.DataLengths() ; ++a)
          {
            txcoords[a] = tcpos[a](0);
            tycoords[a] = tcpos[a](1);
            tzcoords[a] = tcpos[a](2);
          }
          tcoords[0] = txcoords.data();
          tcoords[1] = tycoords.data();
          tcoords[2] = tzcoords.data();

          siloFile.WriteBeamMesh(contactmeshName, m_discreteElementManager.DataLengths(),
                                 tcoords, tnodelist, cycleNum, problemTime);

          DiscreteElementContact::WriteSilo(siloFile, realFields, R1Fields, "DiscreteElementContactFields",
                                            contactmeshName, DB_ZONECENT, cycleNum, problemTime, isRestart );
        }
      }
    } //end DE write


    if ( isRestart )
    {
      this->m_discreteElementSurfaceNodes.WriteSilo( siloFile, "DiscreteElementNodalFields", meshName.c_str(),
                                                     DB_NODECENT, cycleNum, problemTime, isRestart);

      // face mesh
      const std::string facemeshName("de_face_mesh");

      array<lArray1d> elementToNodeMap(1);
      for (array<lArray1d>::const_iterator it = m_discreteElementManager.m_discreteElementToExternalNodesMap.begin() ;
           it != m_discreteElementManager.m_discreteElementToExternalNodesMap.end() ; ++it)
      {
        elementToNodeMap.begin()->push_back(it->size());
        for(lArray1d::const_iterator it1 = it->begin() ; it1 != it->end() ; ++it1)
          elementToNodeMap.begin()->push_back(*it1);
      }

      array<localIndex*> meshConnectivity(1);
      ivector shapecnt(1);
      ivector shapetype(1);
      ivector shapesize(1);

      meshConnectivity[0] = elementToNodeMap.begin()->data();

      shapecnt[0] = m_discreteElementManager.DataLengths();
      shapetype[0] = DB_ZONETYPE_POLYGON;
      shapesize[0] = -elementToNodeMap.begin()->size();

      siloFile.WriteMeshObject(facemeshName, m_discreteElementSurfaceNodes.DataLengths(), coords, NULL, 1,
                               shapecnt.data(), meshConnectivity.data(), NULL, NULL,
                               shapetype.data(), shapesize.data(), cycleNum, problemTime);

      this->m_discreteElementSurfaceFaces.WriteSilo( siloFile, "DiscreteElementFaceFields", facemeshName.c_str(),
                                                     DB_ZONECENT, cycleNum, problemTime, isRestart );

      this->m_externalFaces.WriteSilo( siloFile, "DiscreteElementExternalFaceFields", facemeshName.c_str(),
                                       DB_ZONECENT, cycleNum, problemTime, isRestart );
    }


    this->m_discreteElementManager.WriteSilo( siloFile, "DiscreteElementFields", meshName.c_str(),
                                              DB_ZONECENT, cycleNum, problemTime, isRestart );
  }
}

void PhysicalDomainT::WriteXFEMElements( SiloFile& siloFile,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const bool writeXFEM )
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE DISCRETE ELEMENT DATA-----------------
  //NOTE: CURRENTLY VISIT DELETES EDGES ON SINGLE ELEMENTS, BUT IF
  //THAT IS FIXED WE CAN SWITCH BACK TO ELEMENT RENDERING ... UNTIL THEN, USE
  // POLYGONS
  //ALSO CHANGE LAGRANGEEXPLICITDYNAMICSSOLVERT.CPP TO _NOT_ DUPLICATE DE FIELDS
  // IN EXTERNAL FACES
  const std::string meshName = "xfem_mesh";

  //set the nodal coordinate data structure

  dvector xcoords;
  dvector ycoords;
  dvector zcoords;

  for (localIndex a = 0 ; a < m_feNodeManager.m_numNodes ; ++a)
  {
    R1Tensor nodePosition = (*m_feNodeManager.m_refposition)[a];
    nodePosition += (*m_feNodeManager.m_displacement)[a];

    xcoords.push_back(nodePosition[0]);
    ycoords.push_back(nodePosition[1]);
    zcoords.push_back(nodePosition(2));
  }

  std::vector<int> nodelist;
  std::vector<int> shapesize;
  std::vector<int> shapecounts;
  std::vector<int> shapetype;
  int nshapetypes = 0;
  int nzones = 0;
  int nnodes = m_feNodeManager.m_numNodes;
  int ndims = 2;

  for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator elementRegionIter =
         m_feElementManager.m_ElementRegions.begin() ; elementRegionIter != m_feElementManager.m_ElementRegions.end() ; ++elementRegionIter) //Loop
                                                                                                                                             // over
                                                                                                                                             // regions
  {
    ElementRegionT& elemRegion = elementRegionIter->second;
    localIndex numEle = elemRegion.DataLengths();

    auto& isPhysical = elemRegion.GetFieldData<int>("isPhysical");
    auto& plotElement = elemRegion.GetFieldData<int>("plotElement");
    const auto& elementToPhysicalNodes = elemRegion.GetVariableOneToManyMap("ElementToPhysicalNodes");

    int numTri = 0, numQuad = 0, numPent = 0;
    /*Connectivity for triangles*/
    for (localIndex iElm=0 ; iElm<numEle ; ++iElm)
    {
      if(isPhysical(iElm) == 1 && plotElement(iElm) == 0)
      {
        if(elementToPhysicalNodes[iElm].size()==3)
        {
          if(numTri==0)
            nshapetypes = nshapetypes+1;
          for (localIndex j = 0 ; j < 3 ; j++)
          {
            if(elementToPhysicalNodes[iElm][j] == 0)
              std::cout<<"This seems like a problem"<<std::endl;
            int nodeIndex = elementToPhysicalNodes[iElm][j];
            nodelist.push_back(nodeIndex);
          }
          numTri = numTri+1;
          nzones = nzones + 1;
        }
      }
    }

    if(numTri!=0)
    {
      shapecounts.push_back(numTri);
      shapesize.push_back(3);
      shapetype.push_back(DB_ZONETYPE_TRIANGLE);
    }

    /*Connectivity for quadrilaterals*/
    for (localIndex iElm=0 ; iElm<numEle ; ++iElm)
    {
      if(isPhysical(iElm) == 1 && plotElement(iElm) == 0)
      {
        if(elementToPhysicalNodes[iElm].size()==4)
        {
          if(numQuad==0)
            nshapetypes = nshapetypes+1;
          for (localIndex j = 0 ; j < 4 ; j++)
          {
            if(elementToPhysicalNodes[iElm][j] == 0)
              std::cout<<"This seems like a problem"<<std::endl;

            int nodeIndex = elementToPhysicalNodes[iElm][j];
            nodelist.push_back(nodeIndex);
          }
          numQuad = numQuad+1;
          nzones = nzones + 1;
        }
      }
    }
    if(numQuad!=0)
    {
      shapecounts.push_back(numQuad);
      shapesize.push_back(4);
      shapetype.push_back(DB_ZONETYPE_QUAD);
    }

    /*Connectivity for pentagons*/
    for (localIndex iElm=0 ; iElm<numEle ; ++iElm)
    {
      if(isPhysical(iElm) == 1 && plotElement(iElm) == 0)
      {
        if(elementToPhysicalNodes[iElm].size()==5)
        {
          if(numPent==0)
            nshapetypes = nshapetypes+1;
          for (localIndex j = 0 ; j < 5 ; j++)
          {
            if(elementToPhysicalNodes[iElm][j] == 0)
              std::cout<<"This seems like a problem"<<std::endl;

            int nodeIndex = elementToPhysicalNodes[iElm][j];
            nodelist.push_back(nodeIndex);
          }
          numPent = numPent+1;
          nzones = nzones + 1;
        }
      }
    }
    if(numPent!=0)
    {
      shapecounts.push_back(numPent);
      shapesize.push_back(5);
      shapetype.push_back(DB_ZONETYPE_POLYGON);
    }
  }

  int lnodelist = nodelist.size();

  siloFile.XFEMMesh(xcoords, ycoords, zcoords, nodelist, lnodelist, shapesize, shapecounts, shapetype, nshapetypes, nnodes, nzones, ndims, cycleNum,
                    problemTime);
  m_feNodeManager.WriteSilo(siloFile, "NodalFieldsXFEM", meshName, DB_NODECENT, cycleNum, problemTime, isRestart);

  /*
     if (1)
     {
     const std::string meshName = "xfem_mesh";

     //set the nodal coordinate data structure
     realT* coords[3];
     dvector xcoords;
     dvector ycoords;
     dvector zcoords;

     auto& isNodePhysical = m_feNodeManager.GetFieldData<int>("isPhysical");
     int numPhyNodes = 0;
     for (localIndex a = 0; a < m_feNodeManager.m_numNodes; ++a)
     {
      if(isNodePhysical(a))
      {
        R1Tensor nodePosition, disp; disp(0) = 0; disp(1) = 1.0; disp(2) = 0;
        nodePosition = (*m_feNodeManager.m_refposition)[a];
        if(a==7)
        {
          nodePosition += disp;
        }
        else if(a==8)
        {
          nodePosition += disp;
        }
        else if(a==9)
        {
          nodePosition += 0;
        }
        else if(a==10)
        {
          nodePosition += 0;
        }
        else if(a==12)
        {
          nodePosition += disp;
        }
        else if(a==13)
        {
          nodePosition += 0;
        }
        else
        {
          nodePosition += (*m_feNodeManager.m_displacement)[a];
        }

        xcoords.push_back(nodePosition(0));
        ycoords.push_back(nodePosition(1));
        zcoords.push_back(nodePosition(2));
        numPhyNodes = numPhyNodes + 1;
      }
     }
     coords[0] = xcoords.data();
     coords[1] = ycoords.data();
     coords[2] = zcoords.data();

     {
      //fill arrays
      //nfaces Total number of unique faces
      //nodecnts Number of nodes per face (size nfaces)
      //sumnodecnts Sum of nodecnts for all i
      //nodelist List of node ids associated with faces (size sumnodects)
      auto& isFacePhysical = m_feFaceManager.GetFieldData<int>("isPhysical");
      int sumnodecnts = 0, ntmp;
      std::vector<int> nodecnts, nodelist;
      int numfaces = 0;
         //this->m_discreteElementSurfaceFaces.m_toNodesRelation.size();
      for (localIndex iFace = 0; iFace < m_feFaceManager.DataLengths(); iFace++)
      {
        if(isFacePhysical(iFace))
        {
          numfaces = numfaces + 1;
          ntmp = m_feFaceManager.m_toNodesRelation[iFace].size();
             //this->m_discreteElementSurfaceFaces.m_toNodesRelation[i].size();
          sumnodecnts += ntmp;
          nodecnts.push_back(ntmp);
          for (int jNod = 0; jNod < ntmp; jNod++)
            nodelist.push_back(this->m_feFaceManager.m_toNodesRelation[iFace][jNod]);
        }
      }

      //facecnts Number of faces per DE (nDiscreteElements)
      //sumfacecnts Sum of facecnts for all i
      //facelist List of face ids associated with DE (size sumfacecnts)
      std::vector<int> facecnts, facelist;
      int sumfacecnts = 0;
      int numXFEMElements = 0;

      for (std::map<ElementManagerT::RegKeyType, ElementRegionT>::iterator
         elementRegionIter =
          m_feElementManager.m_ElementRegions.begin(); elementRegionIter !=
             m_feElementManager.m_ElementRegions.end(); ++elementRegionIter)
             //Loop over regions
      {
        ElementRegionT& elemRegion = elementRegionIter->second;
        localIndex numEle = elemRegion.DataLengths();

        auto& isPhysical = elemRegion.GetFieldData<int>("isPhysical");
        auto& plotElement = elemRegion.GetFieldData<int>("plotElement");
        auto& elementToPhysicalFaces =
           elemRegion.GetVariableOneToManyMap("ElementToPhysicalFaces");

        for (localIndex iElm=0 ; iElm<numEle ; ++iElm)
        {
          if(isPhysical(iElm) == 1 && plotElement(iElm) == 0)
          {
            //        size_t nfcs =
               this->m_discreteElementManager.m_discreteElementToExternalFacesMap[i].size();
            size_t nfcs = elementToPhysicalFaces[iElm].size();
            sumfacecnts += nfcs;
            facecnts.push_back(nfcs);
            for (localIndex j = 0; j < nfcs; j++)
            {
              facelist.push_back(elementToPhysicalFaces[iElm][j]);
            }
            numXFEMElements = numXFEMElements+1;
          }
        }
      }

      siloFile.WriteDiscreteElementMeshObject(
          meshName.c_str(), numPhyNodes, coords,
          0,
          numXFEMElements, numfaces, nodecnts.data(), sumnodecnts,
          nodelist.data(), facecnts.data(), sumfacecnts, facelist.data(), 0,
             //(int*)this->m_discreteElementSurfaceNodes.m_localToGlobalMap.data(),
          NULL, cycleNum, problemTime);

     //      siloFile.TestPolyhedralCells();

     } //end DE write
     }
   */
}

void PhysicalDomainT::ReadDiscreteElements( const SiloFile& siloFile,
                                            const int cycleNum,
                                            const realT problemTime,
                                            const bool isRestart )
{
  const std::string meshName = "de_mesh";
  const std::string facemeshName = "de_face_mesh";

  int err = this->m_discreteElementSurfaceNodes.ReadSilo( siloFile, "DiscreteElementNodalFields", meshName.c_str(),
                                                          DB_NODECENT, cycleNum, problemTime, isRestart);
  if(err)
    return;

  err = this->m_discreteElementSurfaceFaces.ReadSilo( siloFile, "DiscreteElementFaceFields", facemeshName.c_str(),
                                                      DB_ZONECENT, cycleNum, problemTime, isRestart );
  if(err)
    return;

  err = this->m_discreteElementManager.ReadSilo( siloFile, "DiscreteElementFields", meshName.c_str(),
                                                 DB_ZONECENT, cycleNum, problemTime, isRestart );
  if(err)
    return;

  if(this->m_externalFaces.DataLengths() != 0 && this->m_discreteElementManager.DataLengths() != 0)
    throw GPException("We do not currently support reading a restart file for mixed FE-DE interface problems");

  err = this->m_externalFaces.ReadSilo( siloFile, "DiscreteElementExternalFaceFields", facemeshName.c_str(),
                                        DB_ZONECENT, cycleNum, problemTime, isRestart );
  if(err)
    return;

  this->m_externalFaces.BuildExternalFacesDiscreteElements();
}

void PhysicalDomainT::WriteEllipsoidalDiscreteElements( SiloFile& siloFile,
                                                        const int cycleNum,
                                                        const realT problemTime,
                                                        const bool isRestart,
                                                        const bool writeEDE )
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE ELLIPSOIDAL DISCRETE ELEMENT DATA-----------------

  if(writeEDE && this->m_ellipsoidalDiscreteElementManager.DataLengths() > 0)
  {
    const std::string meshName = "ede_mesh";
    const localIndex numPoints = this->m_ellipsoidalDiscreteElementManager.DataLengths();

    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(numPoints);
    dvector ycoords(numPoints);
    dvector zcoords(numPoints);
    array<R1Tensor>& cpos = m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::currentPosition> ();
    for (localIndex a = 0 ; a < numPoints ; ++a)
    {
      xcoords[a] = cpos[a](0);
      ycoords[a] = cpos[a](1);
      zcoords[a] = cpos[a](2);
    }
    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();
    siloFile.WritePointMesh( meshName, numPoints, coords, cycleNum, problemTime);



    const array<real64>& q = m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::rotationMagnitude>();
    const array<R1Tensor>& xyz = m_ellipsoidalDiscreteElementManager.GetFieldData<FieldInfo::rotationAxis>();
    const array<R1Tensor>& r = m_ellipsoidalDiscreteElementManager.GetFieldData<R1Tensor>("principalRadii");

    array<R2Tensor> R(m_ellipsoidalDiscreteElementManager.DataLengths());
    array<array<real64> > RComponents(9);


    for( int a=0 ; a<9 ; ++a )
    {
      RComponents[a].resize( m_ellipsoidalDiscreteElementManager.DataLengths() );
    }


    for (localIndex a = 0 ; a < m_ellipsoidalDiscreteElementManager.DataLengths() ; ++a)
    {
      DiscreteElementManagerBaseT::QuaternionToRotation(q[a], xyz[a], R[a]);


//      std::cout<<R[a]<<std::endl;

      R2Tensor PRad;
      PRad(0,0) = 2*r[a][0];
      PRad(1,1) = 2*r[a][1];
      PRad(2,2) = 2*r[a][2];

      R2Tensor temp0,temp;

      temp0.AijBjk( R[a], PRad );
      temp.AijBkj( temp0, R[a] );

      RComponents[0][a] = temp(0,0);
      RComponents[1][a] = temp(0,1);
      RComponents[2][a] = temp(0,2);
      RComponents[3][a] = temp(1,0);
      RComponents[4][a] = temp(1,1);
      RComponents[5][a] = temp(1,2);
      RComponents[6][a] = temp(2,0);
      RComponents[7][a] = temp(2,1);
      RComponents[8][a] = temp(2,2);
    }

    m_ellipsoidalDiscreteElementManager.WriteSilo( siloFile, "EllipsoidDiscreteElementFields",
                                                   meshName, DB_ZONECENT, cycleNum,
                                                   problemTime, isRestart);

    std::string subDirectory =   "EllipsoidDiscreteElementFields";
    std::string rootDirectory = "/EllipsoidDiscreteElementFields";
    DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

    array<string> RNames(9);
    RNames[0] = "R_11";
    RNames[1] = "R_12";
    RNames[2] = "R_13";
    RNames[3] = "R_21";
    RNames[4] = "R_22";
    RNames[5] = "R_23";
    RNames[6] = "R_31";
    RNames[7] = "R_32";
    RNames[8] = "R_33";

    for( int a=0 ; a<9 ; ++a )
    {
      siloFile.WriteDataField<realT>( meshName, RNames[a], RComponents[a], 0, cycleNum, problemTime, rootDirectory, "none" );
    }


    DBSetDir(siloFile.m_dbFilePtr, "..");


    {

      char pwd[256];
      DBGetDir(siloFile.m_dbFilePtr, pwd);
      DBSetDir(siloFile.m_dbFilePtr,"/");
      const char name[] = "DefVars";
      const int ndefs = 1;

      array<string> vnames;
      vnames.push_back("RotationTensor");

      array<string> vdefns;
      vdefns.push_back("{{</EllipsoidDiscreteElementFields/R_11>,</EllipsoidDiscreteElementFields/R_12>,</EllipsoidDiscreteElementFields/R_13>}, "
                       "{</EllipsoidDiscreteElementFields/R_21>,</EllipsoidDiscreteElementFields/R_22>,</EllipsoidDiscreteElementFields/R_23>},"
                       "{</EllipsoidDiscreteElementFields/R_31>,</EllipsoidDiscreteElementFields/R_32>,</EllipsoidDiscreteElementFields/R_33>}}");

      int types[ndefs] = { DB_VARTYPE_TENSOR };

      const char *names[ndefs];
      const char *defns[ndefs];
      DBoptlist *optlist[ndefs];

      for( int i=0 ; i<ndefs ; ++i )
      {
        names[i] = vnames[i].c_str();
        defns[i] = vdefns[i].c_str();

        optlist[i] = DBMakeOptlist(2);
        DBAddOption(optlist[i], DBOPT_CYCLE, const_cast<int*> (&cycleNum));
        DBAddOption(optlist[i], DBOPT_DTIME, const_cast<realT*> (&problemTime));
      }


      DBPutDefvars( siloFile.m_dbFilePtr, name, ndefs,
                    const_cast<char**>(names), types,
                    const_cast<char**>(defns), optlist );

      DBSetDir(siloFile.m_dbFilePtr,pwd);
    }

    //m_ellipsoidalDiscreteElementManager.WriteVTK(cycleNum,
    // m_ellipsoidalContactManager);

    if(m_ellipsoidalContactManager.DataLengths() > 0)
    {
      const std::string contactmeshName = "ede_contact_mesh";
      siloFile.WriteBeamMesh(contactmeshName, numPoints, coords,
                             m_ellipsoidalContactManager.GetFieldData<localIndex>("face1"),
                             m_ellipsoidalContactManager.GetFieldData<localIndex>("face2"),
                             cycleNum, problemTime);
      m_ellipsoidalContactManager.WriteSilo(siloFile, "EllipsoidContactFields",
                                            contactmeshName, DB_ZONECENT, cycleNum,
                                            problemTime, isRestart);
    }
  }
}

void PhysicalDomainT::ReadEllipsoidalDiscreteElements( const SiloFile& siloFile,
                                                       const int cycleNum,
                                                       const realT problemTime,
                                                       const bool isRestart )
{
  const std::string meshName = "ede_mesh";
  m_ellipsoidalDiscreteElementManager.ReadSilo( siloFile, "EllipsoidDiscreteElementFields",
                                                meshName, DB_ZONECENT, cycleNum,
                                                problemTime, isRestart);

  const std::string contactmeshName = "ede_contact_mesh";
  m_ellipsoidalContactManager.ReadSilo(siloFile, "EllipsoidContactFields",
                                       contactmeshName, DB_ZONECENT, cycleNum,
                                       problemTime, isRestart);

}


void PhysicalDomainT::WriteCommonPlanes( SiloFile& siloFile,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const bool writeCP )
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE COMMON PLANE DATA-----------------

  const localIndex numContacts = m_contactManager.m_contactToIntersectionPolygonPointsMap.size();

  //  std::cout << "number of contacts: " << numContacts;
  if (writeCP && m_contactManager.m_contact)
  {
    array<lArray1d> elementToNodeMap(1);
    array<localIndex> mask;

    const int count = 0;
    int numCommonPlanes = 0;

    for (localIndex k = 0 ; k < numContacts ; ++k)
    {
      const lArray1d::size_type numPointsInCP = m_contactManager.m_contactToIntersectionPolygonPointsMap[k].size();
      if (numPointsInCP > 0)
      {
        mask.push_back(k);
        elementToNodeMap[count].push_back(numPointsInCP);
        for (lArray1d::size_type a = 0 ; a < numPointsInCP ; ++a)
        {
          elementToNodeMap[count].push_back( m_contactManager.m_contactToIntersectionPolygonPointsMap[k][a]);
        }
        ++numCommonPlanes;
      }
    }

    // FOR DEFAULT VISUALIZATION ... ELIMINATE ANNOYING VISIT MESSAGES
    bool defaultPlane = false;
    if(numCommonPlanes==0)
    {
      defaultPlane = true;
      if(isRestart)
        return;
      mask.push_back(m_contactManager.DataLengths());
      elementToNodeMap[count].push_back(3);
      for (lArray1d::size_type a = 0 ; a < 3 ; ++a)
      {
        elementToNodeMap[count].push_back(m_contactManager.m_intersectionPolygonPoints.size());
        m_contactManager.m_intersectionPolygonPoints.push_back(m_contactManager.DefaultPolygonCenter());
        m_contactManager.m_intersectionPolygonPoints.back()(a) += m_contactManager.DefaultPolygonDimension();
      }
      ++numCommonPlanes;
      m_contactManager.resize(m_contactManager.DataLengths() + 1); //see line
                                                                   // 1639 for
                                                                   // the
                                                                   // reversal
                                                                   // of this
                                                                   // operation
    }

    // FOR DEFAULT VISUALIZATION : NUMBER OF COMMON PLANES ALWAYS AT LEAST 1
    {
      const localIndex numCommonPlanePoints = m_contactManager.m_intersectionPolygonPoints.size();
      const std::string commonPlaneMeshName("cp_mesh");

      //set the nodal coordinate data structure
      realT* coords[3];
      dvector xcoords(numCommonPlanePoints);
      dvector ycoords(numCommonPlanePoints);
      dvector zcoords(numCommonPlanePoints);
      array<R1Tensor>& cpos = m_contactManager.m_intersectionPolygonPoints;
      for (localIndex a = 0 ; a < numCommonPlanePoints ; ++a)
      {
        xcoords[a] = cpos[a](0);
        ycoords[a] = cpos[a](1);
        zcoords[a] = cpos[a](2);
      }
      coords[0] = xcoords.data();
      coords[1] = ycoords.data();
      coords[2] = zcoords.data();

      array<localIndex*> meshConnectivity(1);
      ivector shapecnt(1);
      ivector shapetype(1);
      ivector shapesize(1);

      meshConnectivity[count] = elementToNodeMap[count].data();

      shapecnt[count] = numCommonPlanes;
      shapetype[count] = DB_ZONETYPE_POLYGON;
      shapesize[count] = -elementToNodeMap[count].size();
      //++count;

      siloFile.WriteMeshObject(commonPlaneMeshName, numCommonPlanePoints, coords, NULL, 1,
                               shapecnt.data(), meshConnectivity.data(), NULL, NULL,
                               shapetype.data(), shapesize.data(), cycleNum, problemTime);


      this->m_contactManager.WriteSilo( siloFile, "CommonPlaneFields", commonPlaneMeshName.c_str(),
                                        DB_ZONECENT, cycleNum,
                                        problemTime, isRestart, "none", mask );

    }

    // FOR DEFAULT VISUALIZATION : UNDO THE DEFAULT VISUALIZATION
    if(defaultPlane)
    {
      for (lArray1d::size_type a = 0 ; a < 3 ; ++a)
        m_contactManager.m_intersectionPolygonPoints.pop_back();
      m_contactManager.resize(((int)m_contactManager.DataLengths()) + 1); //see
                                                                          // line
                                                                          // 1587
                                                                          // for
                                                                          // the
                                                                          // initiation
                                                                          // of
                                                                          // this
                                                                          // operation
    }
  }
}

void PhysicalDomainT::ReadCommonPlanes( const SiloFile& siloFile,
                                        const int cycleNum,
                                        const realT problemTime,
                                        const bool isRestart )
{
  const std::string commonPlaneMeshName("cp_mesh");
  this->m_contactManager.ReadSilo( siloFile, "CommonPlaneFields", commonPlaneMeshName.c_str(),
                                   DB_ZONECENT, cycleNum,
                                   problemTime, isRestart, "none" );
}

void PhysicalDomainT::WriteCartesianGrid( SiloFile& siloFile,
                                          const int cycleNum,
                                          const realT problemTime,
                                          const bool isRestart,
                                          const bool writeCG )
{
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE CARTESIAN GRID DATA-----------------
  if (writeCG && m_cartesianGridManager.DataLengths() > 0)
  {
    const std::string meshName("cartesian_grid");
    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(m_cartesianGridManager.m_nXYZ);
    dvector ycoords(m_cartesianGridManager.m_nXYZ);
    dvector zcoords(m_cartesianGridManager.m_nXYZ);

    localIndex nX = m_cartesianGridManager.m_nX;
    localIndex nY = m_cartesianGridManager.m_nY;
    localIndex nZ = m_cartesianGridManager.m_nZ;

    for (localIndex a = 0 ; a < m_cartesianGridManager.m_nXYZ ; ++a)
    {
      R1Tensor& refPosition = (*m_cartesianGridManager.m_refposition)[a];

      xcoords[a] = refPosition(0);
      ycoords[a] = refPosition(1);
      zcoords[a] = refPosition(2);
    }

    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();


    siloFile.WriteQuadMeshObject(meshName,nX,nY,nZ,coords,cycleNum, problemTime);


    // write data

    this->m_cartesianGridManager.WriteSilo( siloFile, "CartesianGridFields", meshName.c_str(), DB_NODECENT, cycleNum,
                                            problemTime, isRestart);

  }
}

void PhysicalDomainT::ReadCartesianGrid( const SiloFile& siloFile,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart )
{
  const std::string meshName("cartesian_grid");

  this->m_cartesianGridManager.ReadSilo( siloFile, "CartesianGridFields", meshName.c_str(), DB_NODECENT, cycleNum,
                                         problemTime, isRestart);

}

//////////////////////////////////////////////////////////////////

void PhysicalDomainT::RegisterBCFields(){
  RegisterBCFields(m_ellipsoidalDiscreteElementManager);
  RegisterBCFields(m_ellipsoidalContactManager);
  RegisterBCFields(m_discreteElementSurfaceFaces);
  RegisterBCFields(m_discreteElementSurfaceNodes);
  RegisterBCFields(m_discreteElementManager);
  RegisterBCFields(m_contactManager);
  RegisterBCFields(m_feNodeManager);
  RegisterBCFields(m_feElementManager);
  RegisterBCFields(m_feFaceManager);
  RegisterBCFields(m_feEdgeManager);
  RegisterBCFields(m_externalFaces);
#ifdef SRC_EXTERNAL
  RegisterBCFields(m_faultPatchNodes);
  RegisterBCFields(m_faultPatchFaces);
#endif
}

void PhysicalDomainT::RegisterBCFields(ObjectDataStructureBaseT& object){

  for (std::vector<BoundaryConditionBase*>::iterator it_ic = object.m_bcData.begin() ;
       it_ic != object.m_bcData.end() ; ++it_ic)
  {
    (*it_ic)->RegisterFields(*this);
  }

}


void PhysicalDomainT::UpdateEnergy()
{
  const array<real64>& nodalMass = m_feNodeManager.GetFieldData<FieldInfo::mass>();
  const array<R1Tensor>& velocity = m_feNodeManager.GetFieldData<FieldInfo::velocity>();
  const array<integer>& ghostRank = m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

  const array<real64>& work = m_feNodeManager.GetFieldData<realT>("work");

  this->m_energy.Zero();

  realT KE = 0;
  realT workDoneNodes = 0;
  for( localIndex a=0 ; a<m_feNodeManager.m_numNodes ; ++a )
  {
    if( ghostRank[a]<0 )
    {
      KE += 0.5 * nodalMass[a] * Dot(velocity[a],velocity[a]);
      workDoneNodes += work[a];
    }
  }
  this->m_energy.SetKineticEnergy(KE);
  this->m_energy.SetWorkDoneOnNodes(workDoneNodes);



  for( std::map< std::string, ElementRegionT >::const_iterator iter = m_feElementManager.m_ElementRegions.begin() ;
       iter!=m_feElementManager.m_ElementRegions.end() ; ++iter )
  {
    const ElementRegionT& elemRegion = iter->second;

    this->m_energy += elemRegion.m_energy;
  }
}
