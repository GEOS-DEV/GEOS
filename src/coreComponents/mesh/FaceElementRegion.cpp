/*
 * FaceElementRegion.cpp
 *
 *  Created on: May 15, 2019
 *      Author: settgast
 */

#include "FaceElementRegion.hpp"

namespace geosx
{
using namespace dataRepository;

FaceElementRegion::FaceElementRegion( string const & name, ManagedGroup * const parent ):
  ElementRegion( name, parent ),
  m_fractureSetNames(),
  m_edgesToFractureConnectors(),
  m_fractureConnectorsToEdges(),
  m_fractureElementConnectors(),
  m_fractureCellConnectorIndices(),
  m_fractureToCellConnectors()
{
  this->GetGroup(viewKeyStruct::elementSubRegions)->RegisterGroup<FaceElementSubRegion>("default");

  RegisterViewWrapper( viewKeyStruct::fractureSetString, &m_fractureSetNames, false )->
    setInputFlag(InputFlags::OPTIONAL);

  RegisterViewWrapper( viewKeyStruct::edgesTofractureConnectorsMapString, &m_edgesToFractureConnectors, 0 )
    ->setRestartFlags( RestartFlags::NO_WRITE)
    ->setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::fractureConnectorToEdgeMapString, &m_fractureConnectorsToEdges, 0 )
    ->setRestartFlags( RestartFlags::NO_WRITE)
    ->setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::fractureElementConnectorString, &m_fractureElementConnectors, 0 )
    ->setRestartFlags( RestartFlags::NO_WRITE)
    ->setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::fractureCellConnectorIndicesString, &m_fractureCellConnectorIndices, 0 )
    ->setRestartFlags( RestartFlags::NO_WRITE)
    ->setSizedFromParent(0);

  RegisterViewWrapper( viewKeyStruct::fractureToCellConnectorString, &m_fractureToCellConnectors, 0 )
    ->setRestartFlags( RestartFlags::NO_WRITE)
    ->setSizedFromParent(0);


}

FaceElementRegion::~FaceElementRegion()
{
  // TODO Auto-generated destructor stub
}



localIndex FaceElementRegion::AddToFractureMesh( FaceManager const * const faceManager,
                                                 string const & subRegionName,
                                                 localIndex const faceIndices[2]  )
{
  localIndex rval = -1;

  set<localIndex> connectedEdges;

  array2d<localIndex > const & faceToElementRegion = faceManager->elementRegionList();
  array2d<localIndex > const & faceToElementSubRegion = faceManager->elementSubRegionList();
  array2d<localIndex > const & faceToElementIndex = faceManager->elementList();

  ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);

  FaceElementSubRegion * subRegion = elementSubRegions->GetGroup<FaceElementSubRegion>(subRegionName);
  subRegion->resize( subRegion->size() + 1 );
  rval = subRegion->size() - 1;

  m_fractureToCellConnectors.resize( subRegion->size(), 2 );

  FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
  FaceElementSubRegion::EdgeMapType & edgeMap = subRegion->edgeList();
  FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

  OrderedVariableOneToManyRelation const & facesToNodesMap = faceManager->nodeList();
  OrderedVariableOneToManyRelation const & facesToEdgesMap = faceManager->edgeList();

  localIndex const kfe = subRegion->size() - 1;

  faceMap[kfe][0] = faceIndices[0];
  faceMap[kfe][1] = faceIndices[1];

  arrayView1d<localIndex const> const & faceToNodesMap0 = facesToNodesMap[faceIndices[0]];
  arrayView1d<localIndex const> const & faceToNodesMap1 = facesToNodesMap[faceIndices[1]];
  nodeMap[kfe].resize( faceToNodesMap0.size() * 2 );
  for( localIndex a=0 ; a<faceToNodesMap0.size() ; ++a )
  {
    const localIndex aa = a == 0 ? a : faceToNodesMap0.size() - a;

    // TODO HACK need to generalize to something other than quads
    nodeMap[kfe][a]   = faceToNodesMap0[a];
    nodeMap[kfe][a+4] = faceToNodesMap1[aa];
  }

  arrayView1d<localIndex const> const & faceToEdgesMap = facesToEdgesMap[faceIndices[0]];
  edgeMap[kfe].resize( faceToEdgesMap.size() );
  for( localIndex a=0 ; a<faceToEdgesMap.size() ; ++a )
  {
    edgeMap[kfe][a] = faceToEdgesMap[a];
    connectedEdges.insert( faceToEdgesMap[a] );
  }

  for( localIndex ke=0 ; ke<2 ; ++ke )
  {
    m_fractureToCellConnectors.m_toElementRegion[kfe][ke] = faceToElementRegion[faceIndices[ke]][0];
    m_fractureToCellConnectors.m_toElementSubRegion[kfe][ke] = faceToElementSubRegion[faceIndices[ke]][0];
    m_fractureToCellConnectors.m_toElementIndex[kfe][ke] = faceToElementIndex[faceIndices[ke]][0];
  }


  for( auto const & edge : connectedEdges )
  {
    if( m_edgesToFractureConnectors[edge] > 0 )
    {
      localIndex const connectorIndex = m_edgesToFractureConnectors[edge];
      m_fractureElementConnectors[connectorIndex].resize( m_fractureElementConnectors[connectorIndex].size() + 1 );
      m_fractureElementConnectors[connectorIndex][m_fractureElementConnectors[connectorIndex].size()] = kfe;
    }
  }


  for( auto const & setIter : faceManager->sets()->wrappers() )
  {
    set<localIndex> const & faceSet = faceManager->sets()->getReference<set<localIndex> >( setIter.first );
    set<localIndex> & faceElementSet = subRegion->sets()->RegisterViewWrapper< set<localIndex> >( setIter.first )->reference();
    for( localIndex a=0 ; a<faceMap.size(0) ; ++a )
    {
      localIndex const faceIndex = faceMap[a][0];
      if( faceSet.count( faceIndex ) )
      {
        faceElementSet.insert( a );
      }
    }
  }

  return rval;
}

 void FaceElementRegion::GenerateFractureMesh( FaceManager const * const faceManager )
 {

   if( this->m_fractureSetNames.empty() )
   {
     return;
   }

   // key is edge index, value is faceElementIndex....this only works for a single fracture Region with a single subregion!!
   map< localIndex, set<localIndex> > edgeToFractureElementMap;

   array2d<localIndex > const & faceToElementRegion = faceManager->elementRegionList();
   array2d<localIndex > const & faceToElementSubRegion = faceManager->elementSubRegionList();
   array2d<localIndex > const & faceToElementIndex = faceManager->elementList();

   ManagedGroup * elementSubRegions = this->GetGroup(viewKeyStruct::elementSubRegions);
   for( string const & setName : this->m_fractureSetNames )
   {
     FaceElementSubRegion * const subRegion = elementSubRegions->RegisterGroup<FaceElementSubRegion>("default");
     set<localIndex> const & targetSet = faceManager->sets()->getReference<set<localIndex> >(setName);
     subRegion->resize( targetSet.size() );
     subRegion->numNodesPerElement() = 8;

     m_fractureToCellConnectors.resize( targetSet.size(), 2 );

     FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
     FaceElementSubRegion::EdgeMapType & edgeMap = subRegion->edgeList();
     FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

     OrderedVariableOneToManyRelation const & facesToNodesMap = faceManager->nodeList();
     OrderedVariableOneToManyRelation const & facesToEdgesMap = faceManager->edgeList();

     localIndex kfe = 0;
     for( auto const faceIndex : targetSet )
     {
       faceMap[kfe][0] = faceIndex;
       faceMap[kfe][1] = faceIndex;

       arrayView1d<localIndex const> const & faceToNodesMap = facesToNodesMap[faceIndex];
       nodeMap[kfe].resize( faceToNodesMap.size() * 2 );
       for( localIndex a=0 ; a<faceToNodesMap.size() ; ++a )
       {
         const localIndex aa = a == 0 ? a : faceToNodesMap.size() - a;

         // TODO HACK need to generalize to something other than quads
         nodeMap[kfe][a] = faceToNodesMap[a];
         nodeMap[kfe][a+4] = faceToNodesMap[aa];
       }

       arrayView1d<localIndex const> const & faceToEdgesMap = facesToEdgesMap[faceIndex];
       edgeMap[kfe].resize( faceToEdgesMap.size() );
       for( localIndex a=0 ; a<faceToEdgesMap.size() ; ++a )
       {
         edgeMap[kfe][a] = faceToEdgesMap[a];
         edgeToFractureElementMap[ faceToEdgesMap[a] ].insert( kfe );
       }

       for( localIndex ke=0 ; ke<2 ; ++ke )
       {
         m_fractureToCellConnectors.m_toElementRegion[kfe][ke] = faceToElementRegion[faceIndex][ke];
         m_fractureToCellConnectors.m_toElementSubRegion[kfe][ke] = faceToElementSubRegion[faceIndex][ke];
         m_fractureToCellConnectors.m_toElementIndex[kfe][ke] = faceToElementIndex[faceIndex][ke];
       }
       ++kfe;
     }
   }

   m_fractureConnectorsToEdges.resize( edgeToFractureElementMap.size() );
   m_fractureElementConnectors.resize( edgeToFractureElementMap.size() );
   localIndex connectorIndex=0;
   for( auto const & connector : edgeToFractureElementMap )
   {
     if( connector.second.size() > 1 )
     {
       m_edgesToFractureConnectors[connector.first] = connectorIndex;
       m_fractureConnectorsToEdges[connectorIndex] = connector.first;
       m_fractureElementConnectors[connectorIndex].resize( connector.second.size() );
       localIndex fractureElementCounter = -1;
       for( auto const fractureElementIndex : connector.second )
       {
         m_fractureElementConnectors[connectorIndex][++fractureElementCounter] = fractureElementIndex;
       }
       ++connectorIndex;
     }
   }
   m_fractureConnectorsToEdges.resize(connectorIndex);
   m_fractureElementConnectors.resize(connectorIndex);


   forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion  * const subRegion )
   {
     FaceElementSubRegion::FaceMapType const & faceMap = subRegion->faceList();
     for( auto const & setIter : faceManager->sets()->wrappers() )
     {
       set<localIndex> const & faceSet = faceManager->sets()->getReference<set<localIndex> >( setIter.first );
       set<localIndex> & faceElementSet = subRegion->sets()->RegisterViewWrapper< set<localIndex> >( setIter.first )->reference();
       for( localIndex a=0 ; a<faceMap.size(0) ; ++a )
       {
         localIndex const faceIndex = faceMap[a][0];
         if( faceSet.count( faceIndex ) )
         {
           faceElementSet.insert( a );
         }
       }
     }
   });

 }

REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceElementRegion, std::string const &, ManagedGroup * const )

} /* namespace geosx */
