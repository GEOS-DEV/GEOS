
/**
 * @file FractunatorBase.cpp
 * @author settgast1
 * @date Jul 14, 2011
 */

#include "FractunatorBase.h"
#include <limits.h>
#include "IO/BinStream.h"

#include "MPI_Communications/SpatialPartition.h"

#include "Constitutive/CohesiveZone/CohesiveZoneFactory.h"



FractunatorBase::FractunatorBase():
  m_verbose(0),
  m_failstress(0),
  m_checkInterval(1),
  m_rockToughness(1.0e99),
  m_maxTurnAngle(91.0/180.0 * 3.14159265)
{
  // TODO Auto-generated constructor stub


//  m_virtualNodes.AddMap< array<lSet> >("usedFaces");
}

FractunatorBase::~FractunatorBase()
{
  // TODO Auto-generated destructor stub
}


void FractunatorBase::RegisterFieldsAndMaps( NodeManager& nodeManager,
                                             EdgeManagerT& edgeManager,
                                             FaceManagerT& faceManager )
{


  // the virtual FaceManager's rutpureState will be used with the following
  // definitions:
  //   ruptureState = 0 means not reached rupture criteria
  //                = 1 means has reached rupture criteria
  //                = 2 means that face has already been split...but it can
  // still
  //                  be part of a new separation path!!!
  faceManager.AddKeylessDataField<int>( "ruptureState",true, true );

  nodeManager.AddKeylessDataField<int>("numberOfRupturedFaces",true,true);


  nodeManager.AddMap< UnorderedVariableOneToManyRelation >("usedFaces");
  nodeManager.AddMap<OrderedVariableOneToManyRelation>("childIndices");
  nodeManager.AddMap<OneToOneRelation>("parentIndex");
  OneToOneRelation& parentIndexNodes = nodeManager.GetOneToOneMap( "parentIndex" );
  parentIndexNodes = LOCALINDEX_MAX;

  edgeManager.AddMap<OrderedVariableOneToManyRelation>("childIndices");
  edgeManager.AddMap<OneToOneRelation>("parentIndex");
  OneToOneRelation& parentIndexEdge = edgeManager.GetOneToOneMap( "parentIndex" );
  parentIndexEdge = LOCALINDEX_MAX;

  faceManager.AddMap<OrderedVariableOneToManyRelation>("childIndices");
  faceManager.AddMap<OneToOneRelation>("parentIndex");
  OneToOneRelation& parentIndexFace = faceManager.GetOneToOneMap( "parentIndex" );
  parentIndexFace = LOCALINDEX_MAX;

  faceManager.AddKeylessDataField<R1Tensor>( "gapVector",true, true );
  faceManager.AddKeylessDataField<realT>( "separationCoeff",true, true );
  nodeManager.AddKeylessDataField<R1Tensor>("cohesiveForce", true, true);

  faceManager.AddKeylessDataField<realT>( "birthTime",true, true );

  faceManager.AddKeylessDataField<int>("isSeparable", true, true );

  array<real64>* toughnessPointer = faceManager.GetFieldDataPointer<realT>("faceToughness");
  if (toughnessPointer == NULL)
  {
    faceManager.AddKeylessDataField<realT>("faceToughness", true, true);
    m_tounessSetByInitialCondition = false;
  }
  else
  {
    m_tounessSetByInitialCondition = true;
  }

  array<real64>* failStressPointer = faceManager.GetFieldDataPointer<realT>("faceFailStress");
  if (failStressPointer == NULL)
  {
    if (m_failCriterion == 0 || m_failCriterion == 2)
    {
      faceManager.AddKeylessDataField<realT>("faceFailStress", true, true);
    }
    m_failStressSetByInitialCondition = false;
  }
  else
  {
    m_failStressSetByInitialCondition = true;
  }

  if (!m_dfnPrefix.empty())
  {
    faceManager.AddKeylessDataField<int>( "DFN_Index", true, true );
  }

  array<real64>* delta0N = faceManager.GetFieldDataPointer<realT>("delta0N");
  if (delta0N != NULL)
  {
    faceManager.AddKeylessDataField<realT>( "faceContactStiffness", true, false );
  }

}


void FractunatorBase::Initialize( NodeManager& nodeManager,
                                  EdgeManagerT& edgeManager,
                                  FaceManagerT& faceManager,
                                  ElementManagerT& elementManager)
{
  std::map< std::string, lSet >::const_iterator setMap = nodeManager.m_Sets.find( m_separableNodeSet );
  array<integer>& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
  isSeparable = 1;


  // process nodes on the interior
  if( setMap!=nodeManager.m_Sets.end() )
  {
    isSeparable = 0;
    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
    {
      if( setMap->second.count(a)==1 )
        isSeparable[a] = 1;
    }
  }

  array<integer>& isFaceSeparable = faceManager.GetFieldData<int>("isSeparable");
  isFaceSeparable = 1;

  std::map< std::string, lSet >::const_iterator faceSetMap = faceManager.m_Sets.find( m_separableFaceSet );

  if( faceSetMap!=faceManager.m_Sets.end() )
  {
    isFaceSeparable = 0;
    for( localIndex a=0 ; a<faceManager.DataLengths() ; ++a )
    {
      if( faceSetMap->second.count(a)==1 )
        isFaceSeparable[a] = 1;
    }
  }



  array<real64>& separationCoeff = faceManager.GetFieldData<realT>("separationCoeff");
  separationCoeff = 0.0;



  if( faceManager.m_cohesiveZone )
    faceManager.m_cohesiveZone->resize( faceManager.DataLengths() );


  if ( !m_tounessSetByInitialCondition)
  {
    array<real64>& faceToughness = faceManager.GetFieldData<realT>("faceToughness");
    faceToughness = m_rockToughness;
  }

  array<real64>* faceFailStress =  faceManager.GetFieldDataPointer<realT>("faceFailStress");

  if ( !m_failStressSetByInitialCondition && faceFailStress != NULL)
  {
    (*faceFailStress) = m_failstress;
  }


  // Set an initial value for contact stress if necessary (this will get
  // over-written later)
  array<real64>* faceContactStiffness = faceManager.GetFieldDataPointer<realT>("faceContactStiffness");
  if (faceContactStiffness != NULL)
  {
    (*faceContactStiffness) = m_kJn;
  }

  array<integer> * const degreeFromCrack = nodeManager.GetFieldDataPointer<int>("degreeFromCrack");
  if( degreeFromCrack!=nullptr )
  {
    (*degreeFromCrack) = 1000;
  }



//  faceManager.m_cohesiveZones->Copy( *faceManager.m_cohesiveZone,
//                                     faceManager.m_cohesiveZones,
//                                     faceManager.m_numFaces );
}


void FractunatorBase::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  m_verbose = hdn.GetAttributeOrDefault<unsigned int>("verbose",0);
  m_failstress = hdn.GetAttributeOrDefault<realT>("failstress",0);
  m_separableNodeSet = hdn.GetAttributeStringOrDefault("separableNodeSet","");
  m_separableFaceSet = hdn.GetAttributeStringOrDefault("separableFaceSet","");
  m_failgap = hdn.GetAttributeOrDefault<realT>("failgap",0);
  if(isZero(m_failgap))
    throw GPException("The failgap parameter must be non-zero!");
  m_checkInterval = hdn.GetAttributeOrDefault<int>("fractureFlag",1);
  m_rockToughness = hdn.GetAttributeOrDefault<realT>("rockToughness",1e6);
  m_maxTurnAngle = hdn.GetAttributeOrDefault<realT>("maxTurnAngle",91.0);
  m_maxTurnAngle = m_maxTurnAngle / 180.0 * 3.1415926535;
  m_allowVacuumFrac = hdn.GetAttributeOrDefault<int>("allowVacuumFrac",0);
  m_thresholdForEvaluateFace = hdn.GetAttributeOrDefault<realT>("thresholdForEvaluateFace",0.5);
  m_displacementBasedSIF = hdn.GetAttributeOrDefault<realT>("displacementBasedSIF", false);
  m_kJn = hdn.GetAttributeOrDefault<realT>("normalJointStiffness", 1.0e10);

  //backward compatibility
  if (m_separableNodeSet.size() == 0)
    m_separableNodeSet = hdn.GetAttributeStringOrDefault("separableSet","");
  if (m_separableFaceSet.size() == 0)
    m_separableFaceSet = m_separableNodeSet;

  m_dfnPrefix = hdn.GetAttributeString("dfnPrefix");

}



// It is implemented in Fractunator3 and Fractunator2D for 3D and 2D
// respectively.
// The two implementations are quite different in logic and it's better keep
// them separate.
//
//int FractunatorBase::SeparationDriver( NodeManagerT& nodeManager,
//                                       EdgeManagerT& edgeManager,
//                                       FaceManagerT& faceManager,
//                                       ExternalFaceManagerT&
// externalFaceManager,
//                                       ElementManagerT& elementManager,
//                                       SpatialPartition& partition,
//                                       const bool prefrac,
//                                       const realT time)
//{
//
//
//  array<lSet> nodesToRupturedFaces;
//  array<lSet> edgesToRupturedFaces;
//
//  UpdateRuptureStates( nodeManager,
//                       edgeManager,
//                       faceManager,
//                       elementManager,
//                       nodesToRupturedFaces,
//                       edgesToRupturedFaces,
//                       prefrac );
//
//  int rval = 0;
//  int rank ;
//  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
////  array<MaterialBaseStateDataT*>&  temp =
// elementManager.m_ElementRegions["PM1"].m_materialStates;
//
//  const array<integer>& isNodeGhost =
// nodeManager.GetFieldData<FieldInfo::ghostRank>();
//  const array<integer>& isSeparable =
// nodeManager.GetFieldData<int>("isSeparable");
//  const array<integer>& layersFromDomainBoundary =
// nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
//
//  // process nodes on the interior
//  {
//    ModifiedObjectLists modifiedObjects;
//    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
//    {
//
////      const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);
//
//      if( layersFromDomainBoundary[a]>1 &&
//          isSeparable[a] &&
//          isNodeGhost[a]<0 &&
//          nodeManager.m_toElementsRelation[a].size()>1 &&
//          CheckSplitability( a, nodeManager, faceManager, edgeManager,
// prefrac) > 0 )
////        &&           nodesToRupturedFaces[parentNodeIndex].size()>0 )
//      {
//        rval += ProcessNode( a, nodeManager, edgeManager, faceManager,
// externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces,
// elementManager, modifiedObjects, prefrac ) ;
//      }
//    }
//  }
//
//
//  for( int color=0 ; color<partition.NumColor() ; ++color )
//  {
//    ModifiedObjectLists modifiedObjects;
//    if( partition.Color() == color )
//    {
//
//      // process "near-boundary" nodes
//      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
//      {
//
////        const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);
//
//        if( layersFromDomainBoundary[a]<=1 &&
//            isSeparable[a] &&
//            isNodeGhost[a]<0 &&
//            nodeManager.m_toElementsRelation[a].size()>1  &&
//            CheckSplitability( a, nodeManager, faceManager, edgeManager,
// prefrac) > 0 )
//            //&&  nodesToRupturedFaces[parentNodeIndex].size()>0 )
//        {
//          rval += ProcessNode( a, nodeManager, edgeManager, faceManager,
// externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces,
// elementManager, modifiedObjects, prefrac ) ;
//        }
//      }
//    }
//
//
//
//    // TODO need to add to rval as a result of this communication
//    partition.ModifyGhostsAndNeighborLists( modifiedObjects );
//
//    // If a face is split by a domains that does not own this face, the
// rupture state for the virtual face will not be communicated to the owner.
//    // The following is to fix this problem.
//    array<integer>* vfaceRuptureState =
// m_virtualFaces.GetFieldDataPointer<int>( "ruptureState" );
//    if (vfaceRuptureState != NULL)
//    {
//      const array<integer>& faceRuptureState = faceManager.GetFieldData<int>(
// "ruptureState" );
//      for( localIndex kf=0 ; kf<(*vfaceRuptureState).size() ; ++kf )
//      {
//        if( faceRuptureState[kf]==2 )
//        {
//          (*vfaceRuptureState)[kf]=2;
//        }
//      }
//
//    }
//
//
//  }
//
//
//
//  /*
//  for( std::map< std::string, ElementRegionT >::iterator
// i=elementManager.m_ElementRegions.begin() ;
//       i != elementManager.m_ElementRegions.end() ; ++i )
//  {
//    i->second.CalculateNodalMasses( nodeManager ) ;
//  }
//*/
//
//  return rval;
//}



bool FractunatorBase::ProcessNode( const localIndex nodeID,
                                   NodeManager& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager,
                                   ExternalFaceManagerT& externalFaceManager,
                                   array<lSet>& nodesToRupturedFaces,
                                   array<lSet>& edgesToRupturedFaces,
                                   ElementManagerT& elementManager,
                                   ModifiedObjectLists& modifiedObjects,
                                   const bool prefrac)
{
  bool didSplit = false;
  bool fracturePlaneFlag = true;

  while( fracturePlaneFlag )
  {
    lSet facialRupturePath;
    std::map<localIndex,int> edgeLocations;
    std::map<localIndex,int> faceLocations;
    std::map< std::pair< ElementRegionT*, localIndex >, int> elemLocations;


    fracturePlaneFlag = FindFracturePlanes(  nodeID,
                                             nodeManager,
                                             edgeManager,
                                             faceManager,
                                             nodesToRupturedFaces,
                                             edgesToRupturedFaces,
                                             facialRupturePath,
                                             edgeLocations,
                                             faceLocations,
                                             elemLocations );
    if( fracturePlaneFlag )
    {
      didSplit = true;
      PerformFracture( nodeID,
                       nodeManager,
                       edgeManager,
                       faceManager,
                       externalFaceManager,
                       elementManager,
                       modifiedObjects,
                       nodesToRupturedFaces,
                       edgesToRupturedFaces,
                       facialRupturePath,
                       edgeLocations,
                       faceLocations,
                       elemLocations );
    }
  }

  return didSplit;
}


//TODO: Need to clean this function.  Already moved its funcationality to
// another function in Fractunator3.
// But I am not sure about why we needed one of the loops.
// Will keep it here for now to be safe.
void FractunatorBase::UpdateRuptureStates( NodeManager& nodeManager,
                                           EdgeManagerT& edgeManager,
                                           FaceManagerT& faceManager,
                                           ElementManagerT& elementManager,
                                           array<lSet>& nodesToRupturedFaces,
                                           array<lSet>& edgesToRupturedFaces,
                                           const bool prefrac )
{

  nodesToRupturedFaces.resize(nodeManager.DataLengths());
  edgesToRupturedFaces.resize(edgeManager.DataLengths());


  if (!prefrac)
  {
    faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, this->m_failstress );
  }
  else
  {
    // During prefrac, we only need this to calculate the stress on faces.  We
    // don't mark rupture faces based on the stress state.
    faceManager.UpdateRuptureStates( elementManager, nodeManager, std::string(), 1.0e100 );
  }

  const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );

  array<integer>& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );

  // assign the values of the nodeToRupturedFaces and edgeToRupturedFaces
  // arrays.
  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    if( faceRuptureState[kf] == 1 && childFaceIndex[kf].size()==0 )
    {
      for( localIndex a=0 ; a<faceManager.m_toNodesRelation[kf].size() ; ++a )
      {
        const localIndex nodeIndex = faceManager.m_toNodesRelation[kf][a];
        nodesToRupturedFaces[nodeIndex].insert( kf );
      }

      for( localIndex a=0 ; a<faceManager.m_toEdgesRelation[kf].size() ; ++a )
      {
        const localIndex edgeIndex = faceManager.m_toEdgesRelation[kf][a];
        edgesToRupturedFaces[edgeIndex].insert( kf );
      }
    }
  }

  for( array<integer>::iterator i=edgeManager.m_isExternal.begin() ; i!=edgeManager.m_isExternal.end() ; ++i )
  {
    if( *i == -1 )
    {
      *i = 1;
      throw GPException("edgeManager.m_isExternal=-1. Call Pengcheng if you see this error");
      // I couldn't figure out why we need this loop. I am putting this
      // exception here so that when it is actually invoked, we will know why we
      // use it.
    }
  }


  array<lSet>::iterator i=nodesToRupturedFaces.begin();
  array<integer>::iterator j=nodeManager.GetFieldData<int>("numberOfRupturedFaces").begin();

  for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a, ++i, ++j )
  {
    *j = i->size();
  }



}


void FractunatorBase::WriteSilo( SiloFile& siloFile,
                                 const int cycleNum,
                                 const realT problemTime,
                                 const bool isRestart )
{

  std::string subDirectory = "Fractunator";
  std::string rootDirectory = "/Fractunator";
  siloFile.MakeSubDirectory( subDirectory, rootDirectory );
  DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());



  siloFile.DBWriteWrapper("m_verbose",m_verbose);
  siloFile.DBWriteWrapper("m_failstress",m_failstress);
  siloFile.DBWriteWrapper("m_failgap",m_failgap);
  siloFile.DBWriteWrapper("m_separableNodeSet",m_separableNodeSet);
  siloFile.DBWriteWrapper("m_separableFaceSet",m_separableFaceSet);


  WriteSiloDerived( siloFile, cycleNum, problemTime, isRestart );

  DBSetDir(siloFile.m_dbFilePtr, "..");

}

void FractunatorBase::ReadSilo( const SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart )
{

  if( DBSetDir(siloFile.m_dbFilePtr, "Fractunator" ) != -1 )
  {
    siloFile.DBReadWrapper("m_verbose",m_verbose);
    siloFile.DBReadWrapper("m_failstress",m_failstress);
    siloFile.DBReadWrapper("m_separableNodeSet",m_separableNodeSet);
    siloFile.DBReadWrapper("m_separableFaceSet",m_separableFaceSet);
    siloFile.DBReadWrapper("m_failgap",m_failgap);

    ReadSiloDerived( siloFile, cycleNum, problemTime, isRestart );

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }


}

int FractunatorBase::CheckOrphanElement (FaceManagerT& faceManager,
                                         localIndex iFace)
{
  array<integer>& ruptureState = faceManager.GetFieldData<int>("ruptureState");

  localIndex iEle;
  int flagOrphan = 0;
  for (array< std::pair< ElementRegionT*, localIndex > >::iterator iter = faceManager.m_toElementsRelation[iFace].begin() ;
       iter != faceManager.m_toElementsRelation[iFace].end() ; ++iter )
  {
    const ElementRegionT* elementRegion = iter->first;
    iEle= iter->second;
    unsigned int nRuptureFace = 0;
    localIndex* facelist = elementRegion->m_toFacesRelation[iEle];
    for (localIndex a=0 ; a<elementRegion->m_toFacesRelation.Dimension(1) ; ++a )
    {
      localIndex jFace = facelist[a];
      if ( (ruptureState[jFace] == 1 || faceManager.m_isExternal[jFace] == 1) && jFace != iFace)
      {
        nRuptureFace +=1;
      }
    }

    if (nRuptureFace == elementRegion->m_toFacesRelation.Dimension(1) - 1 )
      flagOrphan = 1;

  }
  return flagOrphan;
}

void FractunatorBase::MarkBirthTime( FaceManagerT& faceManager,
                                     ModifiedObjectLists& modifiedObjects,
                                     const realT time)
{
  array<real64>& birthTime = faceManager.GetFieldData<realT>("birthTime");

  for( lSet::const_iterator i=modifiedObjects.newFaces.begin() ; i!=modifiedObjects.newFaces.end() ; ++i )
  {
    birthTime[*i] = time;
    birthTime[faceManager.GetParentIndex(*i)] = time;
  }
}

//Fu: This is a rough correction, but better than no-correction.
void FractunatorBase::CorrectSplitNodalMass (NodeManager& nodeManager,
                                             localIndex node0,
                                             localIndex node1)
{
  array<realT>& mass = nodeManager.GetFieldData<FieldInfo::mass> ();

  realT totalMass = mass[node0] + mass[node1];
  mass[node0] = totalMass * nodeManager.m_toElementsRelation[node0].size() /
                (nodeManager.m_toElementsRelation[node0].size() + nodeManager.m_toElementsRelation[node1].size());
  mass[node1] = totalMass * nodeManager.m_toElementsRelation[node1].size() /
                (nodeManager.m_toElementsRelation[node0].size() + nodeManager.m_toElementsRelation[node1].size());

}

realT FractunatorBase::MinimumToughnessOnEdge( const localIndex edgeID,
                                               EdgeManagerT& edgeManager,
                                               FaceManagerT& faceManager)
{
  realT val = std::numeric_limits<realT>::max();
  array<real64>& faceToughness = faceManager.GetFieldData<realT>("faceToughness");

  for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[edgeID].begin() ;
       iface!=edgeManager.m_toFacesRelation[edgeID].end() ; ++iface )
  {
    val = std::min(val, faceToughness[*iface]);
  }

  return val;
}

realT FractunatorBase::MinimumToughnessOnNode( const localIndex nodeID,
                                               NodeManager& nodeManager,
                                               FaceManagerT& faceManager)
{
  realT val = std::numeric_limits<realT>::max();
  array<real64>& faceToughness = faceManager.GetFieldData<realT>("faceToughness");

  for (lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[nodeID].begin() ;
       iface!=nodeManager.m_nodeToFaceMap[nodeID].end() ; ++iface)
  {
    val = std::min(val, faceToughness[*iface]);
  }

  return val;
}

void FractunatorBase::MarkDiscreteFractureNetworkFaces(FaceManagerT& faceManager,
                                                       SpatialPartition& partition)
{
  if (partition.m_rank == 0)
    std::cout << "Marking DFN ID." << std::endl;
  array<integer>& dfnIndexMap = faceManager.GetFieldData<int>("DFN_Index");
  dfnIndexMap = -1;
  int dfnCount = 0;

  for( std::map<std::string,lSet>::const_iterator ii=faceManager.m_Sets.begin() ; ii!=faceManager.m_Sets.end() ; ++ii )
  {
    std::string setName = ii->first;

    if (setName.find(m_dfnPrefix) == 0)
    {
      lSet set = ii->second;
      dfnCount++;

      for( lSet::const_iterator jj=set.begin() ; jj!=set.end() ; ++jj)
      {
        if (dfnIndexMap[*jj] != -1)
          std::cout << "Warning: In rank " << partition.m_rank << " face " << *jj << " has been tagged with DFN ID " << dfnIndexMap[*jj] <<
            " but now we are retagging it with DFN ID " << dfnCount << std::endl;
        dfnIndexMap[*jj] = dfnCount;
      }
    }
  }
}
