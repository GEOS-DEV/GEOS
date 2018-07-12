/*
 * SurfaceGenerator.cpp
 *
 *  Created on: Jul 3, 2018
 *      Author: settgast
 */

#include "SurfaceGenerator.hpp"

namespace geosx
{

SurfaceGenerator::SurfaceGenerator( const std::string& name,
                                    ManagedGroup * const parent ):
  SolverBase(name,parent)
{
  // TODO Auto-generated constructor stub

}

SurfaceGenerator::~SurfaceGenerator()
{
  // TODO Auto-generated destructor stub
}



void SurfaceGenerator::TimeStep( real64 const & time_n,
                                 real64 const & dt,
                                 const int cycleNumber,
                                 ManagedGroup * domain )
{

}



int SurfaceGenerator::SeparationDriver( NodeManager & nodeManager,
                                        EdgeManager & edgeManager,
                                        FaceManager & faceManager,
                                        ExternalFaceManager & externalFaceManager,
                                        ElementRegionManager & elementManager,
                                        SpatialPartition& partition,
                                        bool const prefrac,
                                        realT const time )
{
  array<lSet> nodesToRupturedFaces;
  array<lSet> edgesToRupturedFaces;

  //  for( map< ElementRegionManager::RegKeyType, ElementRegion >::iterator elementRegionIter = elementManager.m_ElementRegions.begin() ;
  //      elementRegionIter != elementManager.m_ElementRegions.end() ;
  //      ++elementRegionIter )
  //  {
  //    ElementRegion& elementRegion = elementRegionIter->second;
  //
  //    elementRegion.CalculateVelocityGradients(nodeManager);
  //  }
  //  integer_array& ruptureState = faceManager.GetFieldData<int>("ruptureState");

  if (!prefrac)
  {

    if (m_failCriterion >0 )  // Stress intensity factor based criterion and mixed criterion.
    {
      if (m_failCriterion == 1)
      {
        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, std::numeric_limits<realT>::max()); //this->m_failstress );
      }
      else
      {
        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, m_failstress);
      }

      real64_array& SIFonFace = faceManager.GetFieldData<realT>("SIFonFace");
      SIFonFace = std::numeric_limits<double>::min();

      if (m_failCriterion == 1)
      {
        integer_array& ruptureState = faceManager.GetFieldData<int>("ruptureState");
        for (localIndex a=0; a<faceManager.DataLengths(); ++a)
        {
          if (faceManager.m_toElementsRelation[a].size() == 2) ruptureState[a]=0;
        }
      }


      IdentifyRupturedFaces( nodeManager,
                             edgeManager,
                             faceManager,
                             elementManager,
                             partition,
                             prefrac );

    }

    else
    {
      UpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces,
                           prefrac );

    }
  }
  else  // In the prefrac call, we need this to get the stressNOnFace, which will be used in the initialization of contacts for preexisting fractures.
  {
    for (localIndex kf = 0; kf < faceManager.DataLengths(); ++kf) faceManager.CalculateStressOnFace(elementManager, nodeManager, kf);
  }


  if (prefrac)
  {
    ModifiedObjectLists modifiedObjects;
    CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, prefrac);
  }

  // We do this here to get the nodesToRupturedFaces etc.
  // The fail stress check inside has been disabled
  PostUpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces);

  int rval = 0;
  int rank ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //  array<MaterialBaseStateDataT*>&  temp = elementManager.m_ElementRegions["PM1"].m_materialStates;

  const integer_array& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();
  const integer_array& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
  const integer_array& layersFromDomainBoundary = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
  integer_array& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  const integer_array& isFaceGhost = faceManager.GetFieldData<FieldInfo::ghostRank>();



  // process nodes on the interior
  {
    ModifiedObjectLists modifiedObjects;
    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
    {
      //      const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);
      if( layersFromDomainBoundary[a]>1 &&
          (isSeparable[a] || prefrac)&&
          isNodeGhost[a]<0 &&
          nodeManager.m_toElementsRelation[a].size()>1 &&
          CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 ) //&&
        //          nodesToRupturedFaces[a].size()>0 )
      {
        rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
      }
    }
    if (m_failCriterion == 1)
    {
      const localIndex_array& primaryCandidateFace = faceManager.GetFieldData<localIndex>("primaryCandidateFace");


      for (localIndex a=0; a<faceManager.DataLengths(); ++a)
      {
        if (isFaceGhost[a]<0 && ruptureState[a] == -1 && ruptureState[primaryCandidateFace[a]] !=2 )
        {
          ruptureState[a] = 1;
        }
      }

      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
      {

        //      const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);

        if( layersFromDomainBoundary[a]>1 &&
            (isSeparable[a] || prefrac) &&
            isNodeGhost[a]<0 &&
            nodeManager.m_toElementsRelation[a].size()>1 &&
            CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 ) //&&
          //            nodesToRupturedFaces[a].size()>0 )
        {
          rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
        }
      }

    }

    CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);

    MarkBirthTime(faceManager, modifiedObjects, time);
  }


  for( int color=0 ; color<partition.NumColor() ; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( partition.Color() == color )
    {

      // process "near-boundary" nodes
      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
      {

        //        const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);

        if( layersFromDomainBoundary[a]<=1 &&
            (isSeparable[a] || prefrac) &&
            isNodeGhost[a]<0 &&
            nodeManager.m_toElementsRelation[a].size()>1 &&
            CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 )// &&
          //           nodesToRupturedFaces[a].size()>0 )
        {
          rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
        }
      }

      CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);

      MarkBirthTime(faceManager, modifiedObjects, time);

    }




    // TODO need to add to rval as a result of this communication
    partition.ModifyGhostsAndNeighborLists( modifiedObjects );

    // If a face is split by a domains that does not own this face, the rupture state for the virtual face will not be communicated to the owner.
    // The following is to fix this problem.
    integer_array* vfaceRuptureState = m_virtualFaces.GetFieldDataPointer<int>( "ruptureState" );
    if (vfaceRuptureState != NULL)
    {
      const integer_array& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );
      for( localIndex kf=0 ; kf<(*vfaceRuptureState).size() ; ++kf )
      {
        if( faceRuptureState[kf]==2 )
        {
          (*vfaceRuptureState)[kf]=2;
        }
      }

    }


  }

  if (m_failCriterion == 1)
  {
    const localIndex_array& primaryCandidateFace = faceManager.GetFieldData<localIndex>("primaryCandidateFace");

    {
      ModifiedObjectLists modifiedObjects;

      // Turn on rupture state for secondary fracture faces
      for (localIndex a=0; a<faceManager.DataLengths(); ++a)
      {
        if (isFaceGhost[a] < 0 && ruptureState[a] == -1 && ruptureState[primaryCandidateFace[a]] !=2 )
        {
          ruptureState[a] = 1;
          modifiedObjects.modifiedFaces.insert(a);
        }
      }
      partition.ModifyGhostsAndNeighborLists( modifiedObjects );
    }

    for( int color=0 ; color<partition.NumColor() ; ++color )
    {
      ModifiedObjectLists modifiedObjects;
      if( partition.Color() == color )
      {

        // process "near-boundary" nodes
        for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
        {

          if( layersFromDomainBoundary[a]<=1 &&
              (isSeparable[a] || prefrac) &&
              isNodeGhost[a]<0 &&
              nodeManager.m_toElementsRelation[a].size()>1 &&
              CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 ) //&&
            //              nodesToRupturedFaces[a].size()>0 )
          {
            rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
          }
        }
        CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);

        MarkBirthTime(faceManager, modifiedObjects, time);
      }



      // TODO need to add to rval as a result of this communication
      partition.ModifyGhostsAndNeighborLists( modifiedObjects );

      // If a face is split by a domains that does not own this face, the rupture state for the virtual face will not be communicated to the owner.
      // The following is to fix this problem.
      integer_array* vfaceRuptureState = m_virtualFaces.GetFieldDataPointer<int>( "ruptureState" );
      if (vfaceRuptureState != NULL)
      {
        const integer_array& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );
        for( localIndex kf=0 ; kf<(*vfaceRuptureState).size() ; ++kf )
        {
          if( faceRuptureState[kf]==2 )
          {
            (*vfaceRuptureState)[kf]=2;
          }
        }

      }
    }
  }


  if (prefrac && !m_dfnPrefix.empty())
  {
    MarkDiscreteFractureNetworkFaces(faceManager, partition);
  }


  /*
  for( map< std::string, ElementRegion >::iterator i=elementManager.m_ElementRegions.begin() ;
       i != elementManager.m_ElementRegions.end() ; ++i )
  {
    i->second.CalculateNodalMasses( nodeManager ) ;
  }
   */

  return rval;
}


REGISTER_CATALOG_ENTRY( SolverBase,
                        SurfaceGenerator,
                        std::string const &, ManagedGroup * const )

} /* namespace geosx */

