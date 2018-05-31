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
 * @file Fractunator.cpp
 * @author settgast1
 * @date Jul 14, 2011
 */

#include "Fractunator.h"
#include <limits.h>

Fractunator::Fractunator():
  m_verbose(0),
  m_failstress(0)
{
  // TODO Auto-generated constructor stub

}

Fractunator::~Fractunator()
{
  // TODO Auto-generated destructor stub
}


void Fractunator::RegisterFieldsAndMaps( NodeManager& nodeManager,
                                         EdgeManagerT& edgeManager,
                                         FaceManagerT& faceManager )
{



  nodeManager.AddMap<array<lArray1d> >("nodesToRupturedFaces");

  edgeManager.AddMap<array<lArray1d> >("edgesToRupturedFaces");

  // the faceManager's rutpureState will be used with the following definitions:
  //   ruptureState = 0 means not reached rupture criteria
  //                = 1 means has reached rupture criteria
  //                = 2 means that the face meet conditions to be part of a
  // rupture plane
  faceManager.AddKeylessDataField<int>( "ruptureState",true, true );

  nodeManager.AddKeylessDataField<int>("numberOfRupturedFaces",true,true);



  nodeManager.AddMap<array<lArray1d> >("childIndices");
  nodeManager.AddMap<lArray1d>("parentIndex");
  OneToOneRelation& parentIndexNodes = nodeManager.GetOneToOneMap( "parentIndex" );
  parentIndexNodes = LOCALINDEX_MAX;

  edgeManager.AddMap<array<lArray1d> >("childIndices");
  edgeManager.AddMap<lArray1d>("parentIndex");
  lArray1d& parentIndexEdge = edgeManager.GetOneToOneMap( "parentIndex" );
  parentIndexEdge = LOCALINDEX_MAX;

  faceManager.AddMap<array<lArray1d> >("childIndices");
  faceManager.AddMap<lArray1d>("parentIndex");
  lArray1d& parentIndexFace = faceManager.GetOneToOneMap( "parentIndex" );
  parentIndexFace = LOCALINDEX_MAX;


  faceManager.AddKeylessDataField<R1Tensor>( "maxTraction",true, true );

}


void Fractunator::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  m_verbose = hdn.GetAttributeOrDefault<unsigned int>("verbose",0);
  m_failstress = hdn.GetAttributeOrDefault<realT>("failstress",0);

}


void Fractunator::SeparationDriver( NodeManager& nodeManager,
                                    EdgeManagerT& edgeManager,
                                    FaceManagerT& faceManager,
                                    ExternalFaceManagerT& externalFaceManager,
                                    ElementManagerT& elementManager)
{
  UpdateRuptureStates( nodeManager,
                       edgeManager,
                       faceManager,
                       elementManager );


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  const OrderedVariableOneToManyRelation& nodeToRupturedFaces = nodeManager.GetVariableOneToManyMap( "nodesToRupturedFaces" );

  const array<lArray1d>& childNodeIndex = nodeManager.GetVariableOneToManyMap( "childIndices" );

  const array<integer>& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();
  for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
  {

    if( nodeToRupturedFaces[a].size() && childNodeIndex[a].size()==0 && isNodeGhost[a]<0 )
    {

      lSet facialRupturePath;
      std::map<localIndex,int> edgeLocations;
      std::map<localIndex,int> faceLocations;
      std::map< std::pair< ElementRegionT*, localIndex >, int> elemLocations;



      if( FindFracturePlanes(  a,
                               nodeManager,
                               edgeManager,
                               faceManager,
                               facialRupturePath,
                               edgeLocations,
                               faceLocations,
                               elemLocations ) )
      {
        PerformFracture( a,
                         nodeManager,
                         edgeManager,
                         faceManager,
                         elementManager,
                         facialRupturePath,
                         edgeLocations,
                         faceLocations,
                         elemLocations );
      }
    }
  }

  /*
     for( std::map< std::string, ElementRegionT >::iterator
        i=elementManager.m_ElementRegions.begin() ;
       i != elementManager.m_ElementRegions.end() ; ++i )
     {
     i->second.CalculateNodalMasses( nodeManager ) ;
     }
   */
}


void Fractunator::UpdateRuptureStates( NodeManager& nodeManager,
                                       EdgeManagerT& edgeManager,
                                       FaceManagerT& faceManager,
                                       ElementManagerT& elementManager )
{
  faceManager.UpdateRuptureStates( elementManager, nodeManager, std::string(), this->m_failstress );

  OrderedVariableOneToManyRelation& nodeToRupturedFaces = nodeManager.GetVariableOneToManyMap( "nodesToRupturedFaces" );
  OrderedVariableOneToManyRelation& edgesToRupturedFaces = edgeManager.GetVariableOneToManyMap( "edgesToRupturedFaces" );
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
        nodeToRupturedFaces[nodeIndex].push_back( kf );
      }

      for( localIndex a=0 ; a<faceManager.m_toEdgesRelation[kf].size() ; ++a )
      {
        const localIndex edgeIndex = faceManager.m_toEdgesRelation[kf][a];
        edgesToRupturedFaces[edgeIndex].push_back( kf );
      }
    }
  }

  for( array<integer>::iterator i=edgeManager.m_isExternal.begin() ; i!=edgeManager.m_isExternal.end() ; ++i )
  {
    if( *i == -1 )
      *i = 1;
  }


  array<lArray1d>::iterator i=nodeToRupturedFaces.begin();
  array<integer>::iterator j=nodeManager.GetFieldData<int>("numberOfRupturedFaces").begin();

  for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a, ++i, ++j )
  {
    *j = i->size();
  }


}


bool Fractunator::FindFracturePlanes( const localIndex nodeID,
                                      const NodeManager& nodeManager,
                                      const EdgeManagerT& edgeManager,
                                      const FaceManagerT& faceManager,
                                      lSet& separationPathFaces,
                                      std::map<localIndex,int>& edgeLocations,
                                      std::map<localIndex,int>& faceLocations,
                                      std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{

  const lArray1d& nodeToRupturedFaces = nodeManager.GetVariableOneToManyMap( "nodesToRupturedFaces" )[nodeID];
  const lSet& nodeToFaces = nodeManager.GetUnorderedVariableOneToManyMap( "nodeToFaceMap" )[nodeID];
  const lSet& nodeToEdges = nodeManager.GetUnorderedVariableOneToManyMap( "nodeToEdgeMap" )[nodeID];

  const set< std::pair<ElementRegionT*,localIndex> >& nodesToElements = nodeManager.m_toElementsRelation[nodeID] ;

  const array<lArray1d>& edgesToRupturedFaces = edgeManager.GetVariableOneToManyMap( "edgesToRupturedFaces" );
  const array<integer>& isEdgeExternal = edgeManager.m_isExternal;

  const array<lArray1d>& faceToEdges = faceManager.m_toEdgesRelation;


  // **** local working arrays ****

  // array to hold the faces ready for rupture
  set<localIndex> nodeToRuptureReadyFaces;
  for( lArray1d::const_iterator i=nodeToRupturedFaces.begin() ;
       i!=nodeToRupturedFaces.end() ; ++i )
  {
    nodeToRuptureReadyFaces.insert(*i);
  }



  if( nodeToRuptureReadyFaces.size() == 0 )
    return false;



  //
  std::map< localIndex, set<localIndex> > edgesToRuptureReadyFaces;
  for( lSet::const_iterator edgeIndex=nodeToEdges.begin() ; edgeIndex!=nodeToEdges.end() ; ++edgeIndex )
  {
    if( !(edgesToRupturedFaces[*edgeIndex].empty()) )
      edgesToRuptureReadyFaces[*edgeIndex].insert( edgesToRupturedFaces[*edgeIndex].begin(), edgesToRupturedFaces[*edgeIndex].end() );
  }

  // need a map from faces to edges attached to the node
  std::map< localIndex, std::pair<localIndex,localIndex> > localFacesToEdges;


  for( lSet::const_iterator kf=nodeToFaces.begin() ; kf!=nodeToFaces.end() ; ++kf )
  {
    localIndex edge[2] = { INT_MAX,INT_MAX };
    int count = 0;
    for( lArray1d::const_iterator ke=faceToEdges[*kf].begin() ; ke!=faceToEdges[*kf].end() ; ++ke )
    {
      if( edgeManager.hasNode( *ke, nodeID ) )
      {
        edge[count++] = *ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      throw GPException("Fractunator::FindFracturePlanes: invalid edge.");
    }


    localFacesToEdges[*kf] = std::make_pair<localIndex,localIndex>(edge[0],edge[1]);

    if( m_verbose ==2 )
      std::cout<<"localFacesToEdges["<<*kf<<"] = ( "<<localFacesToEdges[*kf].first<<", "<<localFacesToEdges[*kf].second<<" )"<<std::endl;
  }


  // ****
  // if the edge is not external, and the size of edgesToRupturedFaces is less
  // than 2, then the edge is a dead-end
  // as far as a rupture plane is concerned. The face associated with the edge
  // should be removed from the working
  // list of ruptured faces.

  // loop over all the edges
  for( lSet::const_iterator edgeIndex=nodeToEdges.begin() ; edgeIndex!=nodeToEdges.end() ; ++edgeIndex )
  {

    localIndex thisEdge = *edgeIndex;

//    std::cout<<thisEdge<<std::endl;
    // if the edge is internal and the edge is only attached to one ruptured
    // face
    while( isEdgeExternal[thisEdge]!=1 && edgesToRuptureReadyFaces[thisEdge].size()==1 )
    {
      // the index for the face that is a "dead end"
      localIndex deadEndFace = *(edgesToRuptureReadyFaces[thisEdge].begin());

      // get the edge on the other side of the face
      localIndex nextEdge = INT_MAX;
      if( localFacesToEdges[deadEndFace].first == thisEdge )
        nextEdge = localFacesToEdges[deadEndFace].second;
      else if( localFacesToEdges[deadEndFace].second == thisEdge )
        nextEdge = localFacesToEdges[deadEndFace].first;
      else
      {
//        return false;
        std::cout<<"nodeID, thisEdge = "<<nodeID<<", "<<thisEdge<<std::endl;
        std::cout<<"deadEndFace, localFacesToEdges[deadEndFace] = "<<", "<<deadEndFace<<", ( "<<localFacesToEdges[deadEndFace].first<<", "<<
          localFacesToEdges[deadEndFace].second<<" )"<<std::endl;

        std::cout<<"  nodeToEdges["<<nodeID<<"] = ( ";
        for( lSet::const_iterator i=nodeToEdges.begin() ; i!=nodeToEdges.end() ; ++i )
        {
          std::cout<<*i<<", ";
        }
        std::cout<<" )"<<std::endl;

        std::cout<<"  nodeToFaces["<<nodeID<<"] = ( ";
        for( lSet::const_iterator i=nodeToFaces.begin() ; i!=nodeToFaces.end() ; ++i )
        {
          std::cout<<*i<<", ";
        }
        std::cout<<" )"<<std::endl;

        std::cout<<"  faceToEdges["<<deadEndFace<<"] = ( ";
        for( lArray1d::const_iterator i=faceToEdges[deadEndFace].begin() ; i!=faceToEdges[deadEndFace].end() ; ++i )
        {
          std::cout<<*i<<", ";
        }
        std::cout<<" )"<<std::endl;

//        std::cout<<"  edgeToNodes[76922] = (
// "<<edgeManager.m_edgesToNodes(76922,0)<<",
// "<<edgeManager.m_edgesToNodes(76922,1)<<" )"<<std::endl;


//        throw GPException("Fractunator::FindFracturePlanes: Could not find the
// next edge when removing dead end faces.");
      }

      // delete the face from the working arrays
      edgesToRuptureReadyFaces[*edgeIndex].erase( deadEndFace );
      edgesToRuptureReadyFaces[nextEdge].erase(deadEndFace);
      nodeToRuptureReadyFaces.erase(deadEndFace);

      // if all the faces have been deleted, then go ahead and delete the top
      // level entry
      if( edgesToRuptureReadyFaces[*edgeIndex].empty() )
        edgesToRuptureReadyFaces.erase(*edgeIndex);
      if( edgesToRuptureReadyFaces[nextEdge].empty() )
        edgesToRuptureReadyFaces.erase(nextEdge);

      // now increment the "thisEdge" to point to the other edge on the face
      // that was just deleted
      thisEdge = nextEdge;
    }
  }

  if( nodeToRuptureReadyFaces.empty() )
  {
    return false;
  }



  // so now the working arrays have been purged of any faces that are on a
  // dead-end path. All remaining faces
  // are part of a separation plane...of course, there can be more than
  // one...which is bad. We will just take the first
  // path we find, and call this function again after the selected path is
  // processed. That will happen automatically
  // since the new nodes that are created will have higher node indices than the
  // current node, and will be checked for
  // separation prior to completion of the separation driver.



  // We now have to define the separation plane over which a node/face/edge will
  // be split, and all elements on one side
  // of the plane get one set of objects, and all elements on the other side get
  // the other set.


  // these are the edge and face where the separation surface starts.
  localIndex startingEdge = INT_MAX;
  localIndex startingFace = INT_MAX;


  // the startingEdge and startingFace needs to be set. It is best to start with
  // an external face if we have one.
  // loop over all edges that have an entry in edgesToRuptureReadyFaces (i.e.
  // all edges that are attached to a ruptured
  // node.

  for( std::map< localIndex, set<localIndex> >::const_iterator ke=edgesToRuptureReadyFaces.begin() ;
       ke!=edgesToRuptureReadyFaces.end() ; ++ke )
  {

    // make sure there is a face still attached to the edge, as it could have
    // been removed when we got rid of dead ends
    // ...actually, this shouldn't ever happen as we have already removed such
    // edges from the map.
    if( ke->second.size() > 0 )
    {

      startingEdge = ke->first;
      startingFace = *(ke->second.begin());
    }
    // if the size is 1, then the edge is only attached to one ruptured face,
    // which means that it is external. This is
    // the case that we want. The starting edge and face already were set above,
    // so just break at this point
    if( ke->second.size() == 1 && isEdgeExternal[ke->first]==1 )
    {
      break;
    }
  }

  // if the starting face was not set, then we don't have a rupture
  // surface....so just quit.
  if( startingFace==INT_MAX )
    return false;


  // now we start the process of setting the separation path. Begin by
  localIndex thisEdge = startingEdge;
  localIndex thisFace = startingFace;

  localIndex nextEdge;

//  separationPathFaces.insert( startingFace );

  // the seprationPath is used to hold combinations of edge and face
  std::map<localIndex,int> facesInPath;
  std::map<localIndex,int> edgesInPath;

  int numFacesInPath = 0;
  edgesInPath[thisEdge] = numFacesInPath;
  facesInPath[thisFace] = numFacesInPath++;

  // now walk from face->edge->face->edge etc. until we get to another external
  // edge, or back to the startingEdge.
  bool breakFlag = false;
  while ( !breakFlag )
  {
    // assign the other edge on the face as the next edge
    if( localFacesToEdges[thisFace].first == thisEdge )
    {
      nextEdge = localFacesToEdges[thisFace].second;
    }
    else if( localFacesToEdges[thisFace].second == thisEdge )
    {
      nextEdge = localFacesToEdges[thisFace].first;
    }
    else
    {
      throw GPException("Fractunator::FindFracturePlanes breakpoint 2");
    }

    // if we have reached an external face, or the edge we started with, then we
    // are done
    if( isEdgeExternal[nextEdge]==1 || edgesInPath.count(nextEdge)==1 )
    {
      const int startingIndex = edgesInPath[nextEdge];
      for( std::map<localIndex,int>::const_iterator kf=facesInPath.begin() ; kf!=facesInPath.end() ; ++kf )
      {
//        std::cout<<kf->first<<", "<<kf->second<<std::endl;
        if( kf->second >= startingIndex )
        {
          separationPathFaces.insert( kf->first );
        }
      }
      breakFlag = true;
    }
    else if( edgesToRuptureReadyFaces[nextEdge].size() )
    {
      // we need to pick another face attached to the "next edge"
      // increment the face and edge, and add to the separationPathFaces
      set<localIndex>::const_iterator iter_edgeToFace = edgesToRuptureReadyFaces[nextEdge].begin();
      if( *iter_edgeToFace == thisFace )
      {
        thisFace=*(++iter_edgeToFace);
      }
      else
      {
        bool pathFound = false;
        for( ; iter_edgeToFace!=edgesToRuptureReadyFaces[nextEdge].end() ; ++iter_edgeToFace )
        {
          if( *iter_edgeToFace == thisFace )
          {
            thisFace=*(--iter_edgeToFace);
            pathFound = true;
            break;
          }
        }
        if( pathFound == false )
          throw GPException("Fractunator::FindFracturePlanes breakpoint 3");
      }

      thisEdge = nextEdge;
//      separationPathFaces.insert( thisFace );
      edgesInPath[thisEdge] = numFacesInPath;
      facesInPath[thisFace] = numFacesInPath++;
    }
    else
    {
      std::cout<<"you are screwed"<<std::endl;
      throw GPException("Fractunator::FindFracturePlanes breakpoint 3b");
    }

  }



  // now we want to identify the objects on either side of the separation plane.
  // First we assign an array to indicate
  // whether a face/edge is on the fracture plane.
  for( lSet::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    edgeLocations[*ke] = INT_MIN;
  }

  for( lSet::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
  {
    faceLocations[*ke] = INT_MIN;
  }

  // now loop over the separation plane, and set the location arrays.
  for( lSet::const_iterator kf=separationPathFaces.begin() ; kf!=separationPathFaces.end() ; ++kf )
  {
    faceLocations[*kf] = -1;
    edgeLocations[localFacesToEdges[*kf].first] = -1;
    edgeLocations[localFacesToEdges[*kf].second] = -1;
  }


  for( set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodesToElements.begin() ; k!=nodesToElements.end() ; ++k )
  {
    elemLocations[*k] = INT_MIN;
  }

  SetLocations( 0, separationPathFaces, faceManager, nodesToElements, localFacesToEdges,
                edgeLocations, faceLocations, elemLocations );

  if( !(SetLocations( 1, separationPathFaces, faceManager, nodesToElements, localFacesToEdges,
                      edgeLocations, faceLocations, elemLocations )) )
  {
    return false;
  }



  bool fail = false;

  for( lSet::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    if( edgeLocations[*ke] == INT_MIN )
    {
      fail = true;
    }
  }
  for( lSet::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
  {
    if( faceLocations[*ke] == INT_MIN )
    {
      fail = true;
    }
  }

  if( fail )
  {
    return false;
    std::cout<<"Unset face or edge\n";
    std::cout<<"Separating node "<<nodeID<<" on separation plane consisting of faces: ";
    for( lSet::const_iterator i=separationPathFaces.begin() ; i!=separationPathFaces.end() ; ++i )
    {
      std::cout<<*i<<", ";
    }
    std::cout<<std::endl;


    for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator iter_elem=elemLocations.begin() ; iter_elem!=elemLocations.end() ;
         ++iter_elem )
    {
      const std::pair< ElementRegionT*, localIndex >& elem = iter_elem->first;

      ElementRegionT& elemRegion = *(elem.first);
      const localIndex elemIndex = elem.second;

      const int location = iter_elem->second;


      std::cout<<"Element "<<elemIndex<<" at location "<<location<<"\n";


      std::cout<<" faces->edges->nodes = ";
      for( int a=0 ; a<6 ; ++a )
      {
        localIndex faceIndex = elemRegion.m_toFacesRelation(elemIndex,a);

        if( a>0 )
          std::cout<<"                      = ";



        if( faceManager.m_parentIndex[faceIndex] != LOCALINDEX_MAX )
        {
          std::cout<<faceManager.m_parentIndex[faceIndex]<<"->";
        }
        std::cout<<faceIndex<<"[ ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
          std::cout<<edgeIndex<<", ";
        }
        std::cout<<" ] \n";

      }
      std::cout<<std::endl;

    }

    for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator i=elemLocations.begin() ; i!=elemLocations.end() ; ++i )
    {
      std::cout<<"( "<<i->first.second<<","<<i->second<<") , ";
    }
    std::cout<<std::endl;


    std::cout<<"faceLocations = ";
    for( std::map<localIndex,int>::const_iterator i=faceLocations.begin() ; i!=faceLocations.end() ; ++i )
    {
      std::cout<<"( "<<i->first<<","<<i->second<<") , ";
    }
    std::cout<<std::endl;

    std::cout<<"edgeLocations = ";
    for( std::map<localIndex,int>::const_iterator i=edgeLocations.begin() ; i!=edgeLocations.end() ; ++i )
    {
      std::cout<<"( "<<i->first<<","<<i->second<<") , ";
    }
    std::cout<<std::endl;
    throw GPException("Fractunator::FindFracturePlanes breakpoint 6...smart guy");
  }

  return true;
}



bool Fractunator::SetLocations( const int location,
                                const lSet& separationPathFaces,
                                const FaceManagerT& faceManager,
                                const set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                                std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                                std::map<localIndex,int>& edgeLocations,
                                std::map<localIndex,int>& faceLocations,
                                std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{
  bool rval = true;
  const localIndex separationFace = *(separationPathFaces.begin());
  set< std::pair<ElementRegionT*,localIndex> > elem0 ;
  set< std::pair<ElementRegionT*,localIndex> > processedElements ;
  const OneToOneRelation& parentFaceIndex = faceManager.GetOneToOneMap("parentIndex");

  // insert an element attached to the separation face
  elem0.insert( faceManager.m_toElementsRelation[separationFace][location] );
//  elemLocations[faceManager.m_FaceToElementMap[separationFace][location]] =
// location;

  if( m_verbose )
    std::cout<<"Setting Location "<<location<<std::endl;


  bool addedElem = false;
  do
  {
    addedElem = false;
    for( set< std::pair<ElementRegionT*,localIndex> >::iterator k=elem0.begin() ; k!=elem0.end() ; ++k )
    {
      // make sure that we have not already processed the element
      if( processedElements.count(*k)==0 )
      {
        processedElements.insert(*k);
        //      if( m_verbose ) std::cout<<"  processing Element
        // "<<k->second<<std::endl;
        for( localIndex kf=0 ; kf<k->first->m_toFacesRelation.Dimension(1) ; ++kf )
        {

          const localIndex faceIndex = k->first->m_toFacesRelation(k->second,kf);

          //        if( m_verbose ) std::cout<<"    processing Face
          // "<<faceIndex<<std::endl;

          std::map<localIndex,int>::iterator iterFace = faceLocations.find(faceIndex);

          if( iterFace != faceLocations.end() )
          {
            // make sure the face is not on the separation plane
            if( iterFace->second != -1 )
            {
              iterFace->second = location;

              if( edgeLocations[localFacesToEdges[faceIndex].first] == INT_MIN )
                edgeLocations[localFacesToEdges[faceIndex].first] =location;
              if( edgeLocations[localFacesToEdges[faceIndex].second] == INT_MIN )
                edgeLocations[localFacesToEdges[faceIndex].second] = location;

              // if the face has a parent face, then we need to set the location
              // for
              // the parent face as well
              localIndex linkingFaceIndex = faceIndex;

              if( parentFaceIndex[faceIndex] != LOCALINDEX_MAX )
              {
                // make sure the parent face is not on the sepration plane
                if( separationPathFaces.count(parentFaceIndex[faceIndex])==0 )
                {
                  linkingFaceIndex = parentFaceIndex[faceIndex];

                  // find the iterator to the parent face in the faceLocations
                  // map
                  std::map<localIndex,int>::iterator iterFaceP = faceLocations.find(linkingFaceIndex);
                  if( iterFaceP != faceLocations.end() )
                  {
                    iterFaceP->second = location;
                    if( edgeLocations[localFacesToEdges[linkingFaceIndex].first] == INT_MIN )
                      edgeLocations[localFacesToEdges[linkingFaceIndex].first] =location;
                    if( edgeLocations[localFacesToEdges[linkingFaceIndex].second] == INT_MIN )
                      edgeLocations[localFacesToEdges[linkingFaceIndex].second] = location;

                  }
                }
              }

              // now we add the element that is a neighbor to the "linkingFace"
              // of course, this only happens if there are more than one element
              // attached to the face.
              if( faceManager.m_toElementsRelation[linkingFaceIndex].size() > 1 )
              {
                const std::pair<ElementRegionT*,localIndex>& elemIndex0 = faceManager.m_toElementsRelation[linkingFaceIndex][0];
                const std::pair<ElementRegionT*,localIndex>& elemIndex1 = faceManager.m_toElementsRelation[linkingFaceIndex][1];

                // if the first element is the one we are on, and the element is
                // attached
                // to the splitting node, then add the second element to the
                // list.
                if( ( elemIndex0 == *k ) && ( nodesToElements.find(elemIndex1)!=nodesToElements.end() ) )
                {
                  elem0.insert(elemIndex1);
                  addedElem = true;
                }
                // if the second element is the one we are on, and the element
                // is attached
                // to the splitting node, then add the first element to the
                // list.
                else if( ( elemIndex1 == *k ) && ( nodesToElements.find(elemIndex0)!=nodesToElements.end() ) )
                {
                  elem0.insert(elemIndex0);
                  addedElem = true;
                }
              }
            }
          }


        }
        if( addedElem == true)
        {
          break;
        }
      }
    }
  }
  while(addedElem==true);


  for( set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=elem0.begin() ; k!=elem0.end() ; ++k  )
  {
    if( m_verbose )
      std::cout<<"  Setting Element "<<k->second<<" to location "<<location<<std::endl;

    if( elemLocations.find(*k) != elemLocations.end() )
    {
//      std::cout<<k->second<<", "<<elemLocations[*k]<<std::endl;
      if( elemLocations[*k]==INT_MIN )
      {
        elemLocations[*k] = location;
      }
      else
      {
        rval = false;
      }
    }

  }

  return rval;
}

void Fractunator::PerformFracture( const localIndex nodeID,
                                   NodeManager& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager,
                                   ElementManagerT& elementManager,
                                   const lSet& separationPathFaces,
                                   const std::map<localIndex,int>& edgeLocations,
                                   const std::map<localIndex,int>& faceLocations,
                                   const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{

  const OrderedVariableOneToManyRelation& childNodeIndex = nodeManager.GetVariableOneToManyMap( "childIndices" );
  const OrderedVariableOneToManyRelation& childEdgeIndex = edgeManager.GetVariableOneToManyMap( "childIndices" );
  const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );



  // generate a new pair of nodes, keeping the old one around for sentimental
  // purposes...or to be used in the
  // flow element.
  if( childNodeIndex[nodeID].size()==0 )
  {

    // ***** split all the objects first *****

    // split node
    localIndex newNodeIndex[2];
    if( m_verbose )
    {
      std::cout<<"\nSplitting node "<<nodeID<<" along separation plane faces ";
      for( lSet::const_iterator i=separationPathFaces.begin() ; i!=separationPathFaces.end() ; ++i )
      {
        std::cout<<*i<<", ";
      }
      std::cout<<std::endl;
    }


    // if the node doesn't have a parent, then it is an original node, and two
    // should be created from it. If it has a parent, then only one additional
    // node should be created, with the original faces and edges getting
    // reattached to the parent node.
    if( nodeManager.m_parentIndex(nodeID) == LOCALINDEX_MAX )
    {
      nodeManager.SplitObject(nodeID, newNodeIndex);
    }
    else
    {
      newNodeIndex[0] = nodeID;
      nodeManager.SplitObject(nodeID,newNodeIndex[1]);
    }

    if( m_verbose )
      std::cout<<"\nDone splitting node "<<nodeID<<" into nodes "<<newNodeIndex[0]<<" and "<<newNodeIndex[1]<<std::endl;
    nodeManager.m_toElementsRelation.resize( nodeManager.m_numNodes );


    // split edges
    array<integer>& flowEdgeType = edgeManager.GetFieldData<int>("flowEdgeType");
    lSet splitEdges;
    // loop over all edges connected to the node
    for( std::map<localIndex,int>::const_iterator iter_edge=edgeLocations.begin() ; iter_edge!=edgeLocations.end() ; ++iter_edge )
    {
      const localIndex& edgeIndex = iter_edge->first;
      const int& location = iter_edge->second;

      // if the edge is on the separation plane, then split it
      if( location == -1  )
      {
        localIndex newEdgeIndex[2];

        if( edgeManager.SplitObject( edgeIndex, newEdgeIndex ) )
        {
          if( m_verbose )
            std::cout<<"  Split edge "<<edgeIndex<<" into edges "<<childEdgeIndex[edgeIndex][0]<<" and "<<childEdgeIndex[edgeIndex][1]<<std::endl;

          splitEdges.insert( edgeIndex );

          for( int a=0 ; a<2 ; ++a )
          {
            edgeManager.m_toNodesRelation(childEdgeIndex[edgeIndex][0],a) = edgeManager.m_toNodesRelation(edgeIndex,a);
            edgeManager.m_toNodesRelation(childEdgeIndex[edgeIndex][1],a) = edgeManager.m_toNodesRelation(edgeIndex,a);
          }

          flowEdgeType[childEdgeIndex[edgeIndex][0]] = 0;
          flowEdgeType[childEdgeIndex[edgeIndex][1]] = 0;
          flowEdgeType[edgeIndex] = 1;
        }
      }
    }


    // split the faces
    array<integer>& ruptureState = faceManager.GetFieldData<int>("ruptureState");
    array<integer>& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
    lSet splitFaces;

    // loop over all faces attached to the nodeID
    for( std::map<localIndex,int>::const_iterator iter_face=faceLocations.begin() ; iter_face!=faceLocations.end() ; ++iter_face )
    {
      const localIndex& faceIndex = iter_face->first;
      const int& location = iter_face->second;
      // if the face is on the separation plane, then split it
      if( location == -1 )
      {
        localIndex newFaceIndex[2];

        if( faceManager.SplitObject( faceIndex, newFaceIndex ) )
        {
          if( m_verbose )
            std::cout<<"  Split face "<<faceIndex<<" into faces "<<childFaceIndex[faceIndex][0]<<" and "<<childFaceIndex[faceIndex][1]<<std::endl;

          splitFaces.insert( faceIndex );

          ruptureState[childFaceIndex[faceIndex][0]] = 0;
          ruptureState[childFaceIndex[faceIndex][1]] = 0;

          flowFaceType[childFaceIndex[faceIndex][0]] = 0;
          flowFaceType[childFaceIndex[faceIndex][1]] = 0;
          flowFaceType[faceIndex] = 1;

          faceManager.m_toEdgesRelation[newFaceIndex[0]] = faceManager.m_toEdgesRelation[faceIndex];
          faceManager.m_toEdgesRelation[newFaceIndex[1]] = faceManager.m_toEdgesRelation[faceIndex];

          faceManager.m_toNodesRelation[newFaceIndex[0]] = faceManager.m_toNodesRelation[faceIndex];
          faceManager.m_toNodesRelation[newFaceIndex[1]] = faceManager.m_toNodesRelation[faceIndex];
        }
      }
    }
    faceManager.m_toElementsRelation.resize( faceManager.m_numFaces );



    // ***** now correct all the relations between the objects *****

    /* To accomplish this annoying yet exceedingly important task, we will take
       a "top down"
     * approach. Note that this is a two way correction, i.e. if we are
     * correcting
     * elementToNodes, we also correct nodesToElements. This is summarized as:
     * 1) Loop over elements attached to the split node.
     *     2a) correct all relations between the single  element and the nodes.
     *     2b) Loop over all faces on the element
     *         3a) For each face, correct the face relations with the element
     *         3b) For each face, correct the face relations with the nodes
     *         3c) Loop over all edges on the face
     *             4a) For each edge, correct the face relations
     *             4b) for each edge, correct the node relations
     *
     *  The element location will define which side of the rupture everything
     *  is on.
     */

    // 1) loop over all elements attached to the nodeID
    for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator iter_elem=elemLocations.begin() ;
         iter_elem!=elemLocations.end() ; ++iter_elem )
    {
      const std::pair< ElementRegionT*, localIndex >& elem = iter_elem->first;

      ElementRegionT& elemRegion = *(elem.first);
      const localIndex elemIndex = elem.second;
      const int& location = iter_elem->second;

      if( m_verbose )
        std::cout<<"Element "<<elemIndex<<std::endl;

      // 2a) correct elementToNode and nodeToElement
      if( m_verbose )
        std::cout<<"  Looping over all nodes on element, and correcting node<->element maps:"<<std::endl;
      {
        // loop over all nodes on element
        if( m_verbose )
          std::cout<<"    m_ElementToNodeMap = ( ";
        localIndex* const nodelist = elemRegion.m_toNodesRelation.FirstIndexPointer( elemIndex );
        for( localIndex a=0 ; a<elemRegion.m_toNodesRelation.Dimension(1) ; ++a )
        {
          // if the node was just split
          if( nodelist[a] == nodeID )
          {

            if( m_verbose )
              std::cout<<nodelist[a]<<"->"<<newNodeIndex[ location ]<<", ";

            nodelist[a] = newNodeIndex[ location ];

            nodeManager.m_toElementsRelation[newNodeIndex[location]].insert(elem);
            nodeManager.m_toElementsRelation[nodeID].erase(elem);

          }
          else if( m_verbose )
            std::cout<<nodelist[a]<<", ";
        }
        if( m_verbose )
          std::cout<<")"<<std::endl;

        if( m_verbose )
        {
          for( localIndex a=0 ; a<elemRegion.m_toNodesRelation.Dimension(1) ; ++a )
          {
            if( m_verbose ) std::cout<<"    nodeToElemMaps["<<nodelist[a]<<"] = ( ";
            for( set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodelist[a]].begin() ;
                 k!=nodeManager.m_toElementsRelation[nodelist[a]].end() ; ++k )
            {
              std::cout<<k->second<<", ";
            }
            std::cout<<" )"<<std::endl;

          }
        }
      }



      // 2b) loop over all faces on element...actually it is all faces + any
      // parent faces that used to be attached to the element
      if( m_verbose )
      {
        std::cout<<"  Looping over all faces on element (parent and child):"<<std::endl;
      }


      localIndex* const elemToFaces = elemRegion.m_toFacesRelation.FirstIndexPointer(elemIndex);
      // we need to build a list of faces that is elemToFaces FOLLOWED by any
      // parent face of those indicated in elemToFaces


      lArray1d facelist;
      // first fill the facelist with the faces in elemToFaces
      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {
        facelist.push_back(elemToFaces[kf]);
      }
      // now fill the facelist with the parents of elemToFaces
      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {
        const localIndex parentFaceIndex = faceManager.m_parentIndex[elemToFaces[kf]];
        if( parentFaceIndex!=LOCALINDEX_MAX )
        {
          const std::map<localIndex,int>::const_iterator iterFace = faceLocations.find(parentFaceIndex);
          if( iterFace != faceLocations.end() )
          {
            if( iterFace->second == location )
              facelist.push_back(parentFaceIndex);
          }
        }
      }

      // Now we do a loop over the facelist and process all the faces
      for( localIndex kf=0 ; kf<facelist.size() ; ++kf )
//      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {
        // set both faceID and newFaceID to the parent face.
//        localIndex faceID = elemToFaces[kf];
        localIndex faceID = facelist[kf];
        localIndex newFaceID = faceID;

        // 3a) check to see if the face was split. If so, then we will need
        // to alter the face relation with the elements in both directions.
        if( splitFaces.count( faceID ) > 0 && kf<static_cast<unsigned int>(elemRegion.m_numFacesPerElement) )
        {

          // replace the parent face with the child face in elementToFace. Now
          // faceID is the parent face, and newFaceID is the child face.
          elemToFaces[kf] = childFaceIndex[ faceID ][location];
          newFaceID = elemToFaces[kf];

          // add the element to the child faceToElem
          faceManager.m_toElementsRelation[newFaceID].push_back( elem );

          // we are keeping the parent faceToElement map intact. There is no
          // harm
          // and we can know what elements a face originally came from
          //faceManager.m_FaceToElementMap[faceID].clear();

        }



        // 3b) correct faceToNodes and nodeToFaces
        array<lArray1d>& nodeToRupturedFaces = nodeManager.GetVariableOneToManyMap( "nodesToRupturedFaces" );
        array<lArray1d>& edgesToRupturedFaces = edgeManager.GetVariableOneToManyMap( "edgesToRupturedFaces" );
        array<integer>& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );


        if( m_verbose )
        {
          const localIndex parentFace = faceManager.m_parentIndex[newFaceID];
          if( parentFace!=LOCALINDEX_MAX )
          {
            std::cout<<"    m_FaceToNodeMap["<<parentFace<<"->"<<newFaceID<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToNodeMap["<<newFaceID<<"] = ( ";
          }
        }

        // loop over all nodes on the face.
        for( lArray1d::iterator nodeIndex=faceManager.m_toNodesRelation[newFaceID].begin() ; nodeIndex!=faceManager.m_toNodesRelation[newFaceID].end() ;
             ++nodeIndex )
        {
          if( m_verbose )
            std::cout<<*nodeIndex;
          // if the facenode is the one that is being split
          if( *nodeIndex == nodeID )
          {
            *nodeIndex = newNodeIndex[location];
            // if it is not a new face
            if( faceID == newFaceID )
            {
              nodeManager.m_nodeToFaceMap[nodeID].erase(faceID);

              nodeManager.m_nodeToFaceMap[*nodeIndex].insert(faceID);

              if( faceRuptureState[faceID ] )
                nodeToRupturedFaces[*nodeIndex].push_back(faceID);
            }
            else
            {
              nodeManager.m_nodeToFaceMap[*nodeIndex].insert(newFaceID);

              if( faceRuptureState[faceID ] )
                nodeToRupturedFaces[*nodeIndex].push_back(newFaceID);
            }

            if( m_verbose )
              std::cout<<"->"<<*nodeIndex<<", ";

          }
          else
          {
            nodeManager.m_nodeToFaceMap[*nodeIndex].insert(newFaceID);

            if( faceRuptureState[newFaceID ] )
              nodeToRupturedFaces[*nodeIndex].push_back(newFaceID);

            if( m_verbose )
              std::cout<<", ";
          }
        }
        if( m_verbose )
          std::cout<<")"<<std::endl;



        // faceToEdges
        if( m_verbose )
        {
          const localIndex parentFace = faceManager.m_parentIndex[newFaceID];
          if( parentFace!=LOCALINDEX_MAX )
          {
            std::cout<<"    m_FaceToEdgeMap["<<parentFace<<"->"<<newFaceID<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToEdgeMap["<<newFaceID<<"] = ( ";
          }
        }
        // loop over all edges on face
        for( lArray1d::iterator edgeIndex=faceManager.m_toEdgesRelation[newFaceID].begin() ; edgeIndex!=faceManager.m_toEdgesRelation[newFaceID].end() ;
             ++edgeIndex )
        {
          if( splitEdges.count( *edgeIndex ) > 0 )
          {
            if( faceID == newFaceID )
              edgeManager.m_toFacesRelation[*edgeIndex].erase(faceID);
            *edgeIndex = childEdgeIndex[*edgeIndex][location];
          }
          edgeManager.m_toFacesRelation[*edgeIndex].insert(newFaceID);
          if( m_verbose )
            std::cout<<*edgeIndex;

          if( faceRuptureState[newFaceID ] )
            edgesToRupturedFaces[*edgeIndex].push_back(newFaceID);


          //edgesToNodes
          if( m_verbose )
          {
            std::cout<<"(";
          }

          {
            localIndex* const nodeIndex = edgeManager.m_toNodesRelation.FirstIndexPointer(*edgeIndex);
            for( unsigned int a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a )
            {
              if( nodeIndex[a] == nodeID )
              {

                if( m_verbose )
                  std::cout<<nodeIndex[a];

                nodeIndex[a] = newNodeIndex[location];
                nodeManager.m_nodeToEdgeMap[nodeID].erase(*edgeIndex);

                if( m_verbose )
                  std::cout<<"->"<<nodeIndex[a]<<", ";

              }
              else if( m_verbose )
                std::cout<<nodeIndex[a]<<", ";

              nodeManager.m_nodeToEdgeMap[nodeIndex[a]].insert(*edgeIndex);
            }
            if( m_verbose )
              std::cout<<")";
          }
          if( m_verbose )
            std::cout<<", ";
        }
        if( m_verbose )
          std::cout<<")"<<std::endl;
      }
    }



    //**************************************************************************
    // THIS IS ALL JUST CONSISTENCY CHECKING
    //**************************************************************************


    if( m_verbose == 1 )
    {
      std::cout<<"CONSISTENCY CHECKING OF THE MAPS"<<std::endl;

      for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator iter_elem=elemLocations.begin() ; iter_elem!=elemLocations.end() ;
           ++iter_elem )
      {
        const std::pair< ElementRegionT*, localIndex >& elem = iter_elem->first;

        ElementRegionT& elemRegion = *(elem.first);
        const localIndex elemIndex = elem.second;


        lSet elemNodes;


        std::cout<<"Element "<<elemIndex<<"\n";
        std::cout<<" elementToNodes = ";
        for( int a=0 ; a<8 ; ++a )
        {
          elemNodes.insert(elemRegion.m_toNodesRelation(elemIndex,a));
          std::cout<<elemRegion.m_toNodesRelation(elemIndex,a)<<", ";
        }
        std::cout<<std::endl;



        std::cout<<" elementToFaces->edges->nodes = ";

        localIndex* const elemToFaces = elemRegion.m_toFacesRelation.FirstIndexPointer(elemIndex);

        lArray1d facelist;
        // first fill the facelist with the faces in elemToFaces
        for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
        {
          facelist.push_back(elemToFaces[kf]);
        }
        for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
        {
          const localIndex parentFaceIndex = faceManager.m_parentIndex[elemToFaces[kf]];
          if( parentFaceIndex!=LOCALINDEX_MAX )
          {
            facelist.push_back(parentFaceIndex);
          }
        }

        // Now we do a loop over the facelist and process all the faces
        for( int kf=0 ; kf<int(facelist.size()) ; ++kf )
        {
          lSet faceNodes;

          localIndex faceIndex;
          if( kf< elemRegion.m_numFacesPerElement )
            faceIndex = elemRegion.m_toFacesRelation(elemIndex,kf);
          else
            faceIndex = facelist[kf];

          if( kf>0 )
            std::cout<<"                              = ";

          if( kf>=elemRegion.m_numFacesPerElement )
            std::cout<<"P-";

          std::cout<<faceIndex<<"( ";
          for( int b=0 ; b<4 ; ++b )
          {
            localIndex faceNodeID = faceManager.m_toNodesRelation[faceIndex][b];
            faceNodes.insert(faceNodeID);
            if( elemNodes.count(faceNodeID) == 0 && kf<elemRegion.m_numFacesPerElement )
              std::cout<<"*";
            std::cout<<faceNodeID<<",";
          }
          std::cout<<" )      ";



          std::cout<<faceIndex<<"[ ";
          for( int b=0 ; b<4 ; ++b )
          {
            localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
            std::cout<<edgeIndex<<"( ";
            for( int c=0 ; c<2 ; ++c )
            {
              localIndex edgeNodeID = edgeManager.m_toNodesRelation(edgeIndex,c);
              if( elemNodes.count(edgeNodeID) == 0  && kf<elemRegion.m_numFacesPerElement )
                std::cout<<"*";
              if( faceNodes.count(edgeNodeID) == 0 )
                std::cout<<"#";
              std::cout<<edgeNodeID<<",";
            }
            std::cout<<" ), ";
          }
          std::cout<<" ] \n";

        }
        std::cout<<std::endl;

      }

    }

    if( m_verbose == 2 )
    {
      std::cout<<" elementToFaces->edges->nodes = ";
//      for( int a=0; a<6 ; ++a )
      {
        lSet faceNodes;
        localIndex faceIndex = 73847;

//        if( a>0 )
//          std::cout<<"                              = ";

        std::cout<<faceIndex<<"( ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex faceNodeID = faceManager.m_toNodesRelation[faceIndex][b];
          faceNodes.insert(faceNodeID);
//          if( elemNodes.count(nodeID) == 0 )
//            std::cout<<"*";
          std::cout<<faceNodeID<<",";
        }
        std::cout<<" )      ";



        std::cout<<faceIndex<<"[ ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
          std::cout<<edgeIndex<<"( ";
          for( int c=0 ; c<2 ; ++c )
          {
            localIndex edgeNodeID = edgeManager.m_toNodesRelation(edgeIndex,c);
//            if( elemNodes.count(nodeID) == 0 )
//              std::cout<<"*";
            if( faceNodes.count(edgeNodeID) == 0 )
              std::cout<<"#";
            std::cout<<edgeNodeID<<",";
          }
          std::cout<<" ), ";
        }
        std::cout<<" ] \n";

      }
      std::cout<<std::endl;

    }


    if( m_verbose == 2 )
    {
      // nodeToEdge
      array<lSet> tempNodesToEdges( nodeManager.m_numNodes );

      for( localIndex ke=0 ; ke<edgeManager.DataLengths() ; ++ke )
      {
        for( localIndex b= 0 ; b<edgeManager.m_toNodesRelation.Dimension(1) ; ++b )
        {
          localIndex nodeIndex = edgeManager.m_toNodesRelation(ke,b);
          tempNodesToEdges[nodeIndex].insert(ke);
        }
      }
      std::cout<<"Check NodeToEdge "<<std::endl;
      for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
      {
        std::cout<<"m_nodesToEdges["<<a<<"] = ( ";
        for( lSet::const_iterator iedge=nodeManager.m_nodeToEdgeMap[a].begin() ;
             iedge!=nodeManager.m_nodeToEdgeMap[a].end() ; ++iedge )
        {
          if( tempNodesToEdges[a].count(*iedge) == 0 )
            std::cout<<"*";

          std::cout<<*iedge<<", ";
        }
        std::cout<<")    (";

        for( lSet::const_iterator iedge=tempNodesToEdges[a].begin() ;
             iedge!=tempNodesToEdges[a].end() ; ++iedge )
        {
          if( nodeManager.m_nodeToEdgeMap[a].count(*iedge) == 0 )
            std::cout<<"*";

          std::cout<<*iedge<<", ";
        }
        std::cout<<")"<<std::endl;
      }


    }

    if( m_verbose == 2 )
    {
      // nodeToFace
      array<lSet> tempNodesToFaces( nodeManager.m_numNodes );
      for( localIndex kf=0 ; kf<faceManager.m_numFaces ; ++kf )
      {
        for( lArray1d::const_iterator b=faceManager.m_toNodesRelation[kf].begin() ;
             b!=faceManager.m_toNodesRelation[kf].end() ; ++b )
        {
          tempNodesToFaces[*b].insert(kf);
        }
      }
      std::cout<<"Check NodeToFace "<<std::endl;
      for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
      {
        std::cout<<"m_nodeToFaceMap["<<a<<"] = ( ";
        for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[a].begin() ;
             iface!=nodeManager.m_nodeToFaceMap[a].end() ; ++iface )
        {
          if( tempNodesToFaces[a].count(*iface) == 0 )
            std::cout<<"*";

          std::cout<<*iface<<", ";
        }
        std::cout<<")    (";

        for( lSet::const_iterator iface=tempNodesToFaces[a].begin() ;
             iface!=tempNodesToFaces[a].end() ; ++iface )
        {
          if( nodeManager.m_nodeToFaceMap[a].count(*iface) == 0 )
            std::cout<<"*";

          std::cout<<*iface<<", ";
        }
        std::cout<<")"<<std::endl;
      }

    }



    if( m_verbose == 2 )
    {

      // nodeToElement
      array<set<std::pair< ElementRegionT*, localIndex > > > tempNodesToElems( nodeManager.m_numNodes );
      for( std::map< std::string, ElementRegionT >::iterator ielem=elementManager.m_ElementRegions.begin() ;
           ielem!=elementManager.m_ElementRegions.end() ; ++ielem )
      {
        ElementRegionT& elemRegion = ielem->second;
        for( localIndex k=0 ; k<elemRegion.m_numElems ; ++k )
        {
          std::pair< ElementRegionT*, localIndex > elem = std::make_pair<ElementRegionT*, localIndex>(&elemRegion,k);

          for( localIndex a=0 ; a<elemRegion.m_toNodesRelation.Dimension(1) ; ++a )
          {
            tempNodesToElems[elemRegion.m_toNodesRelation(k,a)].insert(elem);
          }
        }
      }
      std::cout<<"Check NodeToElem "<<std::endl;
      for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
      {
        std::cout<<"m_NodeToElementMap["<<a<<"] = ( ";
        for( set<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=nodeManager.m_toElementsRelation[a].begin() ;
             ielem!=nodeManager.m_toElementsRelation[a].end() ; ++ielem )
        {
          if( tempNodesToElems[a].count(*ielem) == 0 )
            std::cout<<"*";

          std::cout<<ielem->second<<", ";
        }
        std::cout<<")    (";

        for( set<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=tempNodesToElems[a].begin() ;
             ielem!=tempNodesToElems[a].end() ; ++ielem )
        {
          if( nodeManager.m_toElementsRelation[a].count(*ielem) == 0 )
            std::cout<<"*";
          std::cout<<ielem->second<<", ";
        }
        std::cout<<")"<<std::endl;
      }


      // edgeToFace
      array<lSet> tempEdgeToFaces( edgeManager.DataLengths() );
      for( localIndex kf=0 ; kf<faceManager.m_numFaces ; ++kf )
      {
        for( lArray1d::const_iterator b=faceManager.m_toEdgesRelation[kf].begin() ;
             b!=faceManager.m_toEdgesRelation[kf].end() ; ++b )
        {
          tempEdgeToFaces[*b].insert(kf);
        }
      }
      std::cout<<"Check EdgeToFace "<<std::endl;
      for( localIndex ke=0 ; ke<edgeManager.DataLengths() ; ++ke )
      {
        std::cout<<"m_edgesToFaces["<<ke<<"] = ( ";
        for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[ke].begin() ;
             iface!=edgeManager.m_toFacesRelation[ke].end() ; ++iface )
        {
          if( tempEdgeToFaces[ke].count(*iface) == 0 )
            std::cout<<"*";

          std::cout<<*iface<<", ";
        }
        std::cout<<")    (";

        for( lSet::const_iterator iface=tempEdgeToFaces[ke].begin() ;
             iface!=tempEdgeToFaces[ke].end() ; ++iface )
        {
          if( edgeManager.m_toFacesRelation[ke].count(*iface) == 0 )
            std::cout<<"*";

          std::cout<<*iface<<", ";
        }
        std::cout<<")"<<std::endl;
      }

      // faceToElement
      OneToOneRelation& parentFaceIndex = faceManager.GetOneToOneMap("parentIndex");
      array<set<std::pair< ElementRegionT*, localIndex > > > tempFacesToElems( faceManager.m_numFaces );
      for( std::map< std::string, ElementRegionT >::iterator ielem=elementManager.m_ElementRegions.begin() ;
           ielem!=elementManager.m_ElementRegions.end() ; ++ielem )
      {
        ElementRegionT& elemRegion = ielem->second;
        for( localIndex k=0 ; k<elemRegion.m_numElems ; ++k )
        {
          std::pair< ElementRegionT*, localIndex > elem = std::make_pair<ElementRegionT*, localIndex>(&elemRegion,k);

          for( localIndex a=0 ; a<elemRegion.m_toFacesRelation.Dimension(1) ; ++a )
          {
            const localIndex faceID = elemRegion.m_toFacesRelation(k,a);
            tempFacesToElems[faceID].insert(elem);

            if( parentFaceIndex[faceID] != LOCALINDEX_MAX )
            {
              tempFacesToElems[parentFaceIndex[faceID]].insert(elem);
            }
          }
        }
      }
      std::cout<<"Check FacesToElem "<<std::endl;
      for( localIndex a=0 ; a<faceManager.m_numFaces ; ++a )
      {
        std::cout<<"m_FaceToElementMap["<<a<<"] = ( ";

        for( array<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=faceManager.m_toElementsRelation[a].begin() ;
             ielem!=faceManager.m_toElementsRelation[a].end() ; ++ielem )
        {
          if( tempFacesToElems[a].count(*ielem) == 0 )
            std::cout<<"*";

          std::cout<<ielem->second<<", ";
        }
        std::cout<<")    (";

        for( set<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=tempFacesToElems[a].begin() ;
             ielem!=tempFacesToElems[a].end() ; ++ielem )
        {

          if( faceManager.m_toElementsRelation[a].size() == 2 )
          {
            if( (faceManager.m_toElementsRelation[a][0] != *ielem) && (faceManager.m_toElementsRelation[a][1] != *ielem) )
              std::cout<<"*";
          }
          else if ( faceManager.m_toElementsRelation[a].size() )
          {
            if( (faceManager.m_toElementsRelation[a][0] != *ielem)  )
              std::cout<<"*";
          }
          else
          {
            std::cout<<"****";
          }


          std::cout<<ielem->second<<", ";
        }
        std::cout<<")"<<std::endl;
      }
    }

  }

}
