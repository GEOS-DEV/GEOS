//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file EdgeManagerT.cpp
 * @author settgast1
 * @date Jun 22, 2011
 */

#include "EdgeManagerT.h"

#include "DataStructures/VectorFields/NodeManagerT.h"
#include "FaceManager.hpp"

EdgeManagerT::EdgeManagerT():
ObjectDataStructureBaseT(ObjectDataStructureBaseT::EdgeManager),
m_toNodesRelation(this->m_FixedOneToManyMaps["edgesToNodes"]),
m_toFacesRelation(this->m_UnorderedVariableOneToManyMaps["edgesToFaces"])
{
  m_toNodesRelation.resize2( 0, 2 );
  // TODO Auto-generated constructor stub

}

EdgeManagerT::~EdgeManagerT()
{
  // TODO Auto-generated destructor stub
}



void EdgeManagerT::BuildEdges( const FaceManagerT& faceManager, const NodeManager& nodeManager )
{
  if (faceManager.DataLengths() == 0 || nodeManager.DataLengths() == 0)
    return;

  localIndex numMultiNodeEdges = 0;
  Array1dT<lArray1d>& faceToEdgeMap = faceManager.m_toEdgesRelation;

  // this will be used to hold a list pairs that store the 2nd node in the edge, and the edge number.
  // they are keyed the edges lowest node.
  std::map< localIndex , Array1dT<std::pair<localIndex, localIndex> > > edgesByLowestNode;

  lSet singleNodeEdges;

  // For 2D problems, there is a one-to-one relationship between nodes and edges.
  if (faceManager.m_toNodesRelation[0].size() == 2)
  {
    m_toNodesRelation.resize2( 0, 1 );
    for( localIndex kn=0 ; kn<nodeManager.DataLengths() ; ++kn )
    {
      singleNodeEdges.insert(kn);
    }
  }


  // loop over all the faces
  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    const lArray1d::size_type numNodesInFace = faceManager.m_toNodesRelation[kf].size();

    if( numNodesInFace > 2 )
    {
      const lArray1d& nodeList = faceManager.m_toNodesRelation[kf];


      localIndex node0, node1, temp;

      // loop over all the nodes in the face. there will be an edge for each node.
      for( lArray1d::size_type a=0 ; a<numNodesInFace ; ++a )
      {
        // sort the nodes in order of index value
        node0 = nodeList[a];
        if( a == numNodesInFace-1 )
        {
          node1 = nodeList[0];
        }
        else
        {
          node1 = nodeList[a+1];
        }

        if( node0 > node1 )
        {
          temp = node0;
          node0 = node1;
          node1 = temp;
        }


        // check to see if the edge is in the edgesByLowestNode array (i.e. it is already registered)


        bool duplicate = false;
        Array1dT<std::pair<localIndex, localIndex> >& edgesWithSameFirstNode = edgesByLowestNode[node0];
        Array1dT<std::pair<localIndex, localIndex> >::iterator i=edgesWithSameFirstNode.begin();
        for(  ; i!=edgesWithSameFirstNode.end() ; ++i )
        {
          localIndex existingNode1 = i->first;
          localIndex existingEdgeIndex = i->second;
          if( existingNode1 == node1 ) // node1 is has already been entered...thus the edge already had been processed
          {
            duplicate = true;
            faceToEdgeMap(kf).push_back(existingEdgeIndex);

            break;
          }
        }

        // if the edge is not duplicate, then we will assign a new pair into the edgesByLowestNode array
        if( !duplicate )
        {
          edgesByLowestNode[node0].push_back( std::make_pair(node1,numMultiNodeEdges) );
          faceToEdgeMap(kf).push_back(numMultiNodeEdges);
          ++numMultiNodeEdges;
        }
      }
    }
    else if( numNodesInFace == 2 )
    {
      const lArray1d& nodeList = faceManager.m_toNodesRelation[kf];

      const localIndex node0 = nodeList[0];
      const localIndex node1 = nodeList[1];

      faceToEdgeMap(kf).push_back(node0);
      faceToEdgeMap(kf).push_back(node1);

    }
  }

  // all edge data is stored in the edgesByLowest node and faceManager::faceToEdgeMap arrays. So now we can
  // just extract the data into the edgeManager structures
  this->resize( numMultiNodeEdges + singleNodeEdges.size() );

  if (faceManager.m_toNodesRelation[0].size() > 2)
  {
    // fill edgesToNodes array by using data in the edgesByLowestNode array
    for( std::map< localIndex , Array1dT<std::pair<localIndex, localIndex> > >::iterator a=edgesByLowestNode.begin() ;
        a!=edgesByLowestNode.end() ; ++a )
    {
      Array1dT<std::pair<localIndex, localIndex> >& edgesWithSameFirstNode = a->second;
      const localIndex& node0 = a->first;
      for( Array1dT<std::pair<localIndex, localIndex> >::iterator i=edgesWithSameFirstNode.begin() ;
          i!=edgesWithSameFirstNode.end() ; ++i )
      {
        const localIndex& node1 = i->first;
        const localIndex& edgeIndex = i->second;

        m_toNodesRelation(edgeIndex,0) = node0;
        m_toNodesRelation(edgeIndex,1) = node1;

        nodeManager.m_nodeToEdgeMap(node0).insert(edgeIndex);
        nodeManager.m_nodeToEdgeMap(node1).insert(edgeIndex);

        //      std::cout<<"m_edgesToNodes("<<edgeIndex<<",*) = "<<node0<<' '<<node1<<std::endl;
      }
    }
  }
  else
  {
    for( localIndex edgeIndex=0 ; edgeIndex<singleNodeEdges.size() ; ++edgeIndex )
    {
      m_toNodesRelation(edgeIndex,0) = edgeIndex;
      nodeManager.m_nodeToEdgeMap(edgeIndex).insert(edgeIndex);
    }

  }


  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    const lArray1d::size_type numEdgesInFace = faceToEdgeMap(kf).size();

    // loop over all the nodes in the face. there will be an edge for each node.
    for( lArray1d::size_type a=0 ; a<numEdgesInFace ; ++a )
    {
      localIndex edgeIndex = faceToEdgeMap(kf)(a);
      m_toFacesRelation(edgeIndex).insert(kf);
    }
  }
  

  // fill edgesToNodes array by using data in the for single node edges
//  lSet::size_type edgeIndex = numMultiNodeEdges;
//  for( lSet::const_iterator i=singleNodeEdges.begin() ; i!=singleNodeEdges.end() ; ++i )
//  {
//    const localIndex& nodeIndex = *i;
//    m_toNodesRelation(edgeIndex,0) = nodeIndex;
//    nodeManager.m_nodeToEdgeMap(nodeIndex).insert(edgeIndex);
//
//    ++edgeIndex;
//  }
//

  // make sets from nodesets
  for( std::map< std::string, lSet >::const_iterator i = nodeManager.m_Sets.begin() ;
       i != nodeManager.m_Sets.end() ; ++i )
  {
    const std::string& setname = i->first;
    const lSet& set = i->second;
    this->ConstructSetFromSetAndMap( set, this->m_toNodesRelation, setname );
  }

}


/*

void EdgeManagerT::BuildEdges( const ElementManagerT& elementManager, const NodeManagerT& nodeManager )
{

  localIndex numEdges = 0;
  Array1dT<lArray1d>& faceToEdgeMap = faceManager.m_ToEdgesRelation;

  // this will be used to hold a list pairs that store the 2nd node in the edge, and the edge number.
  // they are keyed the edges lowest node.
  std::map< localIndex , Array1dT<std::pair<localIndex, localIndex> > > edgesByLowestNode;


  // loop over elements
  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator i=elementManager.m_ElementRegions.begin() ;
       i!=elementManager.m_ElementRegions.end() ; ++i )
  {
    ElementRegionT& elementRegion = i->second;

    FixedOneToManyRelation& elementToEdges = elementRegion.m_toEdgesMap;

    for( localIndex ke=0 ; ke<elementRegion.DataLengths() ; ++ke )
    {
      const lArray1d::size_type numNodesInElement = elementRegion.m_numNodesPerElem;
//      const localIndex* nodeList = elementRegion.m_.m_ToNodesRelation[kf];

    }
  }


  // loop over all the faces
  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    const lArray1d::size_type numNodesInFace = faceManager.m_ToNodesRelation[kf].size();
    const lArray1d& nodeList = faceManager.m_ToNodesRelation[kf];


    localIndex node0, node1, temp;

    // loop over all the nodes in the face. there will be an edge for each node.
    for( lArray1d::size_type a=0 ; a<numNodesInFace ; ++a )
    {
      // sort the nodes in order of index value
      node0 = nodeList[a];
      if( a == numNodesInFace-1 )
      {
        node1 = nodeList[0];
      }
      else
      {
        node1 = nodeList[a+1];
      }

      if( node0 > node1 )
      {
        temp = node0;
        node0 = node1;
        node1 = temp;
      }


      // check to see if the edge is in the edgesByLowestNode array (i.e. it is already registered)


      bool duplicate = false;
      Array1dT<std::pair<localIndex, localIndex> >& edgesWithSameFirstNode = edgesByLowestNode[node0];
      Array1dT<std::pair<localIndex, localIndex> >::iterator i=edgesWithSameFirstNode.begin();
      for(  ; i!=edgesWithSameFirstNode.end() ; ++i )
      {
        localIndex existingNode1 = i->first;
        localIndex existingEdgeIndex = i->second;
        if( existingNode1 == node1 ) // node1 is has already been entered...thus the edge already had been processed
        {
          duplicate = true;
          faceToEdgeMap(kf).push_back(existingEdgeIndex);

          break;
        }
      }

      // if the edge is not duplicate, then we will assign a new pair into the edgesByLowestNode array
      if( !duplicate )
      {
        edgesByLowestNode[node0].push_back( std::make_pair(node1,numEdges) );
        faceToEdgeMap(kf).push_back(numEdges);
        ++numEdges;
      }
    }
  }

  // all edge data is stored in the edgesByLowest node and faceManager::faceToEdgeMap arrays. So now we can
  // just extract the data into the edgeManager structures
  this->resize(numEdges);


  // fill edgesToNodes array by using data in the edgesByLowestNode array
  for( std::map< localIndex , Array1dT<std::pair<localIndex, localIndex> > >::iterator a=edgesByLowestNode.begin() ;
       a!=edgesByLowestNode.end() ; ++a )
  {
    Array1dT<std::pair<localIndex, localIndex> >& edgesWithSameFirstNode = a->second;
    const localIndex& node0 = a->first;
    for( Array1dT<std::pair<localIndex, localIndex> >::iterator i=edgesWithSameFirstNode.begin() ;
         i!=edgesWithSameFirstNode.end() ; ++i )
    {
      const localIndex& node1 = i->first;
      const localIndex& edgeIndex = i->second;

      m_toNodesRelation(edgeIndex,0) = node0;
      m_toNodesRelation(edgeIndex,1) = node1;

      nodeManager.m_nodeToEdgeMap(node0).insert(edgeIndex);
      nodeManager.m_nodeToEdgeMap(node1).insert(edgeIndex);

//      std::cout<<"m_edgesToNodes("<<edgeIndex<<",*) = "<<node0<<' '<<node1<<std::endl;
    }
  }

  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    const lArray1d::size_type numEdgesInFace = faceToEdgeMap(kf).size();

    // loop over all the nodes in the face. there will be an edge for each node.
    for( lArray1d::size_type a=0 ; a<numEdgesInFace ; ++a )
    {
      localIndex edgeIndex = faceToEdgeMap(kf)(a);
      m_toFacesRelation(edgeIndex).insert(kf);

    }
  }





  // make sets from nodesets
  for( std::map< std::string, lSet >::const_iterator i = nodeManager.m_Sets.begin() ;
       i != nodeManager.m_Sets.end() ; ++i )
  {
    const std::string& setname = i->first;
    const lSet& set = i->second;
    this->ConstructSetFromSetAndMap( set, this->m_toNodesRelation, setname );
  }

}
*/









/// Calculates the midpoint of the edge
void EdgeManagerT::EdgeCenter(const NodeManager& nodeManager, localIndex edge, R1Tensor& center)const{

  const Array1dT< R1Tensor >& refPosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT< R1Tensor >& displacement = nodeManager.GetFieldData<FieldInfo::displacement>();

  if (m_toNodesRelation.Dimension(1) >= 2)
  {
    const localIndex& node0 = m_toNodesRelation(edge,0);
    center =  refPosition[node0];
    center += displacement[node0];
    const localIndex& node1 = m_toNodesRelation(edge,1);
    center += refPosition[node1];
    center += displacement[node1];
    center *= 0.5;
  }
  else
  {
    const localIndex& node0 = m_toNodesRelation(edge,0);
    center =  refPosition[node0];
  }
}


/// Calculates the vector from node 0 to node 1
void EdgeManagerT::EdgeVector(const NodeManager& nodeManager, localIndex edge, R1Tensor& v) const{
  const Array1dT< R1Tensor >& refPosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT< R1Tensor >& displacement = nodeManager.GetFieldData<FieldInfo::displacement>();

  const localIndex& node1 = m_toNodesRelation(edge,1);
  v =  refPosition[node1];
  v += displacement[node1];
  const localIndex& node0 = m_toNodesRelation(edge,0);
  v -= refPosition[node0];
  v -= displacement[node0];
}

/// Returns the length of the edge
realT EdgeManagerT::EdgeLength(const NodeManager& nodeManager, localIndex edge) const{
  const Array1dT< R1Tensor >& refPosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT< R1Tensor >& displacement = nodeManager.GetFieldData<FieldInfo::displacement>();

  if (m_toNodesRelation.Dimension(1) >= 2)
  {
  const localIndex& node0 = m_toNodesRelation(edge,0);
  const localIndex& node1 = m_toNodesRelation(edge,1);
  R1Tensor v =  refPosition[node0];
  v += displacement[node0];
  v -= refPosition[node1];
  v -= displacement[node1];
  return v.L2_Norm();
  }
  else
  {
    return 1.0;
  }
}


void EdgeManagerT::SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject )
{
  // make sure that the reference object is a faceManger object
  referenceObject->CheckObjectType( ObjectDataStructureBaseT::FaceManager );

  // cast the referenceObject into a faceManager
  const FaceManagerT& faceManager = static_cast<const FaceManagerT&>(*referenceObject);

  // get the "isDomainBoundary" field from the faceManager...This should have been set already!
  const iArray1d& isFaceOnDomainBoundary = faceManager.GetFieldData<FieldInfo::isDomainBoundary>();

  // get the "isDomainBoudnary" field from for *this, and set it to zero
  iArray1d& isNodeOnDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
  isNodeOnDomainBoundary = 0;

  // loop through all faces
  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isFaceOnDomainBoundary[kf] == 1 )
    {
      lArray1d& faceToEdges = faceManager.m_toEdgesRelation(kf);

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( lArray1d::const_iterator a=faceToEdges.begin() ; a!=faceToEdges.end() ; ++a )
      {
        isNodeOnDomainBoundary(*a) = 1;
      }
    }
  }
}

bool EdgeManagerT::hasNode( const localIndex edgeID, const localIndex nodeID ) const
{
  if( m_toNodesRelation(edgeID,0) == nodeID || m_toNodesRelation(edgeID,1) == nodeID )
  {
    return true;
  }
  else
    return false;
}

localIndex EdgeManagerT::FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex nodeB, const NodeManager& nodeManager)
{
  localIndex val = std::numeric_limits<localIndex>::max();

  if (nodeA == nodeB) return (val);

  for( lSet::const_iterator iedge=nodeManager.m_nodeToEdgeMap[nodeA].begin() ;
      iedge!=nodeManager.m_nodeToEdgeMap[nodeA].end() ; ++iedge )
  {
    if (hasNode(*iedge, nodeB)) val = *iedge;
  }
  return(val);
}

void EdgeManagerT::SetIsExternal( const ObjectDataStructureBaseT* const referenceObject )
{
  // make sure that the reference object is a faceManger object
  referenceObject->CheckObjectType( ObjectDataStructureBaseT::FaceManager );

  // cast the referenceObject into a faceManager
  const FaceManagerT& faceManager = static_cast<const FaceManagerT&>(*referenceObject);

  // get the "isExternal" field from the faceManager...This should have been set already!
  const iArray1d& isExternalFace = faceManager.m_isExternal;

  // get the "isExternal" field from for *this, and set it to zero
  m_isExternal = 0;

  // loop through all faces
  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isExternalFace[kf] == 1 )
    {
      lArray1d& faceToEdges = faceManager.m_toEdgesRelation(kf);

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( lArray1d::const_iterator a=faceToEdges.begin() ; a!=faceToEdges.end() ; ++a )
      {
        m_isExternal(*a) = 1;
      }
    }
  }
}

#if 0
// It seems like this function is not used anywhere
void EdgeManagerT::SplitEdge( const localIndex indexToSplit,
                              const localIndex parentNodeIndex,
                              const localIndex childNodeIndex[2],
                              Array1dT<lSet>& nodesToEdges )
{

  localIndex newEdgeIndex[2] ;

  bool didSplit = SplitObject( indexToSplit, newEdgeIndex );

  // copy the parent edges edgeToNodes relation, replacing the parentNode with one of the new nodes.

  // loop over each new edge
  for( int ke=0 ; ke<2 ; ++ke )
  {


    // loop over each node on the edge
    for( int a=0 ; a<2 ; ++a )
    {
      // nodeIndex is the node on the parent edge
      const localIndex& nodeIndex = m_toNodesRelation(indexToSplit,a);

      // modify the edgesToNodes for new edges. They should point at a combination of one of the parents nodes,
      // one or two of the new nodes that were split.

      // if the nodeIndex==parentNodeIndex then the parent node is the target node
      if( nodeIndex == parentNodeIndex )
      {
        // adding the child node to the edgesToNodes map of the child edge
        m_toNodesRelation(newEdgeIndex[ke],a) = childNodeIndex[ke];
        std::cout<<"    m_edgesToNodes("<<newEdgeIndex[ke]<<","<<a<<") = "<<childNodeIndex[ke]<<std::endl;

        nodesToEdges(childNodeIndex[ke]).insert(newEdgeIndex[ke]);
//        std::cout<<"    nodesToEdges("<<childNodeIndex[ke]<<").insert("<<newEdgeIndex[ke]<<")"<<std::endl;
      }
      else
      {
        if( didSplit )
        {
          m_toNodesRelation(newEdgeIndex[ke],a) = nodeIndex;
          std::cout<<"    m_edgesToNodes("<<newEdgeIndex[ke]<<","<<a<<") = "<<nodeIndex<<std::endl;

          nodesToEdges(nodeIndex).insert(newEdgeIndex[ke]);
//          std::cout<<"    nodesToEdges("<<nodeIndex<<").insert("<<newEdgeIndex[ke]<<")"<<std::endl;
        }
      }
    }
  }



}

#endif

template< typename T_indices >
unsigned int EdgeManagerT::PackEdges( const T_indices& sendedges,
                                      const NodeManager& nodeManager,
                                      const FaceManagerT& faceManager,
                                      bufvector& buffer,
                                      const bool packConnectivityToGlobal,
                                      const bool packFields,
                                      const bool packMaps,
                                      const bool packSets  ) const
{

  unsigned int sizeOfPacked = 0;

  const std::string label = "EdgeData";
  sizeOfPacked += buffer.Pack(label);

  // pack data in the base
  sizeOfPacked += ObjectDataStructureBaseT::PackBaseObjectData( buffer, sendedges, packFields, packMaps, packSets, packConnectivityToGlobal );

  // pack the edge specific data
  for( typename T_indices::const_iterator edgeIndex=sendedges.begin() ; edgeIndex!=sendedges.end() ; ++edgeIndex )
  {
    const localIndex* const nodelist = m_toNodesRelation[*edgeIndex];

    int numnodes =m_toNodesRelation.Dimension(1);
    sizeOfPacked += buffer.Pack(numnodes);

    for( int a=0 ; a<numnodes ; ++a )
    {
      globalIndex gnode = GLOBALINDEX_MAX;
      if( packConnectivityToGlobal )
      {
        gnode = nodeManager.m_localToGlobalMap(nodelist[a]) ;
      }
      else
      {
        gnode = nodelist[a];
      }
      sizeOfPacked += buffer.Pack(gnode);
    }
  }


  const Array1dT<lSet>* const edgesToFlowFaces = GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
  if( edgesToFlowFaces != NULL )
  {
    for( typename T_indices::const_iterator edgeIndex=sendedges.begin() ; edgeIndex!=sendedges.end() ; ++edgeIndex )
    {
      const lSet& edgeToFlowFaces = (*edgesToFlowFaces)[*edgeIndex];
//      std::cout<<"pack edgeIndex, size = "<<this->m_localToGlobalMap(*edgeIndex)<<", "<<edgeToFlowFaces.size()<<" ";


      if( packConnectivityToGlobal )
      {
        sizeOfPacked += buffer.PackGlobal( edgeToFlowFaces, faceManager.m_localToGlobalMap );
//        std::cout<<"global"<<std::endl;
      }
      else
      {
        sizeOfPacked += buffer.Pack( edgeToFlowFaces );
//        std::cout<<"local"<<std::endl;
      }
    }
  }


  return sizeOfPacked;
}
template unsigned int EdgeManagerT::PackEdges( const lSet&, const NodeManager&, const FaceManagerT&, bufvector&, const bool, const bool, const bool, const bool ) const;
template unsigned int EdgeManagerT::PackEdges( const lArray1d&, const NodeManager&, const FaceManagerT&, bufvector&, const bool, const bool, const bool, const bool ) const;



unsigned int EdgeManagerT::UnpackEdges( const char*& buffer,
                                        const NodeManager& nodeManager,
                                        const FaceManagerT& faceManager,
                                        lArray1d& edgeReceiveLocalIndices,
                                        const bool unpackConnectivityToLocal,
                                        const bool unpackFields,
                                        const bool unpackMaps,
                                        const bool unpackSets  )
{
  unsigned int sizeOfUnpacked = 0;

  const std::string label = "EdgeData";
  std::string temp;

  sizeOfUnpacked += bufvector::Unpack( buffer, temp );
  if( label.compare(temp)!=0 )
  {
    throw GPException("EdgeManagerT::UnpackEdges: buffer location incorrect\n");
  }

  lArray1d newEdgeLocalIndices;
  // unpack data from base object
  sizeOfUnpacked += ObjectDataStructureBaseT::UnpackBaseObjectData( buffer, edgeReceiveLocalIndices, newEdgeLocalIndices, unpackFields, unpackMaps, unpackSets, unpackConnectivityToLocal );

  const lArray1d::size_type numReceivedEdges = edgeReceiveLocalIndices.size();
  // unpack face specific data
  for( lArray1d::size_type ke=0 ; ke<numReceivedEdges ; ++ke )
  {
    int numnodes;
    sizeOfUnpacked += bufvector::Unpack( buffer, numnodes );
    for( int a=0 ; a<numnodes ; ++a )
    {
      globalIndex gnode;
      sizeOfUnpacked += bufvector::Unpack( buffer, gnode );

      if( unpackConnectivityToLocal )
      {
        const localIndex lnode = stlMapLookup( nodeManager.m_globalToLocalMap, gnode );
        m_toNodesRelation(edgeReceiveLocalIndices(ke),a) = lnode;
      }
      else
      {
        m_toNodesRelation(edgeReceiveLocalIndices(ke),a) = gnode;
      }
    }

  }


  Array1dT<lSet>* const edgesToFlowFaces = GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
  if( edgesToFlowFaces != NULL )
  {
    for( lArray1d::iterator edgeIndex=edgeReceiveLocalIndices.begin() ; edgeIndex!=edgeReceiveLocalIndices.end() ; ++edgeIndex )
    {
      lSet newEdgeToFlowFaces;
      lSet& edgeToFlowFaces = (*edgesToFlowFaces)[*edgeIndex];


      if( unpackConnectivityToLocal )
      {
//        std::cout<<"global unpack edgeIndex, size = "<<this->m_localToGlobalMap(*edgeIndex)<<", ";
        sizeOfUnpacked += bufvector::UnpackGlobal( buffer, faceManager.m_globalToLocalMap, newEdgeToFlowFaces );
      }
      else
      {
//        std::cout<<"local unpack edgeIndex, size = "<<this->m_localToGlobalMap(*edgeIndex)<<", ";
        sizeOfUnpacked += bufvector::Unpack( buffer, newEdgeToFlowFaces );

        gSet globals;
        for( lSet::iterator i=edgeToFlowFaces.begin() ; i!=edgeToFlowFaces.end() ; ++i )
        {
          globalIndex gi = faceManager.m_localToGlobalMap[*i];
          globals.insert(gi);
        }
        edgeToFlowFaces.clear();
        edgeToFlowFaces.insert( globals.begin(), globals.end() );
      }
//      std::cout<<edgeToFlowFaces.size()<<std::endl;

      edgeToFlowFaces.insert( newEdgeToFlowFaces.begin(), newEdgeToFlowFaces.end() );
    }
  }


  return sizeOfUnpacked;

}



void EdgeManagerT::ConnectivityFromGlobalToLocal( const lSet& indices,
                                                  const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                                  const std::map<globalIndex,localIndex>& faceGlobalToLocal )
{
  for( lSet::const_iterator ke=indices.begin() ; ke!=indices.end() ; ++ke )
  {
    for( unsigned int a=0 ; a<m_toNodesRelation.Dimension(1) ; ++a )
    {
      const globalIndex gnode = m_toNodesRelation(*ke,a);
      const localIndex lnode = stlMapLookup( nodeGlobalToLocal, gnode );
      m_toNodesRelation(*ke,a) = lnode;
    }
  }

  Array1dT<lSet>* const edgesToFlowFaces = GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
  if( edgesToFlowFaces != NULL )
  {
    for( lSet::const_iterator ke=indices.begin() ; ke!=indices.end() ; ++ke )
    {
      lSet& edgeToFlowFaces = (*edgesToFlowFaces)[*ke];
      lSet newSet;
      for( lSet::iterator faceIndex=edgeToFlowFaces.begin() ; faceIndex!=edgeToFlowFaces.end() ; ++faceIndex )
      {

        std::map<globalIndex,localIndex>::const_iterator MapIter = faceGlobalToLocal.find( static_cast<globalIndex>(*faceIndex) );
        if( MapIter!=faceGlobalToLocal.end()  )
        {
//          const localIndex faceLocalIndex = stlMapLookup( faceGlobalToLocal, static_cast<globalIndex>(*faceIndex) );
          newSet.insert( MapIter->second );

        }
      }
      edgeToFlowFaces = newSet;
    }

  }

}

// Fu note on 20130416: Looks like this was temporary.  This function is now taken care of by element region. 
// We will keep it here for a while and delete it later.
//void EdgeManagerT::UpdateEdgeExternalityFromSplit( const FaceManagerT& faceManager,
//                                                 const lSet& newEdgeIndices,
//                                                 const lSet& modifiedEdgeIndices )
//{
//  lSet allEdges;
//  allEdges.insert( newEdgeIndices.begin(), newEdgeIndices.end() );
//  allEdges.insert( modifiedEdgeIndices.begin(), modifiedEdgeIndices.end() );
//
//
//  for( lSet::const_iterator edgeIndex=allEdges.begin() ; edgeIndex!=allEdges.end() ; ++edgeIndex )
//  {
//
//    for( lSet::const_iterator iface=m_toFacesRelation[*edgeIndex].begin() ;
//        iface!=m_toFacesRelation[*edgeIndex].end() ; ++iface )
//    {
//      if (faceManager.m_isExternal[*iface] == 1)
//      {
//        m_isExternal[*edgeIndex] =1;
//      }
//    }
//
//  }
//}


void EdgeManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                     Array1dT<gArray1d>& objectToCompositionObject )
{
  compositionObjectManager.CheckObjectType( ObjectDataStructureBaseT::NodeManager );

  iArray1d& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();

  for( localIndex kf=0 ; kf<DataLengths() ; ++kf )
  {
    if( isDomainBoundary(kf) != 0 )
    {
      gArray1d temp;

      for( lArray2d::size_type a=0 ; a<m_toNodesRelation.Dimension(1) ; ++a )
      {
        const localIndex lnode = m_toNodesRelation(kf,a);
        const globalIndex gnode = compositionObjectManager.m_localToGlobalMap(lnode);
        temp.push_back( gnode );
      }
      std::sort( temp.begin(), temp.end() );
      temp.insert( temp.begin(), this->m_localToGlobalMap[kf]);
      objectToCompositionObject.push_back(temp);
    }
  }
}


void EdgeManagerT::AddToEdgeToFaceMap( const FaceManagerT& faceManager,
                                       const lArray1d& newFaceIndices )
{
  // loop over all faces in list
  for( size_t kf=0 ; kf<newFaceIndices.size() ; ++kf )
  {
    localIndex lfi = newFaceIndices(kf);
    // now iterate over the faceToNodeMap (i.e. all nodes in the faceToNodeMap)
    for( lArray1d::const_iterator a=faceManager.m_toEdgesRelation[lfi].begin() ;
         a!=faceManager.m_toEdgesRelation[lfi].end() ; ++a )
    {
      // enter the value of the face index into the nodeToFace map
      m_toFacesRelation[*a].insert(lfi);
    }
  }
}

void EdgeManagerT::SetLayersFromDomainBoundary(const NodeManager& nodeManager)
{
  iArray1d& layersEdge = this->GetFieldData<int>("LayersFromDomainBoundary");
  const iArray1d& layersNode = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");

  for (localIndex ie = 0; ie!= this->DataLengths(); ++ie)
  {
    for( unsigned int a=0 ; a<m_toNodesRelation.Dimension(1) ; ++a )
    {
      if (a==0)
      {
        layersEdge[ie] = layersNode[m_toNodesRelation(ie,a)];
      }
      else
      {
        layersEdge[ie] = std::max(layersNode[m_toNodesRelation(ie,a)], layersEdge[ie]);
      }
    }
  }
}
