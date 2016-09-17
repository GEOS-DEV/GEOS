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
/*
 * NodeManagerT.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "NodeManager.hpp"

#include "IO/BinStream.h"
#include "ObjectManagers/FaceManagerT.h"
#include "ObjectManagers/EdgeManagerT.h"
#include "ObjectManagers/ElementManagerT.h"
#include "Utilities/Utilities.h"
#include <fstream>
#include "ElementRegionT.hpp"


// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager():
ObjectDataStructureBaseT( ObjectDataStructureBaseT::NodeManager ),
m_toElementsRelation(),
m_nodeToFaceMap(m_UnorderedVariableOneToManyMaps["nodeToFaceMap"]),
m_nodeToEdgeMap(m_UnorderedVariableOneToManyMaps["nodeToEdgeMap"])
//m_toCrackSurfacesRelation(m_UnorderedVariableOneToManyMaps["nodeToCrackSurfaceMap"])
{
  this->AddKeylessDataField<int>("LayersFromDomainBoundary", true, true );

  this->AddKeylessDataField<int>("isSeparable", true, true );
  iArray1d& isSeparable = this->GetFieldData<int>("isSeparable");
  isSeparable = 1;


  Initialize();
}



// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
/*
NodeManagerT::NodeManagerT( const NodeManagerT& init ):
ObjectDataStructureBaseT(init),
DataLengths()(this->m_DataLengths),
m_refposition(NULL),
m_displacement(NULL),
m_incrementalDisplacement(NULL),
m_velocity(NULL),
m_acceleration(NULL),
m_force(NULL),
m_mass(NULL),
m_toElementsRelation(init.m_toElementsRelation),
m_nodeToFaceMap(m_UnorderedVariableOneToManyMaps["nodeToFaceMap"]),
m_nodeToEdgeMap(m_UnorderedVariableOneToManyMaps["nodeToEdgeMap"])
{}
*/

// *********************************************************************************************************************
/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::~NodeManager()
{
}


void NodeManager::Initialize()
{
  this->AddKeyedDataField<FieldInfo::referencePosition>();

  this->AddKeyedDataField<FieldInfo::displacement>();

  this->AddKeyedDataField<FieldInfo::incrementalDisplacement>();

  this->AddKeyedDataField<FieldInfo::velocity>();

  this->AddKeyedDataField<FieldInfo::acceleration>();

  this->AddKeyedDataField<FieldInfo::force>();

  this->AddKeyedDataField<FieldInfo::mass>();
}

/**
 * @author R.R. Settgast
 * @param referenceObject optional object that will be used by the function to determine boundary status.
 *
 */
void NodeManager::SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject )
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
      lArray1d& faceToNodes = faceManager.m_toNodesRelation(kf);

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( lArray1d::const_iterator a=faceToNodes.begin() ; a!=faceToNodes.end() ; ++a )
      {
        isNodeOnDomainBoundary(*a) = 1;
      }
    }
  }
}

void NodeManager::SetLayersFromDomainBoundary( const int layer )
{
  const iArray1d& isNodeOnDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
  iArray1d& layersFromDomainBoundary = this->GetFieldData<int>("LayersFromDomainBoundary");


  if( layer == 0 )
  {
    layersFromDomainBoundary = INT_MAX;
    for( localIndex a=0 ; a<this->DataLengths() ; ++a )
    {
      if( isNodeOnDomainBoundary[a] == 1 )
      {
        layersFromDomainBoundary[a] = 0;
      }
    }
  }
  else
  {
    SetLayersFromDomainBoundary( layer-1 );

    nodeToElemType elementsOnLayer;

    for( localIndex a=0 ; a<this->DataLengths() ; ++a )
    {
      if( layersFromDomainBoundary[a] == (layer-1) )
      {

        elementsOnLayer.insert( m_toElementsRelation[a].begin(), m_toElementsRelation[a].end() );
      }
    }

    for( nodeToElemType::const_iterator iter_elem=elementsOnLayer.begin() ; iter_elem!=elementsOnLayer.end() ; ++iter_elem )
    {
      const ElementRegionT& elemRegion = *(iter_elem->first);
      const localIndex elemIndex = iter_elem->second;
      const localIndex numNodes = elemRegion.m_toNodesRelation.Dimension(1);
      const localIndex* const nodelist = elemRegion.m_toNodesRelation[elemIndex];

      for( localIndex a=0 ; a<numNodes ; ++a )
      {
        if( layersFromDomainBoundary[nodelist[a]] > (layer-1) )
        {
          layersFromDomainBoundary[nodelist[a]] = layer;
        }

      }

    }
  }


}




void NodeManager::SetIsExternal( const ObjectDataStructureBaseT* const referenceObject )
{
  // make sure that the reference object is a faceManger object
  referenceObject->CheckObjectType( ObjectDataStructureBaseT::FaceManager );

  // cast the referenceObject into a faceManager
  const FaceManagerT& faceManager = static_cast<const FaceManagerT&>(*referenceObject);

  // get the "isExternal" field from the faceManager...This should have been set already!
  const iArray1d& isExternalFace = faceManager.m_isExternal;

  // loop through all faces
  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isExternalFace(kf) == 1 )
    {
      lArray1d& faceToNodes = faceManager.m_toNodesRelation(kf);

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( lArray1d::const_iterator a=faceToNodes.begin() ; a!=faceToNodes.end() ; ++a )
      {
        m_isExternal(*a) = 1;
      }
    }
  }
}


// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @param destination local node number of destination node
 * @param source local node number of source node
 */
void NodeManager::CopyNode( const int destination, const int source )
{
  Array1dT< R1Tensor >& refPosition = GetFieldData<FieldInfo::referencePosition>();
  Array1dT< R1Tensor >& displacement = GetFieldData<FieldInfo::displacement>();
  Array1dT< R1Tensor >& incDisplacement = GetFieldData<FieldInfo::incrementalDisplacement>();
  Array1dT< R1Tensor >& velocity = GetFieldData<FieldInfo::velocity>();
  Array1dT< R1Tensor >& acceleration = GetFieldData<FieldInfo::acceleration>();
  Array1dT< R1Tensor >& force = GetFieldData<FieldInfo::force>();
  rArray1d& mass = GetFieldData<FieldInfo::mass>();

  refPosition[destination] = refPosition[source];
  displacement[destination] = displacement[source];
  incDisplacement[destination] = incDisplacement[source];
  velocity[destination] = velocity[source];
  acceleration[destination] = acceleration[source];
  force[destination] = force[source];
  mass[destination] = mass[source];

}

/**
 * @author R. Settgast
 * @param elementManager the element manager
 * This function constructs the nodeToElementMap by using the data in each element regions elementToNode map to
 * invert the relation.
 */
void NodeManager::ConstructNodeToElementMap( const ElementManagerT& elementManager )
{
  // because the nodeToElementMap is an odd creature, it is not managed by ObjectDataStructureBaseT...so we must
  // resize.
  m_toElementsRelation.resize(DataLengths());

  // iterate over all element regions
  for( std::map<ElementManagerT::RegKeyType, ElementRegionT>::const_iterator ielemRegion=elementManager.m_ElementRegions.begin() ;
       ielemRegion!=elementManager.m_ElementRegions.end() ; ++ielemRegion )
  {
    // the element region is the mapped value of the iterator
    const ElementRegionT& elemRegion = ielemRegion->second;

    // loop over all elements in the region
    for( localIndex k=0 ; k<elemRegion.m_numElems; ++k)
    {
      // get the elementToNodeMap for element k
      const localIndex* const elementToNodeMap = elemRegion.ElementToNodeMap(k);

      // loop over all nodes in elementToNodeMap
      for( unsigned int a=0 ; a<elemRegion.m_numNodesPerElem ; ++a )
      {
        // get local index of the node from elementToNodeMap
        const localIndex localNodeIndex = elementToNodeMap[a];

        // now set the NodeToElementMap as a combination of element region and element index.
        m_toElementsRelation(localNodeIndex).insert( std::make_pair(const_cast<ElementRegionT*>(&elemRegion),k));
      }
    }
  }
  // trim the array since the push_back function can allocate too much space.
//  m_NodeToElementMap.Trim();
}


void NodeManager::AddToNodeToElementMap( const ElementManagerT& elementManager,
                                          const std::map<std::string,lArray1d>& newElementIndices )
{
  // because the nodeToElementMap is an odd creature, it is not managed by ObjectDataStructureBaseT...so we must
  // resize.
  m_toElementsRelation.resize(DataLengths());

  // iterate over all element regions
  for( std::map<ElementManagerT::RegKeyType, ElementRegionT>::const_iterator ielemRegion=elementManager.m_ElementRegions.begin() ;
       ielemRegion!=elementManager.m_ElementRegions.end() ; ++ielemRegion )
  {
    // the element region is the mapped value of the iterator
    const std::string& regionName = ielemRegion->first;
    const ElementRegionT& elemRegion = ielemRegion->second;

    std::map<std::string,lArray1d>::const_iterator i=newElementIndices.find( regionName );
    if( i!= newElementIndices.end() )
    {
      const lArray1d& elementIndices = i->second;

      // loop over all elements in list
      for( size_t k=0 ; k<elementIndices.size(); ++k)
      {
        const localIndex elemIndex = elementIndices(k);

        // get the elementToNodeMap for element k
        const localIndex* const elementToNodeMap = elemRegion.ElementToNodeMap(elemIndex);

        // loop over all nodes in elementToNodeMap
        for( unsigned int a=0 ; a<elemRegion.m_numNodesPerElem ; ++a )
        {
          // get local index of the node from elementToNodeMap
          const localIndex localNodeIndex = elementToNodeMap[a];

          // now set the NodeToElementMap as a combination of element region and element index.
          m_toElementsRelation(localNodeIndex).insert( std::make_pair(const_cast<ElementRegionT*>(&elemRegion),elemIndex));
        }
      }
    }
  }
}




/**
 * @author R. Settgast
 * @param faceManager the face manager
 * this function constucts the nodeToFace map using the data in the face manager's facetonode map.
 */
void NodeManager::ConstructNodeToFaceMap( const FaceManagerT& faceManager )
{
  // loop over all faces
  for( localIndex lfi=0 ; lfi<faceManager.DataLengths() ; ++lfi )
  {
    // now iterate over the faceToNodeMap (i.e. all nodes in the faceToNodeMap)
    for( lArray1d::const_iterator a=faceManager.m_toNodesRelation[lfi].begin() ;
         a!=faceManager.m_toNodesRelation[lfi].end() ; ++a )
    {
      // enter the value of the face index into the nodeToFace map
      m_nodeToFaceMap[*a].insert(lfi);
    }
  }
}





void NodeManager::ModifyNodeToEdgeMapFromSplit( const EdgeManagerT& edgeManager,
                                                 const lSet& newEdgeIndices,
                                                 const lSet& modifiedEdgeIndices )
{

  lSet allEdges;
  allEdges.insert( newEdgeIndices.begin(), newEdgeIndices.end() );
  allEdges.insert( modifiedEdgeIndices.begin(), modifiedEdgeIndices.end() );

  // if an edge is new, we need to check who was connected to the parent, and modify accordingly
  for( lSet::const_iterator edgeIndex=allEdges.begin() ; edgeIndex!=allEdges.end() ; ++edgeIndex )
  {
    const localIndex* const nodelist = edgeManager.m_toNodesRelation[ *edgeIndex ];
    for( lArray2d::size_type a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a )
    {
      localIndex nodeIndex = nodelist[a];
      localIndex parentNodeIndex = m_parentIndex[nodeIndex];

      while( parentNodeIndex != LOCALINDEX_MAX )
      {
        m_nodeToEdgeMap[parentNodeIndex].erase(*edgeIndex);

        nodeIndex = parentNodeIndex;
        parentNodeIndex = m_parentIndex[nodeIndex];
      }
    }

    for( lArray2d::size_type a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a )
    {
      localIndex nodeIndex = nodelist[a];
      m_nodeToEdgeMap[nodeIndex].insert(*edgeIndex);
    }


  }

}

// Fu note on 20130416: Looks like this was temporary.  This function is now taken care of by element region. 
// We will keep it here for a while and delete it later.

//void NodeManagerT::UpdateNodeExternalityFromSplit( const FaceManagerT& faceManager,
//                                                 const lSet& newNodeIndices,
//                                                 const lSet& modifiedNodeIndices )
//{
//  lSet allNodes;
//  allNodes.insert( newNodeIndices.begin(), newNodeIndices.end() );
//  allNodes.insert( modifiedNodeIndices.begin(), modifiedNodeIndices.end() );
//
//
////  for( lSet::const_iterator nodeIndex=allNodes.begin() ; nodeIndex!=allNodes.end() ; ++nodeIndex )
//  for (localIndex nodeIndex = 0; nodeIndex != DataLengths(); ++nodeIndex)
//  {
//
//    for( lSet::const_iterator iface=m_nodeToFaceMap[nodeIndex].begin() ;
//        iface!=m_nodeToFaceMap[nodeIndex].end() ; ++iface )
//    {
//      if (faceManager.m_isExternal[*iface] == 1)
//      {
//        m_isExternal[nodeIndex] =1;
//      }
//    }
//  }
//}

void NodeManager::AddToNodeToFaceMap( const FaceManagerT& faceManager,
                                       const lArray1d& newFaceIndices )
{
  // loop over all faces in list
  for( size_t kf=0 ; kf<newFaceIndices.size() ; ++kf )
  {
    localIndex lfi = newFaceIndices(kf);
    // now iterate over the faceToNodeMap (i.e. all nodes in the faceToNodeMap)
    for( lArray1d::const_iterator a=faceManager.m_toNodesRelation[lfi].begin() ;
         a!=faceManager.m_toNodesRelation[lfi].end() ; ++a )
    {
      // enter the value of the face index into the nodeToFace map
      m_nodeToFaceMap[*a].insert(lfi);
    }
  }
}


void NodeManager::AddToNodeToEdgeMap( const EdgeManagerT& edgeManager,
                                       const lArray1d& newEdgeIndices )
{
  // loop over all edges in list
  for( size_t kf=0 ; kf<newEdgeIndices.size() ; ++kf )
  {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    localIndex lfi = newEdgeIndices(kf);

    for (localIndex i = 0; i < edgeManager.m_toNodesRelation.Dimension(1); ++i)
    {
      localIndex node = edgeManager.m_toNodesRelation(lfi,i);
      m_nodeToEdgeMap[node].insert(lfi);
    }
  }
}

/**
 * @brief Sort nodes on a plane.  We need this to sort node on a 2D element, which cannot be handled by the face manager.
 * @author Pengcheng Fu
 */
void NodeManager::SortNodeOnPlane (lArray1d& nodeList) const
{
  //GeometryUtilities::OrderPoints_2DPolygon();
  localIndex numNodes = nodeList.size();
  const Array1dT< R1Tensor >& refPosition = GetFieldData<FieldInfo::referencePosition>();
  const Array1dT< R1Tensor >& displacement = GetFieldData<FieldInfo::displacement>();

  if (numNodes<=2)
    throw GPException("You are trying to sort nodes on a plane but there are only two nodes.");

  R1Tensor fc;
  Array1dT<R1Tensor> nodeCoords(numNodes);
  for( localIndex i =0; i < numNodes; ++i)
  {
    localIndex nd = nodeList[i];
    nodeCoords[i] = refPosition[nd] + displacement[nd];
    fc += nodeCoords[i];
  }
  fc /= realT(numNodes);

  R1Tensor ex, ey, ez;

  {
    // Approximate face normal direction (unscaled)
    ex = nodeCoords[0];
    ex -= fc;
    ex /= ex.L2_Norm();

    ey = nodeCoords[1];
    ey -= fc;
    ey /= ey.L2_Norm();
    ez.Cross(ex, ey);
    if (ez.L2_Norm() < 0.01)  // Node 0, 1, and face center are roughly along one straight line.  We use another node to construct the vectors.
    {
      ey = nodeCoords[2]-fc;
      ey /= ey.L2_Norm();
    }

    ez.Cross(ex, ey); ez /= ez.L2_Norm();
    ey.Cross(ez,ex); ey /= ey.L2_Norm();
  }


  /// Sort nodes counterclockwise around face center
  {
    Array1dT< std::pair<realT,int> > thetaOrder(numNodes);
    for( unsigned int n =0; n < numNodes; ++n){
      R1Tensor v = nodeCoords[n]-fc;
      thetaOrder[n] = std::pair<realT,int>(atan2(v*ey,v*ex),nodeList[n]);
    }

    sort(thetaOrder.begin(), thetaOrder.end());

    // Reorder nodes on face
    for( unsigned int n =0; n < numNodes; ++n)
    {
      nodeList[n] = thetaOrder[n].second;
    }

    lArray1d tempNodes(numNodes);
    localIndex firstIndexIndex = 0;
    for( unsigned int n =0; n < numNodes; ++n)
    {
      tempNodes[n] = thetaOrder[n].second;
      if( tempNodes[n] == nodeList[0] )
      {
        firstIndexIndex = n;
      }
    }
    for( unsigned int n=0; n < numNodes; ++n)
    {
      const localIndex index = firstIndexIndex+n < numNodes ? firstIndexIndex+n : firstIndexIndex+n-numNodes;
      nodeList[n] = tempNodes[index];
    }


  }

}

void NodeManager::ZeroDetachedNodeVelocity()
{
  if (HasField<int>("isDetachedFromSolidMesh"))
  {
    iArray1d& isNodeDetached = GetFieldData<int>("isDetachedFromSolidMesh");
    Array1dT<R1Tensor>& velocity = GetFieldData<FieldInfo::velocity> ();

    for (localIndex i=0; i<DataLengths(); ++i)
    {
      if (isNodeDetached[i] == 1) velocity[i] *= 0;
    }

  }

}
void NodeManager::UpdateDetachedNodeLocationAndVelocity()
{
  iArray1d& isDetachedFromSolidMesh = GetFieldData<int> ("isDetachedFromSolidMesh");
  Array1dT<R1Tensor>& velocity = GetFieldData<FieldInfo::velocity> ();
  Array1dT<R1Tensor>& displacement = GetFieldData<FieldInfo::displacement> ();
  Array1dT<realT>& mass = GetFieldData<FieldInfo::mass> ();

  for (localIndex i=0; i < DataLengths(); ++i)
  {
    if (isDetachedFromSolidMesh[i] == 1)
    {
      lArray1d listOfEffectiveChildren;
      rArray1d weight;
      FindAllEffectiveChildren(i, listOfEffectiveChildren, weight);

      if (mass[listOfEffectiveChildren[0]] > 0)
      {
        realT totalWeight(0);
        for (localIndex j=0; j<listOfEffectiveChildren.size(); ++j)
        {
          weight[j] = mass[listOfEffectiveChildren[j]];
          totalWeight += weight[j];
        }
        weight /= totalWeight;
      }

      velocity[i] *= 0.0;
      displacement[i] *= 0.0;
      for (localIndex j=0; j<listOfEffectiveChildren.size(); ++j)
      {
        velocity[i] += velocity[listOfEffectiveChildren[j]] * weight[j];
        displacement[i] += displacement[listOfEffectiveChildren[j]] * weight[j];
      }
    }

  }

}

//Pre-order depth-first binary tree traversal
void NodeManager::FindAllEffectiveChildren(localIndex& nodeID,
                                            lArray1d& list,
                                            rArray1d& weight)
{
  lArray1d parentStack;
  iArray1d stackLevel;
  localIndex* node(&nodeID);
  int level(0);
  while (!parentStack.empty() || node != NULL)
  {
    if (node != NULL)
    {
      list.push_back(*node);
      weight.push_back(pow(0.5, level));
      parentStack.push_back(*node);
      stackLevel.push_back(level);
      if (m_childIndices[*node].size() > 1)
      {
        node = &(m_childIndices[*node][0]);
        ++level;
      }
      else
      {
        node = NULL;
      }
    }
    else
    {
      node = &(parentStack.back());
      level = stackLevel.back();
      parentStack.pop_back();
      stackLevel.pop_back();
      if (m_childIndices[*node].size() > 1)
      {
        node = &(m_childIndices[*node][1]);
        ++level;
      }
      else
      {
        node = NULL;
      }
    }
  }
  iArray1d& isDetachedFromSolidMesh = GetFieldData<int> ("isDetachedFromSolidMesh");

  rArray1d::iterator iw=weight.end()-1 ;
  for( lArray1d::iterator i=list.end()-1 ; i!=list.begin()-1 ; --i, --iw )
  {
    if (isDetachedFromSolidMesh[*i] == 1)
    {
      list.erase(i);
      weight.erase(iw);
    }
  }

}

/**
 * @author R. Settgast
 * @param sendnodes local indices of nodes to pack and send
 * @param buffer the buffer to pack the nodes into
 * @return size of characters packed into the buffer.
 *
 * This function packs complete nodes into a buffer. this should include all information needed to reconstruct the node
 * on a remote process domain. This does not include maps to other objects, as those are locally indexed relations and
 * must be constructed on the receiving domain.
 */
template< typename T_indices >
unsigned int NodeManager::PackNodes( const T_indices& sendnodes,
                                      const FaceManagerT& faceManager,
                                      bufvector& buffer,
                                      const bool packConnectivityToGlobal,
                                      const bool packFields,
                                      const bool packMaps,
                                      const bool packSets ) const
{
  unsigned int sizeOfPacked = 0;

  const std::string label = "NodeData";
  sizeOfPacked += buffer.Pack(label);

  sizeOfPacked += ObjectDataStructureBaseT::PackBaseObjectData( buffer, sendnodes, packFields, packMaps, packSets, packConnectivityToGlobal );




  const Array1dT<lSet>* const usedFacesPointer = GetUnorderedVariableOneToManyMapPointer("usedFaces");

  if( usedFacesPointer!=NULL )
  {
    const Array1dT<lSet>& usedFaces = *usedFacesPointer;
    for( typename T_indices::const_iterator iter_index=sendnodes.begin() ; iter_index!=sendnodes.end() ; ++iter_index )
    {
      if( packConnectivityToGlobal )
      {
        sizeOfPacked += buffer.PackGlobal( usedFaces[*iter_index], faceManager.m_localToGlobalMap );
      }
      else
      {
        sizeOfPacked += buffer.Pack( usedFaces[*iter_index] );
      }
    }
  }

  return sizeOfPacked;
}
template unsigned int NodeManager::PackNodes<lSet>(  const lSet&, const FaceManagerT&, bufvector&, const bool, const bool, const bool, const bool  ) const;
template unsigned int NodeManager::PackNodes<lArray1d>( const lArray1d&, const FaceManagerT&, bufvector&, const bool, const bool, const bool, const bool  ) const;


/**
 * @author R. Settgast
 * @param[in,out] buffer the buffer to unpack nodes from
 * @param[in,out] nodeReceiveLocalIndices the local indices of the nodes
 * @return
 */
unsigned int NodeManager::UnpackNodes( const char*& buffer,
                                        const FaceManagerT& faceManager,
                                        lArray1d& nodeReceiveLocalIndices,
                                        const bool unpackConnectivityToLocal,
                                        const bool unpackFields,
                                        const bool unpackMaps,
                                        const bool unpackSets  )
{
  unsigned int sizeOfUnpacked = 0;

  const std::string label = "NodeData";
  std::string temp;

  sizeOfUnpacked += bufvector::Unpack( buffer, temp );
  if( label.compare(temp)!=0 )
  {
    throw GPException("NodeManagerT::UnpackNodes: buffer location incorrect\n");
  }

  lArray1d junk;
  sizeOfUnpacked += ObjectDataStructureBaseT::UnpackBaseObjectData( buffer, nodeReceiveLocalIndices, junk, unpackFields, unpackMaps, unpackSets, unpackConnectivityToLocal );




  Array1dT<lSet>* const usedFacesPointer = GetUnorderedVariableOneToManyMapPointer("usedFaces");

  if( usedFacesPointer!=NULL )
  {
    Array1dT<lSet>& usedFaces = *usedFacesPointer;
    for( lArray1d::const_iterator a=nodeReceiveLocalIndices.begin() ; a!=nodeReceiveLocalIndices.end() ; ++a )
    {
      if( unpackConnectivityToLocal )
      {
        sizeOfUnpacked += bufvector::UnpackGlobal( buffer, faceManager.m_globalToLocalMap,
                                                   usedFaces[*a] );
      }
      else
      {
        sizeOfUnpacked += bufvector::Unpack( buffer,
                                             usedFaces[*a] );
      }
    }
  }



  return sizeOfUnpacked;

}


void NodeManager::ConnectivityFromGlobalToLocal( const lSet& indices,
                                                  const lSet& clearIndices,
                                                  const std::map<globalIndex,localIndex>& faceGlobalToLocal )
{
  Array1dT<lSet>* const usedFacesPointer = GetUnorderedVariableOneToManyMapPointer("usedFaces");
  if( usedFacesPointer!=NULL )
  {
    Array1dT<lSet>& usedFaces = *usedFacesPointer;

    for( lSet::const_iterator a=indices.begin() ; a!=indices.end() ; ++a )
    {
      gArray1d temp;

      for( lSet::iterator faceIndex=usedFaces[*a].begin() ; faceIndex!=usedFaces[*a].end() ; ++faceIndex  )
      {
        temp.push_back(stlMapLookup( faceGlobalToLocal, static_cast<globalIndex>(*faceIndex) ) );
      }
      usedFaces[*a].clear();
      for( gArray1d::const_iterator faceIndex=temp.begin() ; faceIndex!=temp.end() ; ++faceIndex  )
      {
        usedFaces[*a].insert( *faceIndex );
      }
    }


    for( lSet::const_iterator a=clearIndices.begin() ; a!=clearIndices.end() ; ++a )
    {
      usedFaces[*a].clear();
    }
  }
}

void NodeManager::WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                                     const std::string&,
                                                     const std::string&,
                                                     const int,
                                                     const int,
                                                     const realT,
                                                     const bool isRestart,
                                                     const std::string&,
                                                     const std::string&,
                                                     const lArray1d&)
{
  if( isRestart )
  {
    siloFile.DBWriteWrapper("m_matchedBoundaryNodes",m_matchedBoundaryNodes);
  }



}

void NodeManager::ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                                      const std::string&,
                                                      const std::string&,
                                                      const int,
                                                      const int,
                                                      const realT,
                                                      const bool isRestart,
                                                      const std::string&,
                                                      const lArray1d&)
{
  if( isRestart )
  {
    siloFile.DBReadWrapper("m_matchedBoundaryNodes",m_matchedBoundaryNodes);

  }



}

void NodeManager::CalculateEffectiveNormal( const localIndex index,
                                             const FaceManagerT& faceManager,
                                             R1Tensor& normal ) const
{
  const lSet& faceList = this->m_nodeToFaceMap[index];

  const iArray1d& isFaceExternal = faceManager.GetFieldData<int>("isExternal");

  normal = 0.0;

  for( lSet::const_iterator kf=faceList.begin() ; kf!=faceList.end() ; ++kf )
  {
    if( isFaceExternal[*kf] > 0 )
    {
      R1Tensor N = faceManager.FaceNormal( *this, *kf );
      normal += N;
    }
  }
  normal.Normalize();
}


void NodeManager::GetDomainExtents(R1Tensor& pmin, R1Tensor& pmax, localIndex Ndims)
{
  const Array1dT<R1Tensor>& referencePosition = this->GetFieldData<FieldInfo::referencePosition>();

  pmin = referencePosition[0];
  pmax = referencePosition[0];
  for (localIndex ii=0; ii<referencePosition.size(); ii++){
    for (localIndex dim=0; dim<Ndims; dim++){
      pmin[dim] = (pmin[dim]<referencePosition[ii][dim])?pmin[dim]:referencePosition[ii][dim];
      pmax[dim] = (pmax[dim]>referencePosition[ii][dim])?pmax[dim]:referencePosition[ii][dim];
    }  
  }
}
