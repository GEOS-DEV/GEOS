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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * ElementManagerT.cpp
 *
 *  Created on: Sep 14, 2010
 *      Author: settgast1
 */

#include "ElementManagerT.h"
#include "FaceManagerT.h"
#include "IO/BinStream.h"
#include <map>
#include <vector>
#include "Constitutive/Material/MaterialFactory.h"



ElementManagerT::ElementManagerT():
  ObjectDataStructureBaseT( ObjectDataStructureBaseT::ElementManager),
  m_numElems(0),
  m_ElementToNodeMap(m_VariableOneToManyMaps["ElementToNode"]),
  m_ElementToFaceMap(m_VariableOneToManyMaps["ElementToFace"]),
  m_ElementToElementMap(m_VariableOneToManyMaps["ElementToElement"]),
  m_ElementIndexToRegionIndex(m_IntegerData["ElementIndexToRegionIndex"]),
  m_ElementIndexToRegionLocalIndex(m_OneToOneMaps["ElementIndexToRegionLocalIndex"]),
  m_ElementRegions()
{
  FieldInfo::AttributesByName["ElementIndexToRegionIndex"]  = new FieldBase( FieldInfo::noKey, "ElementIndexToRegionIndex", true, false );
//  std::cout<<"ElementManagerT constructor
// "<<m_FixedOneToManyMaps.size()<<std::endl;

}

ElementManagerT::~ElementManagerT()
{
  // TODO Auto-generated destructor stub
}


globalIndex ElementManagerT::resize( const lvector& numElements,
                                     const array<string>& elementRegionNames,
                                     const array<string>& elementTypes )
{
  m_numElems = 0;

  if( numElements.size() != elementRegionNames.size() || numElements.size() != elementTypes.size() )
  {
    throw GPException( "ElementManagerT::resize(): size mismatch");
  }

  std::map< std::string, std::pair< std::string,localIndex> > regionData;
  const lvector::size_type size=numElements.size();
  for( lvector::size_type k=0 ; k<size ; ++k )
  {
    localIndex nelem = 0;
    if( regionData.find( elementRegionNames[k] ) != regionData.end() )
      nelem = regionData[ elementRegionNames[k] ].second;
    regionData[ elementRegionNames[k] ] = std::make_pair( elementTypes[k], nelem + numElements[k] );
  }


  for( std::map< std::string, std::pair< std::string,localIndex> >::iterator iterRegion=regionData.begin() ; iterRegion!=regionData.end() ; ++iterRegion )
  {
    const std::string regionName = iterRegion->first;
    const std::string elemType = iterRegion->second.first;
    const localIndex numElems = iterRegion->second.second;
    ElementRegionT* const elemRegion = stlMapLookupPointer(m_ElementRegions,regionName);

    if( elemRegion==NULL )
    {
      m_regionIndexToKey.push_back(regionName);
      const size_t regionNumber = m_ElementRegions.size();
      m_ElementRegions.insert( std::pair<std::string,ElementRegionT>(regionName,ElementRegionT()) );
      m_ElementRegions[regionName].m_regionName = regionName;
      m_ElementRegions[regionName].m_regionNumber = regionNumber;
    }
    m_ElementRegions[regionName].m_elementGeometryID = elemType;
    m_ElementRegions[regionName].SetGeometryBasedVariables();
    m_ElementRegions[regionName].resize(numElems);


  }

  // Fu: Calculate the total number of element separately will allow only
  // resizing one region without passing the other regions in.
  for( std::map< RegKeyType, ElementRegionT >::iterator iter=m_ElementRegions.begin() ; iter!=m_ElementRegions.end() ; ++iter)
  {
    m_numElems += iter->second.m_numElems;
  }


  const globalIndex firstNewGlobalIndex = ObjectDataStructureBaseT::resize(m_numElems, false);

  return firstNewGlobalIndex;
}

void ElementManagerT::SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject )
{
  for( std::map< RegKeyType, ElementRegionT >::iterator elementRegion=m_ElementRegions.begin() ;
       elementRegion!=m_ElementRegions.end() ; ++elementRegion )
  {
    elementRegion->second.SetDomainBoundaryObjects(referenceObject);
  }
}

void ElementManagerT::ResetGlobalToLocalMap( )
{
  for( std::map< RegKeyType, ElementRegionT >::iterator elementRegion=m_ElementRegions.begin() ;
       elementRegion!=m_ElementRegions.end() ; ++elementRegion )
  {
    elementRegion->second.ResetGlobalToLocalMap();
  }
}

/*
   void ElementManagerT::ConstructListOfBoundaryObjects( gArray1d& objectList )
      const
   {
   for( std::map< RegKeyType, ElementRegionT >::const_iterator
      elementRegion=m_ElementRegions.begin() ;
      elementRegion!=m_ElementRegions.end(); ++elementRegion )
   {
   //    elementRegion->second.ConstructListOfBoundaryObjects(objectList);
   }
   }
 */

void ElementManagerT::HACKInitialConditions(  )
{
  for( std::map< RegKeyType, ElementRegionT >::iterator i=m_ElementRegions.begin() ;
       i!=m_ElementRegions.end() ; ++i )
  {
    ElementRegionT& elementRegion = i->second;


    const array<real64>& sigma_x = elementRegion.GetFieldData<realT>("sigma_x");
    const array<real64>& sigma_y = elementRegion.GetFieldData<realT>("sigma_y");
    const array<real64>& sigma_z = elementRegion.GetFieldData<realT>("sigma_z");
    const array<real64>& sigma_xy = elementRegion.GetFieldData<realT>("sigma_xy");
    const array<real64>& sigma_yz = elementRegion.GetFieldData<realT>("sigma_yz");
    const array<real64>& sigma_xz = elementRegion.GetFieldData<realT>("sigma_xz");

    array<real64>& pressure       = elementRegion.GetFieldData<FieldInfo::pressure>();
    array<R2SymTensor>& s = elementRegion.GetFieldData<FieldInfo::deviatorStress>();

    array<R2SymTensor> *refStress = elementRegion.GetFieldDataPointer<R2SymTensor>("referenceStress");

    if (refStress != nullptr)
    {
      for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
      {
        (*refStress)[k](0,0) = sigma_x[k];
        (*refStress)[k](1,1) = sigma_y[k];
        (*refStress)[k](2,2) = sigma_z[k];
        (*refStress)[k](0,1) = sigma_xy[k];
        (*refStress)[k](1,2) = sigma_yz[k];
        (*refStress)[k](0,2) = sigma_xz[k];
      }
    }

    for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
    {

      pressure[k] = ( sigma_x[k] + sigma_y[k] + sigma_z[k] ) / 3.0;

      s[k](0,0) = sigma_x[k] - pressure[k];
      s[k](1,1) = sigma_y[k] - pressure[k];
      s[k](2,2) = sigma_z[k] - pressure[k];
      s[k](0,1) = sigma_xy[k];
      s[k](1,2) = sigma_yz[k];
      s[k](0,2) = sigma_xz[k];

      for( localIndex a=0 ; a<elementRegion.m_numIntegrationPointsPerElem ; ++a )
      {
        elementRegion.m_mat->StateData(k,a)->pressure = pressure[k];
        elementRegion.m_mat->StateData(k,a)->devStress = s[k];
      }
    }
  }
}


/**
 * @author settgast
 * @param toElementMap
 * @param nodeList
 * @param localIndexes
 * @param depth
 *
 * This function takes a list of nodes, and makes a list of all elements
 * connected to those nodes.
 *
 */
void ElementManagerT::ConstructListOfIndexesFromMap( const array< std::set< std::pair<ElementRegionT*,localIndex> > >& toElementMap,
                                                     const lArray1d& nodeList,
                                                     std::map< std::string, lArray1d>& localIndexes,
                                                     const int depth )
{


  // list to hold all attached nodes for depth.
  lSet nodesInDepth( nodeList.begin(), nodeList.end() );



  // we need to all elements at a "depth" (i.e. connectivity distance) away from
  // the list nodes. So we need a set of nodes
  // that are a certain connectivity distance from the original nodeList.

  // get all nodes that are in the depth range
  for( int d=1 ; d<depth ; ++d )
  {
    lSet oldSet(nodesInDepth);

    // loop over all nodes in set
    for( lSet::const_iterator nodeIndex=oldSet.begin() ; nodeIndex!=oldSet.end() ; ++nodeIndex )
    {

      // now that we have the node, we look for every element attached to it
      for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator elementPair=toElementMap[*nodeIndex].begin() ;
           elementPair!=toElementMap[*nodeIndex].end() ;
           ++elementPair )
      {
        const ElementRegionT& elemRegion = *(elementPair->first);
        const localIndex elemIndex = elementPair->second;

        // now we add all nodes connected to that element to the nodesInDepth
        for( unsigned int a=0 ; a<elemRegion.m_numNodesPerElem ; ++a )
        {
          nodesInDepth.insert( elemRegion.ElementToNodeMap(elemIndex,a));
        }
      }
    }
  }



  for( std::map< RegKeyType, ElementRegionT >::const_iterator i=m_ElementRegions.begin() ; i!=m_ElementRegions.end() ; ++i )
  {
    localIndexes[i->first];
  }

  // loop over all nodes in list
  for( lSet::const_iterator nodeIndex=nodesInDepth.begin() ; nodeIndex!=nodesInDepth.end() ; ++nodeIndex )
  {

    // now that we have the node, we look for every element attached to it
    for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator elementPair=toElementMap[*nodeIndex].begin() ;
         elementPair!=toElementMap[*nodeIndex].end() ;
         ++elementPair )
    {
      const ElementRegionT& elemRegion = *(elementPair->first);
      const localIndex elemIndex = elementPair->second;
      // add the element to the map
      localIndexes[elemRegion.m_regionName].push_back( elemIndex );
    }
  }


  // remove duplicate entries
  for( std::map< std::string, lArray1d>::iterator localPair=localIndexes.begin() ; localPair!=localIndexes.end() ; ++localPair )
  {
    lArray1d& objectList = localPair->second;

    // sort the entries
    std::sort(objectList.begin(),objectList.end());
    // now remove the duplicates
    lArray1d::iterator iend = std::unique(objectList.begin(),objectList.end());
    objectList.resize( iend - objectList.begin() );


  }
}

template< typename T_indices >
unsigned int ElementManagerT::PackElements( bufvector& buffer,
                                            lSet& sendnodes,
                                            lSet& sendfaces,
                                            const std::map<std::string,T_indices>& elementList,
                                            const NodeManager& nodeManager,
                                            const FaceManagerT& faceManager,
                                            const bool packConnectivityToGlobal,
                                            const bool packFields,
                                            const bool packMaps,
                                            const bool packSets ) const
{
  unsigned int sizeOfPacked = 0;


  int numRegions = elementList.size();

  sizeOfPacked += buffer.Pack( numRegions );

  for( typename std::map<std::string,T_indices>::const_iterator ilist=elementList.begin() ; ilist!=elementList.end() ; ++ilist )
  {
    const std::string& regionName = ilist->first;
    const T_indices& list = ilist->second;

    const ElementRegionT& elemRegion = stlMapLookup( m_ElementRegions, regionName, "ElementManagerT::PackElements()" );


    sizeOfPacked += buffer.Pack(regionName);



    sizeOfPacked += elemRegion.PackElements( buffer, sendnodes, sendfaces, list, nodeManager, faceManager, packConnectivityToGlobal, packFields, packMaps,
                                             packSets );
  }

  return sizeOfPacked;

}
template unsigned int ElementManagerT::PackElements( bufvector& buffer, lSet& sendnodes, lSet& sendfaces, const std::map<std::string,lSet>& elementList,
                                                     const NodeManager& nodeManager, const FaceManagerT& faceManager, const bool, const bool, const bool,
                                                     const bool ) const;
template unsigned int ElementManagerT::PackElements( bufvector& buffer, lSet& sendnodes, lSet& sendfaces, const std::map<std::string,lArray1d>& elementList,
                                                     const NodeManager& nodeManager, const FaceManagerT& faceManager, const bool, const bool, const bool,
                                                     const bool ) const;


unsigned int ElementManagerT::UnpackElements( const bufvector& buffer,
                                              const NodeManager& nodeManager,
                                              const FaceManagerT& faceManager,
                                              std::map< std::string, lArray1d>& elementRegionReceiveLocalIndices,
                                              const bool unpackConnectivityToLocal,
                                              const bool unpackFields,
                                              const bool unpackMaps,
                                              const bool unpackSets )
{
  const char* pbuffer = buffer.data();


  return UnpackElements( pbuffer, nodeManager, faceManager, elementRegionReceiveLocalIndices, unpackConnectivityToLocal, unpackFields, unpackMaps, unpackSets );

}

unsigned int ElementManagerT::UnpackElements( const char*& pbuffer,
                                              const NodeManager& nodeManager,
                                              const FaceManagerT& faceManager,
                                              std::map< std::string, lArray1d>& elementRegionReceiveLocalIndices,
                                              const bool unpackConnectivityToLocal,
                                              const bool unpackFields,
                                              const bool unpackMaps,
                                              const bool unpackSets  )
{
  unsigned int sizeOfUnpacked = 0;

  int numRegions = 0;

  sizeOfUnpacked += bufvector::Unpack( pbuffer, numRegions );

  for( int r=0 ; r<numRegions ; ++r )
  {
    std::string regionName;
    sizeOfUnpacked += bufvector::Unpack( pbuffer, regionName );

    ElementRegionT& elementRegion = this->m_ElementRegions[regionName];

    elementRegion.UnpackElements( pbuffer, nodeManager, faceManager,
                                  elementRegionReceiveLocalIndices[regionName],
                                  unpackConnectivityToLocal, unpackFields, unpackMaps, unpackSets );
  }

  m_numElems = 0;
  for( std::map< RegKeyType, ElementRegionT >::const_iterator i=m_ElementRegions.begin() ;
       i!=m_ElementRegions.end() ; ++i )
  {
    const ElementRegionT& elementRegion = i->second;
    m_numElems += elementRegion.m_numElems;

  }

  return sizeOfUnpacked;

}



void ElementManagerT::ConnectivityFromGlobalToLocal( const std::map< std::string, lSet>& allReceivedElements,
                                                     const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                                     const std::map<globalIndex,localIndex>& faceGlobalToLocal )
{

  for( std::map<std::string,lSet>::const_iterator iter=allReceivedElements.begin() ; iter!=allReceivedElements.end() ; ++iter )
  {
    const std::string& regionName = iter->first;
    const lSet& list = iter->second;

    ElementRegionT& elemRegion = stlMapLookup( m_ElementRegions, regionName );

    elemRegion.ConnectivityFromGlobalToLocal( list, nodeGlobalToLocal, faceGlobalToLocal );
  }

}


void ElementManagerT::ModifyToElementMapsFromSplit( const std::map< std::string, lSet>& modifiedElements,
                                                    NodeManager& nodeManager,
                                                    FaceManagerT& faceManager )
{
  for( std::map< RegKeyType, ElementRegionT >::iterator i=m_ElementRegions.begin() ; i!=m_ElementRegions.end() ; ++i )
  {
    const std::string& regionName = i->first;
    ElementRegionT& elementRegion = i->second;

    const std::map< std::string, lSet>::const_iterator iter = modifiedElements.find( regionName );

    if( iter!=modifiedElements.end() )
    {
      const lSet& modifiedElementSet = iter->second;
      elementRegion.ModifyToElementMapsFromSplit( modifiedElementSet, nodeManager, faceManager );
    }
  }
}

void ElementManagerT::UpdateExternalityFromSplit( const std::map< std::string, lSet>& modifiedElements,
                                                  NodeManager& nodeManager,
                                                  EdgeManagerT& edgeManager,
                                                  FaceManagerT& faceManager )
{
  for( std::map< RegKeyType, ElementRegionT >::iterator i=m_ElementRegions.begin() ; i!=m_ElementRegions.end() ; ++i )
  {
    const std::string& regionName = i->first;
    ElementRegionT& elementRegion = i->second;

    const std::map< std::string, lSet>::const_iterator iter = modifiedElements.find( regionName );

    if( iter!=modifiedElements.end() )
    {
      const lSet& modifiedElementSet = iter->second;
      elementRegion.UpdateExternalityFromSplit( modifiedElementSet, nodeManager, edgeManager, faceManager );
    }
  }
}



void ElementManagerT::InitializeFlowFaceRegion()
{
  lvector numElements;
  array<string> elementRegionNames;
  array<string> elementTypes;
  std::string type0=" ";
  localIndex nRegion = 0;
  for( std::map< RegKeyType, ElementRegionT >::iterator iter=m_ElementRegions.begin() ; iter!=m_ElementRegions.end() ; ++iter)
  {
    nRegion++;
//    numElements.push_back(0);
//    elementRegionNames.push_back(iter->first);
//    elementTypes.push_back( iter->second.m_elementGeometryID);

    if ( nRegion > 1 && iter->second.m_elementGeometryID != type0 && (type0 == "C3D4" || type0 == "C3D8"))
      throw GPException( "Element types in multiple 3D solid regions are inconsistent.  This is not allowed");
    // For 3D problem, the element regions must have the same type; we can mix
    // 2D regions
    type0 = iter->second.m_elementGeometryID;
  }
  numElements.push_back(0);
  elementRegionNames.push_back("FlowFaceRegion");

  if (type0 == "C3D4")
  {
    elementTypes.push_back("TRSH");
  }
  else if (type0 == "C3D8")
  {
    elementTypes.push_back("S4R");
  }
  else
  {
    throw GPException( "Flow faces other than triangles and quads are currently not implemented. Flow face regions are current not supported for prisms.");
  }

  resize( numElements, elementRegionNames, elementTypes );
}

void ElementManagerT::GenerateFlowFaceRegion(FaceManagerT& faceManager)
{
  lvector numElements;
  array<string> elementRegionNames;
  array<string> elementTypes;

  ElementRegionT& elemRegion = m_ElementRegions["FlowFaceRegion"];
  lSet facesToInclude;

  if( elemRegion.m_parentFaceSetNames.empty() )
  {
    for (localIndex i = 0 ; i < faceManager.DataLengths() ; ++i)
      facesToInclude.insert(i);
  }
  else
  {
    for( array<string>::size_type i =0 ; i < elemRegion.m_parentFaceSetNames.size() ; ++i)
    {
      lSet& set = faceManager.GetSet(elemRegion.m_parentFaceSetNames[i]);
      facesToInclude.insert(set.begin(), set.end());
    }
  }



  numElements.push_back(facesToInclude.size());
  elementRegionNames.push_back("FlowFaceRegion");
  elementTypes.push_back(elemRegion.m_elementGeometryID);
  resize( numElements, elementRegionNames, elementTypes );


  localIndex count = 0;

  for (lSet::const_iterator it = facesToInclude.begin() ; it != facesToInclude.end() ; ++it)
  {
    if (faceManager.m_toNodesRelation[*it].size() != elemRegion.m_numNodesPerElem)
      throw GPException( "The number of nodes on face is inconsistent with the number of nodes per FlowFace element");
    for (localIndex i = 0 ; i < elemRegion.m_numNodesPerElem ; ++i)
    {
      elemRegion.m_toNodesRelation[count][i] = faceManager.m_toNodesRelation[*it][i];
    }
    elemRegion.m_localToGlobalMap[count] = faceManager.m_localToGlobalMap[*it];
    elemRegion.m_toFacesRelation[count][0] = *it;

    // Add the new element to the face to element map of the face itself
    std::pair< ElementRegionT*, localIndex > tempFaceToElemEntry;
    tempFaceToElemEntry = std::make_pair(  &elemRegion, count);
    faceManager.m_toElementsRelation[*it].push_back(tempFaceToElemEntry);
    ++count;
  }
}

void ElementManagerT::WriteSilo( SiloFile& siloFile,
                                 const std::string& meshname,
                                 const int cycleNum,
                                 const realT problemTime,
                                 const bool isRestart,
                                 const std::string& regionName,
                                 const lArray1d& mask )
{

  siloFile.DBWriteWrapper("m_numElems",m_numElems );


  siloFile.WriteRegionSpecifications(*this, meshname, cycleNum, problemTime);

  array<string> regionNames;

  for (std::map<RegKeyType, ElementRegionT>::iterator elementRegionIter=m_ElementRegions.begin() ;
       elementRegionIter!= m_ElementRegions.end() ; ++elementRegionIter)
  {
    const std::string& elementRegionName = elementRegionIter->first;
    ElementRegionT& elementRegion = elementRegionIter->second;

    regionNames.push_back(elementRegionName);

//    if (elementRegion.m_numElems > 0)
    {

      elementRegion.WriteSiloRegionMesh( siloFile, meshname, cycleNum, problemTime, isRestart,
                                         elementRegionName);

    }
  }

  siloFile.DBWriteWrapper( "regionNames", regionNames);
}

void ElementManagerT::ReadSilo( const SiloFile& siloFile,
                                const std::string& meshname,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart )
{

  siloFile.DBReadWrapper("m_numElems",m_numElems );

  array<string> regionNames;
  siloFile.DBReadWrapper( "regionNames", regionNames );

  for( array<string>::const_iterator i=regionNames.begin() ; i!=regionNames.end() ; ++i )
  {
    m_ElementRegions[*i];
  }


  for (std::map<RegKeyType, ElementRegionT>::iterator elementRegionIter=m_ElementRegions.begin() ;
       elementRegionIter!= m_ElementRegions.end() ; ++elementRegionIter)
  {
    const std::string& elementRegionName = elementRegionIter->first;
    ElementRegionT& elementRegion = elementRegionIter->second;

    elementRegion.ReadSiloRegionMesh( siloFile, meshname, cycleNum, problemTime, isRestart,
                                      elementRegionName);

  }
}

void ElementManagerT::ResizeNumberOfElementsAfterSplit ()
{
  this->m_numElems = 0;
  for (auto& i:m_ElementRegions)
  {
    this->m_numElems = this->m_numElems + i.second.m_numElems;
  }
}
