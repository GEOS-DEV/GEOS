/**
 * @file FaceManager.cpp
 * @author settgast1
 */

#include "FaceManager.hpp"

namespace geosx
{

/**
 *
 * @return
 */
FaceManager::FaceManager( string const & , ObjectManagerBase * const parent ):
ObjectManagerBase("FaceManager",parent),
{
  //0-based; note that the following field is ALSO 0
  //for faces that are not external faces, so check isExternal before using
  this->AddKeylessDataField<localIndex>("externalFaceIndex", true, true);

  this->AddKeylessDataField<R1Tensor>("FaceCenter",true,true);
}

/**
 *
 * @return
 */
FaceManager::~FaceManager()
{
#if USECPP11!=1
  if( m_cohesiveZone )
    delete m_cohesiveZone;
#endif

}
//
//void FaceManager::BuildFaces( const NodeManager& nodeManager, const ElementManager& elementManager )
//{
//
//  lArray1d tempNodeList;
//  Array1dT<lArray1d> tempFaceToNodeMap;
//
//  localIndex numFaces = 0;
//  Array1dT<lArray1d> facesByLowestNode;
//
//
//
//  for( std::map< ElementManager::RegKeyType, ElementRegionT >::const_iterator elementRegionIter = elementManager.m_ElementRegions.begin() ;
//       elementRegionIter != elementManager.m_ElementRegions.end() ;
//       ++elementRegionIter )
//  {
//    const ElementRegionT& elementRegion = elementRegionIter->second;
//
////    if( elementRegion.m_ElementDimension == 3 )
//    {
//      for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
//      {
//        // kelf = k'th element local face index
//        for( int kelf=0 ; kelf<elementRegion.m_numFacesPerElement ; ++kelf )
//        {
//          // get the nodes associated with the local face
//          elementRegion.GetFaceNodes( k, kelf, tempNodeList );
//
//          //Special treatment for the triangle faces of prisms.
//          if (tempNodeList[tempNodeList.size()-1] == std::numeric_limits<localIndex>::max()) tempNodeList.pop_back();
//
//          // sort the nodes
//          std::sort(tempNodeList.begin(), tempNodeList.end() );
//
//          // get the lowest node index from the list for simplicity
//          const localIndex& lowNode = tempNodeList[0];
//
//          // now check to see if the lowest node index has an entry in the facesByLowestNode vector
//          if( facesByLowestNode.size() < (lowNode+1) )
//          {
//            // the node has not been entered, so add it.
//            facesByLowestNode.resize(lowNode+1);
//
//            // this a new face, so add it,
//            AddNewFace( k, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, elementRegion );
//          }
//          else
//          {
//            // does the node have an entry? If not, then this has to be a new face.
//            if( facesByLowestNode[lowNode].empty() )
//            {
//              // this a new face, so add it,
//              AddNewFace( k, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, elementRegion );
//            }
//            else
//            {
//              // the node does have an entry, so it is possible that the facet has already be assigned a number
//
//              // make a flag to indicate whether the face is a duplicate...assume that it isn't unless this is disproved.
//              bool duplicate = false;
//
//
//              // there are faces in facesByLowestNode, so lets loop over them and check for duplicates
//              for( lvector::iterator existingFaceIndex = facesByLowestNode[lowNode].begin() ;
//                  existingFaceIndex != facesByLowestNode[lowNode].end() ; ++existingFaceIndex )
//              {
//                // this is the nodelist of the face that we are testing agains
//                const lArray1d& existingFaceNodelist = tempFaceToNodeMap[*existingFaceIndex];
//
//                // test to see if the size of the nodelists are the same....
//                if( existingFaceNodelist.size() == tempNodeList.size() )
//                {
//                  // since the size is the same, then we should test the nodes...they are sorted, so
//                  // the std::equal() algorithm will work for this.
//                  if( std::equal( existingFaceNodelist.begin(), existingFaceNodelist.end(), tempNodeList.begin() ) )
//                  {
//                    // they are equal!
//                    duplicate = true;
//
//                    // add the element to the faceToElement map
//                    m_toElementsRelation[*existingFaceIndex].push_back( std::pair<ElementRegionT*, localIndex>( const_cast<ElementRegionT*>(&elementRegion), k) );
//
//                    // add the face to the elementToFaceMap for the element region.
//                    elementRegion.m_toFacesRelation(k,kelf) = *existingFaceIndex;
//
//                    // now remove the entry from the face that we were checking against from the facesByLowestNode list...
//                    // because it is no longer possible that it will have another element that has this face.
//		    facesByLowestNode[lowNode].erase( existingFaceIndex );
//
//                    // break the loop
//                    break;
//                  }
//                }
//              }
//              if( !duplicate )
//              {
//                // the face is not a duplicate of any in the facesByLowestNode list, so we need to add a new face.
//                AddNewFace( k, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, elementRegion );
//              }
//            }
//          }
//        }
//      }
//    }
//  }
//
//  // resize the data vectors according to the number of faces indicated by the size of tempFaceToNodeMap
//  this->resize(tempFaceToNodeMap.size());
//
//  // set m_FaceToNodeMap
//  this->m_toNodesRelation = tempFaceToNodeMap;
//
//  auto const & nodeSets = nodeManager.GetGroup(string("Sets")).wrappers();
//
//  // make sets from nodesets
//  for( auto const & setWrapper : nodeSets )
//  {
//    const std::string& setname = setWrapper->name();
//    const lSet& set = ( dataRepository::ViewWrapper<lSet>::cast( *setWrapper ) ).reference() ;
//    this->ConstructSetFromSetAndMap( set, this->m_toNodesRelation, setname );
//  }
//
//
//  // sort the face node lists
//  SortAllFaceNodes(nodeManager);
//
//
//  Array1dT<R1Tensor>& faceCenter = this->GetFieldData<R1Tensor>( "FaceCenter" );
//  for( localIndex k=0 ; k<DataLengths() ; ++k )
//  {
//    FaceCenter( nodeManager, k , faceCenter[k] );
//
//  }
//
//
//  // Figure out if we need to write arbitrary polygons to silo.  Have to do this here because cannot do allreduce in the write silo file.
//  int maxNodePerFace(-100), minNodePerFace(1000), writeArbitraryPolygonLocal(0);
//  for (localIndex kf = 0; kf < DataLengths(); ++kf)
//  {
//    maxNodePerFace = std::max(maxNodePerFace, int(m_toNodesRelation[kf].size()));
//    minNodePerFace = std::min(minNodePerFace, int(m_toNodesRelation[kf].size()));
//  }
//  if (maxNodePerFace != minNodePerFace || maxNodePerFace > 4) writeArbitraryPolygonLocal = 1;
//  MPI_Allreduce(&writeArbitraryPolygonLocal, &m_writeArbitraryPolygon, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//
//}



//
//void FaceManager::AddNewFace( const localIndex& k,
//                               const localIndex& kelf,
//                               localIndex& numFaces,
//                               Array1dT<lArray1d>& facesByLowestNode,
//                               lArray1d& tempNodeList,
//                               Array1dT<lArray1d>& tempFaceToNodeMap,
//                               const ElementRegionT& elementRegion )
//{
//  Array1dT< std::pair< ElementRegionT*, localIndex > > tempFaceToElemEntry;
//
//  // and add the face to facesByLowestNode[]
//  facesByLowestNode[tempNodeList[0]].push_back(numFaces);
//
//
//  // add the face to the elementToFaceMap
//  elementRegion.m_toFacesRelation(k,kelf) = numFaces;
//
//  // add the nodes to the faceToNodeMap
//  tempFaceToNodeMap.push_back(tempNodeList);
//
//  // add to the element information to the faceToElementMap
//  tempFaceToElemEntry.push_back( std::pair<ElementRegionT*, localIndex>( const_cast<ElementRegionT*>(&elementRegion), k) );
//  m_toElementsRelation.push_back(tempFaceToElemEntry);
//
//  // now increment numFaces to reflect the number of faces rather than the index of the new face
//  ++numFaces;
//}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceManager, std::string const &, ObjectManagerBase * const )

}

