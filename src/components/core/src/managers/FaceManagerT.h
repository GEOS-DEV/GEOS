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
 * @file FaceManagerT.h
 * @author settgast1
 * @date Feb 4, 2011
 */

#ifndef FACEMANAGERT_H_
#define FACEMANAGERT_H_

#include "math/TensorT/TensorT.h"
//#include "ObjectManagers/ExternalFaceStructs.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ArrayT/bufvector.h"
//#include "IO/silo/SiloFile.h"


#include "dataRepository/SynchronizedGroup.hpp"


struct ExternalFaceStruct;
class SiloFile;

#include <memory>
//#include "NestedRelation.h"



class ElementRegionT;
class NodeManagerT;
class ElementManagerT;
class PhysicalDomainT;
class EdgeManagerT;
class ExternalFaceManagerT;
class CohesiveZoneBase;

class FaceManagerT: public geosx::dataRepository::SynchronizedGroup
{
public:

  FaceManagerT();
  virtual ~FaceManagerT();

  static std::string CatalogName()
  {
    return "FaceManagerT";
  }

  void Initialize(  ){}

  globalIndex resize( const localIndex size, const bool assignGlobals=false )
  {
    const globalIndex firstNewGlobalIndex = ObjectDataStructureBaseT::resize(size,assignGlobals);
    m_toElementsRelation.resize(size);
    m_DataLengths = size;

    return firstNewGlobalIndex;
  }

  void BuildFaces( const NodeManagerT& nodeManager, const ElementManagerT& elemManager );

  void AddToFaceToElementMap( const ElementManagerT& elementManager,
                              const std::map<std::string,lArray1d>& newElementIndices );

//  void DetermineGlobalFaceNumbers( const NodeManagerT& nodeManager );

  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         Array1dT<gArray1d>& objectToCompositionObject );


  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL);
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL);

  realT FaceCenter( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& center ) const;
  R1Tensor FaceNormal( const NodeManagerT& nodeManager, const localIndex faceIndex, const bool referenceFlag = false ) const;
  realT FaceNormal( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& normal) const;

  void FaceTangential ( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& tanA, R1Tensor& tanB) const;

  realT FaceCenterAndNormal( const NodeManagerT& nodeManager,
                             const localIndex faceIndex,
                             R1Tensor& center,
                             R1Tensor& normal,
                             const bool referenceState = false) const;

  /// surface area - sum areas of triangular facets on the face.
  realT SurfaceArea( const NodeManagerT& nodeManager,
                     const localIndex faceIndex,
                     const bool referenceFlag = false ) const;
  
  /// Area of the face projected onto a surface with the given normal
  realT ProjectedArea( const NodeManagerT& nodeManager, const localIndex faceIndex, const R1Tensor& norm) const;
  
  R1Tensor CalculateGapVector( const NodeManagerT& nodeManager, const localIndex faceIndex ) const;

  R1Tensor CalculateGapRateVector( const NodeManagerT& nodeManager, const localIndex faceIndex ) const;

  void CalculateGapVectorDerivative( const NodeManagerT& nodeManager, const localIndex faceIndex,
                                     Array1dT< Array1dT<R1Tensor> >& gapDerivative,
                                     Array1dT<iArray1d>& gapDerivativeIndices ) const;

  void EdgeVectors( const NodeManagerT& nodeManager, const localIndex faceIndex, Array1dT<R1Tensor>& edgeVectors );

  /// Calculates the vector from node 0 to node 1 in 2D.  Note that 2D face is actually an edge.
  void FaceVector2D(const NodeManagerT& nodeManager, localIndex iFace, localIndex iNd, R1Tensor& v) ;

  void InFaceVectorNormalToEdge(const NodeManagerT& nodeManager,
                                const EdgeManagerT& edgeManager,
                                localIndex iFace,
                                localIndex iEdge,
                                R1Tensor& v) ;

  template< typename T_indices >
  unsigned int PackFaces( const T_indices& sendfaces,
                          const NodeManagerT& nodeManager,
                          const EdgeManagerT* const edgeManager,
                          bufvector& buffer,
                          const bool packConnectivityToGlobal,
                          const bool packFields,
                          const bool packMaps,
                          const bool packSets  ) const;

  unsigned int UnpackFaces( const char*& buffer,
                            const NodeManagerT& nodeManager,
                            const EdgeManagerT* const edgeManager,
                            ExternalFaceManagerT* const externalFaceManager,
                            lArray1d& faceReceiveLocalIndices,
                            const bool unpackConnectivityToLocal,
                            const bool unpackFields,
                            const bool unpackMaps,
                            const bool unpackSets  );

  void ConnectivityFromGlobalToLocal( const lSet& indices,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& edgeGlobalToLocal );


  void AddNewFace( const localIndex& k,
                   const localIndex& kelf,
                   localIndex& numFaces,
                   Array1dT<lArray1d>& facesByLowestNode,
                   lArray1d& tempNodeList,
                   Array1dT<lArray1d>& tempFaceToNodeMap,
                   const ElementRegionT& elementRegion );

  /// Calculate the bounding sphere for a face
  void FaceBoundingSphere( const NodeManagerT& nodeManager,
                           const localIndex faceIndex,
                           R1Tensor& center,
                           realT& radius) const;

  void FaceBoundingSphere(
                  const NodeManagerT& nodeManager,
                  const localIndex faceIndex,
                  const realT dt,
                  R1Tensor& center,
                  R1Tensor& velocity,
                  realT& radius,
                  R1Tensor& normal,
                  const bool referenceState = false) const;

  void NodalPositions(  const NodeManagerT& nodeManager,
                        const localIndex faceIndex,
                        Array1dT<R1Tensor>& xs) const;

  void FaceProperties( const NodeManagerT& nodeManager,
                       const localIndex faceIndex,
                       const realT dt,
                       R1Tensor & xfc,
                       R1Tensor & dxfc,
                       ExternalFaceStruct& efs) const;

  void FaceProperties( const NodeManagerT& nodeManager,
                       const localIndex faceIndex,
                       const realT dt,
                       ExternalFaceStruct& efs) const;

#ifdef STATES_ON_CONTACTS
  void FaceProperties( const NodeManagerT& nodeManager,
                       const localIndex faceIndex,
                       const realT dt,
                       R1Tensor & nx,
                       R1Tensor & xmin,
                       R1Tensor & xmax,
                       realT & area,
                       realT & length,
                       Array1dT<R1Tensor> & xs,
                       Array1dT<R1Tensor> & dxs) const;
#endif

  int ParentNormalDirection( const int lfn );
  

  void PreSeparateFaces( const std::string& setname, const int setState );

  void UpdateRuptureStates( const ElementManagerT& elementManager,
                            const NodeManagerT& nodeManager,
                            const std::string& separableSet,
                            const realT failval );

  void UpdateRuptureState( const ElementManagerT& elementManager,
                            const NodeManagerT& nodeManager,
                            const localIndex kf,
                            const realT failval );


  void CalculateStressOnFace( const ElementManagerT& elementManager,
                            const NodeManagerT& nodeManager,
                            const localIndex kf);

  void CalculateStressOnFace( const ElementManagerT& elementManager,
                            const NodeManagerT& nodeManager,
                            const localIndex kf,
                            realT& stressNOnFace,
                            R1Tensor& stressTOnFace);

  /// Check if face elements pass a given test, if one does it returns true and the first element to pass, otherwise it returns false.
  template<typename UnaryElementBooleanFunction>
  bool GetFirstMatchingElement(localIndex face, UnaryElementBooleanFunction criteria, ElementIdPair& match );
  
  /// Compares elements on either side of the face, returns index of element that best matches the criteria
  template<typename BinaryElementComparisonFunction>
  int  GetMatchingElementIndex(localIndex face, BinaryElementComparisonFunction criteria);
  
  /// Rearrange face node order counter-clockwise around face
  void SortFaceNodes(const NodeManagerT& nodeManager, const localIndex faceIndex );
  
  /// Rearrange all face node orders counter-clockwise around face
  void SortAllFaceNodes(const NodeManagerT& nodeManager);
  
  /// Copy face field values to their respective nodes using a weighted average. 
  template<typename T, typename WeightFuncPtr>
  void CopyFieldToNodes(const std::string& faceField, NodeManagerT& nodeManager,const std::string& nodeField, WeightFuncPtr weightFunc) const;
  /// Copy values from a subset of faces to their nodes using a weighted average. 
  template<typename T, typename WeightFuncPtr>
  void CopyFieldToNodes(const std::string& faceField, const lArray1d& set, NodeManagerT& nodeManager, const std::string& nodeField,
                        WeightFuncPtr weightFunctionPtr) const;
                                    
  /// Copy face field to nodes - overloaded functions when field name is the same for node and face fields.
  template<typename T, typename WeightFuncPtr>
  void CopyFieldToNodes(const std::string& fieldName, NodeManagerT& nodeManager, WeightFuncPtr weightFunc) const
     { CopyFieldToNodes<T>(fieldName, nodeManager,fieldName, weightFunc); }
  template<typename T, typename WeightFuncPtr>
  void CopyFieldToNodes(const std::string& fieldName, const lArray1d& set, NodeManagerT& nodeManager, WeightFuncPtr weightFunc) const
     { CopyFieldToNodes<T>(fieldName, set, nodeManager,fieldName, weightFunc); }

  void SplitFace( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  const lSet& splitEdges,
                  const Array1dT<lArray1d>& childEdges,
                  Array1dT<lSet>& nodesToFaces,
                  Array1dT<lSet>& edgesToFaces );

  void ModifyToFaceMapsFromSplit( const lSet& newFaces,
                                  const lSet& modifiedFaces,
                                  NodeManagerT& nodeManager,
                                  EdgeManagerT& edgeManager,
                                  ExternalFaceManagerT& externalFaceManager );


  void SetApertureFromRigidWall( const NodeManagerT& nodeManager );


  void WriteSiloMesh( SiloFile& siloFile,
                      const std::string& meshname,
                      const NodeManagerT& nodeManager,
                      const int cycleNum,
                      const realT problemTime,
                      const bool isRestart );

  bool IsNodeOnFace( const localIndex faceIndex,
                     const localIndex nodeIndex );


protected:
  void WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                         const std::string& siloDirName,
                                         const std::string& meshname,
                                         const int centering,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const std::string& multiRoot,
                                         const std::string& regionName = "none",
                                         const lArray1d& mask = lArray1d());

  void ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                          const std::string& siloDirName,
                                          const std::string& meshname,
                                          const int centering,
                                          const int cycleNum,
                                          const realT problemTime,
                                          const bool isRestart,
                                          const std::string& regionName = "none",
                                          const lArray1d& mask = lArray1d());

public:


  const localIndex& m_numFaces;
  OrderedVariableOneToManyRelation& m_toNodesRelation;
  OrderedVariableOneToManyRelation& m_toEdgesRelation;
  Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > > m_toElementsRelation;
  

  lSet& m_externalFaces;

  std::set< localIndex > m_matchedBoundaryFaces;

#if USECPP11==1
  std::unique_ptr<CohesiveZoneBase> m_cohesiveZone;
#else
  CohesiveZoneBase* m_cohesiveZone;
#endif


  ElementManagerT* m_elemManagerHACK;

  int m_writeArbitraryPolygon;

private:

};

/// returns true if an element matching the criteria in the ElementUnaryComparisonFunc is found
/// the id pair for the matching element is copied into "match"
///
/// The UnaryElementBooleanFunction is a function of the form:
/// bool f(const ElementIdPair& )
template<typename UnaryElementBooleanFunction>
bool FaceManagerT::GetFirstMatchingElement(localIndex face, UnaryElementBooleanFunction criteria, ElementIdPair& match ){
  bool rv = false; 
  if( m_toElementsRelation[face].size() == 1 ) {
     match = m_toElementsRelation[face][0];
     rv = criteria(match);
  } else if( m_toElementsRelation[face].size() == 2 ){
     match = m_toElementsRelation[face][0];
     rv = criteria(match);
     if( !rv ){
       match = m_toElementsRelation[face][1];
       rv = criteria(match);
     } 
  } 
        
  return rv;
}

/// Compares elements on either side of face using the given criteria 
///
/// The BinaryElementComparisonFunction is a function of the form:
/// bool f(constElementIdPair , const ElementIdPair&)
///
/// Returns 0 if criteria( element0, element1) is true or there is only one elment on the face
/// Returns -1 to indicate unintialized face
/// Returns 1 if there are two elements on the face and criteria( element0, element1) is false 
///
/// Use this function if you need to compare two elements and pick the better one
/// Use GetFirstMatchingElement if you need to check if none of the elements meet a given criteria

template<typename BinaryElementComparisonFunction>
int FaceManagerT::GetMatchingElementIndex(localIndex face, BinaryElementComparisonFunction criteria){
  int i = -1; 
  if( m_toElementsRelation[face].size() == 1 ) {
     i = 0;
  } else if( m_toElementsRelation[face].size() == 2 ){
     ElementIdPair& epa = m_toElementsRelation[face][0];
     ElementIdPair& epb = m_toElementsRelation[face][1];
     
     bool isFirst = criteria(epa,epb);
     if( isFirst ){
       i = 0;
     } else {
       i = 1;	
     }
  } 
        
  return i;
}

/*
 * 
 * 
 * @param weightFunctionPtr Pointer to a function of type 
 *         realT (const FaceManagerT&, localIndex faceId, NodeManagerT& , localIndex nodeId) 
 *  The function determines how to weigh each faces's contribution to the 
 *  nodal values.
 * 
 * 
 * Node values are set according to:
 *   nodeValue = \sum_{i=faces} w_{i} faceField_{i} / \sum_{j=faces} w_{j}  

template<typename T, typename WeightFuncPtr>
void FaceManagerT::CopyFieldToNodes(const std::string& faceFieldStr, 
                                    NodeManagerT& nodeManager, const std::string& nodeFieldStr,
                                    WeightFuncPtr weightFunctionPtr) const{

  const Array1dT<T>& faceField = this->GetFieldData<T>(faceFieldStr);
  Array1dT<T>& nodeField = nodeManager.GetFieldData<T>(nodeFieldStr);
  nodeField = 0;
  std::vector<realT> weights(nodeManager.m_numNodes,0.0);
  
  for(int f =0; f < m_numFaces; ++f){
  	const lArray1d& faceNodeMap = m_FaceToNodeMap[f];
    for( localIndex a=0 ; a<faceNodeMap.size() ; ++a )
    {
      const localIndex nd = faceNodeMap[a];
      realT w = (*weightFunctionPtr)(&this,f, nodeManager, nd);
      nodeField[nd] += w*faceField[f];
      weights[nd] += w;
    }
  }
  
  for(int nd =0; nd < nodeManager.m_numNodes; ++nd){
  	if(weights[nd] != 0) nodeField[nd] /= weights[nd];
  }
}


template<typename T, typename WeightFuncPtr>
void FaceManagerT::CopyFieldToNodes(const std::string& faceFieldStr, const lArray1d& set,
                                    NodeManagerT& nodeManager, const std::string& nodeFieldStr,
                                    WeightFuncPtr weightFunctionPtr) const{

  const Array1dT<T>& faceField = this->GetFieldData<T>(faceFieldStr);
  Array1dT<T>& nodeField = nodeManager.GetFieldData<T>(nodeFieldStr);
    
  localIndex setSize = set.size();
  
  std::map<localIndex,realT > weights;
  std::pair<std::map<localIndex,realT>::iterator,bool> ret;
  for(localIndex i =0; i < setSize; ++i){
    const localIndex f = set[i]; 
  	const lArray1d& faceNodeMap = m_FaceToNodeMap[f];
    for( localIndex a=0 ; a<faceNodeMap.size(); ++a )
    {
      const localIndex nd = faceNodeMap[a];
      realT w = (*weightFunctionPtr)(&this,f, nodeManager, nd);
                  
      ret=weights.insert (std::pair<localIndex,realT> (nd,w) );
      if (ret.second) {
         // node did not exist in map (weight was inserted)
         nodeField[nd] = w*faceField[f];
      } else {
      	 // node did exist in map (new weight not yet added)
      	 ret.first->second += w;  // weights[nd] += w;
         nodeField[nd] += w*faceField[f];
      }
    } 
  }
  
  std::map<localIndex,realT>::iterator iend = weights.end();
  for(std::map<localIndex,realT>::iterator itr = weights.begin(); itr != iend; ++itr){
  	if(itr->second != 0.0) nodeField[itr->first] /= itr->second;
  }
}
 */
#endif /* FACEMANAGERT_H_ */
