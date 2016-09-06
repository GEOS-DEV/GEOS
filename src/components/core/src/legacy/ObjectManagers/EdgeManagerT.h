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
 * @file EdgeManagerT.h
 * @author settgast1
 * @date Jun 22, 2011
 */

#ifndef EDGEMANAGERT_H_
#define EDGEMANAGERT_H_

//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "dataRepository/SynchronizedGroup.hpp"


class FaceManagerT;
class NodeManagerT;

class EdgeManagerT: public ObjectDataStructureBaseT
{
public:
  EdgeManagerT();
  ~EdgeManagerT();

  void Initialize() {}

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL) ;
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL) ;
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         Array1dT<gArray1d>& objectToCompositionObject );


  void BuildEdges( const FaceManagerT& faceManager, const NodeManagerT& nodeManager );
  void BuildEdges( const ElementManagerT& elementManager, const NodeManagerT& nodeManager );
  
  template< typename T_indices >
  unsigned int PackEdges( const T_indices& sendedges,
                          const NodeManagerT& nodeManager,
                          const FaceManagerT& faceManager,
                          bufvector& buffer,
                          const bool packConnectivityToGlobal,
                          const bool packFields,
                          const bool packMaps,
                          const bool packSets  ) const;

  unsigned int UnpackEdges( const char*& buffer,
                            const NodeManagerT& nodeManager,
                            const FaceManagerT& faceManager,
                            lArray1d& edgeReceiveLocalIndices,
                            const bool unpackConnectivityToLocal,
                            const bool unpackFields,
                            const bool unpackMaps,
                            const bool unpackSets  );

  void ConnectivityFromGlobalToLocal( const lSet& indices,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );

//  void UpdateEdgeExternalityFromSplit( const FaceManagerT& faceManager,
//                                     const lSet& newEdgeIndices,
//                                     const lSet& modifiedEdgeIndices );

  void EdgeCenter(const NodeManagerT& nodeManager, localIndex edge, R1Tensor& center)const;
  void EdgeVector(const NodeManagerT& nodeManager, localIndex edge, R1Tensor& vector)const;
  realT EdgeLength(const NodeManagerT& nodeManager, localIndex edge) const;

  void AddToEdgeToFaceMap( const FaceManagerT& faceManager,
                           const lArray1d& newFaceIndices );

  void SplitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  Array1dT<lSet>& nodesToEdges );

  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

  localIndex FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex nodeB, const NodeManagerT& nodeManager);

  void SetLayersFromDomainBoundary(const NodeManagerT& nodeManager);


  FixedOneToManyRelation& m_toNodesRelation;
  UnorderedVariableOneToManyRelation& m_toFacesRelation;

};

#endif /* EDGEMANAGERT_H_ */
