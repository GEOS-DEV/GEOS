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
 * @file NodeManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */


#ifndef NODEMANAGERT_H_
#define NODEMANAGERT_H_

#include "ObjectManagerBase.hpp"
#include <string.h>


// *********************************************************************************************************************
// *********************************************************************************************************************
class SiloFile;

namespace geosx
{

class ElementRegionT;
class FaceManager;
class EdgeManagerT;

/**
 * @author Randolph Settgast
 *
 * The NodeManagerT class manages the node data using the ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry corresponds to a node.
 */
class NodeManager : public ObjectDataStructureBaseT
{
public:


  /// default constructor
  NodeManager();



  /// default destructor
  ~NodeManager();


  globalIndex resize( const localIndex size, const bool assignGlobals = false )
  {
    m_toElementsRelation.resize(size);
    return ObjectDataStructureBaseT::resize(size, assignGlobals );
  }

  void Initialize();

  /// pure virtual function that sets what objects are on the boundary of the domain
  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL );
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL) ;
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager ,
                                                         Array1dT<gArray1d>& objectToCompositionObject  )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("NodeManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }

  void SetLayersFromDomainBoundary( const int layer );

  /// pack nodes into a buffer
  template< typename T_indices >
  unsigned int PackNodes( const T_indices& sendnodes,
                          const FaceManager& faceManager,
                          bufvector& buffer,
                          const bool packConnectivityToGlobal,
                          const bool packFields,
                          const bool packMaps,
                          const bool packSets  ) const;

  /// unpack nodes from a buffer
  unsigned int UnpackNodes( const char*& buffer,
                            const FaceManager& faceManager,
                            lArray1d& nodeReceiveLocalIndices,
                            const bool unpackConnectivityToLocal,
                            const bool unpackFields,
                            const bool unpackMaps,
                            const bool unpackSets  );

  void CalculateEffectiveNormal( const localIndex index,
                                 const FaceManager& faceManager,
                                 R1Tensor& normal ) const;


  /// copy fields from one node to another
  void CopyNode( const int destination, const int source );

  /// construct the nodeToElementMap using data in elementManager
  void ConstructNodeToElementMap( const ElementManagerT& elementManager );

  /// construct the nodeToElementMap using data in elementManager
  void AddToNodeToElementMap( const ElementManagerT& elementManager,
                                  const std::map<std::string,lArray1d>& newElementIndices );

  /// construct nodeToFaceMap using data in the faceManager
  void ConstructNodeToFaceMap( const FaceManager& faceManager );

  void ConnectivityFromGlobalToLocal( const lSet& indices,
                                      const lSet& clearIndices,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );


  void ModifyNodeToEdgeMapFromSplit( const EdgeManagerT& edgeManager,
                                     const lSet& newEdgeIndices,
                                     const lSet& modifiedEdgeIndices );

// Fu note on 20130416: Looks like this was temporary.  This function is now taken care of by element region. 
// We will keep it here for a while and delete it later.
//  void UpdateNodeExternalityFromSplit( const FaceManager& faceManager,
//                                     const lSet& newNodeIndices,
//                                     const lSet& modifiedNodeIndices );

  /// construct nodeToFaceMap using data in the faceManager
  void AddToNodeToFaceMap( const FaceManager& faceManager,
                           const lArray1d& newFaceIndices );

  void AddToNodeToEdgeMap( const EdgeManagerT& edgeManager,
                           const lArray1d& newEdgeIndices );

  /// get the current position of the node
  void GetPosition(localIndex nd, R1Tensor& position){
    Array1dT< R1Tensor >& refPosition = this->GetFieldData<FieldInfo::referencePosition>();
    Array1dT< R1Tensor >& displacement = this->GetFieldData<FieldInfo::displacement>();

    position = refPosition[nd];
    position += displacement[nd];
  }


  void SortNodeOnPlane (lArray1d& nostList) const;

  void UpdateDetachedNodeLocationAndVelocity();
  void ZeroDetachedNodeVelocity();
  void FindAllEffectiveChildren(localIndex& nodeID,
                                lArray1d& list,
                                rArray1d& weight);

  void GetDomainExtents(R1Tensor& pmin, R1Tensor& pmax, localIndex Ndims);

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



  /** @name Maps
   * The Maps
   */
  ///@{

  typedef std::set< std::pair<ElementRegionT*,localIndex> > nodeToElemType;
  Array1dT< nodeToElemType >  m_toElementsRelation;

  UnorderedVariableOneToManyRelation&  m_nodeToFaceMap;
  UnorderedVariableOneToManyRelation&  m_nodeToEdgeMap;

//  UnorderedVariableOneToManyRelation& m_toCrackSurfacesRelation;

  ///@}

  std::set< localIndex > m_matchedBoundaryNodes;


protected:

private:
  /// copy constructor
  NodeManager( const NodeManager& init );
  NodeManager& operator=( const NodeManager&);


};
// *********************************************************************************************************************
// *********************************************************************************************************************


// *********************************************************************************************************************


}


#endif /* NODEMANAGERT_H_ */
