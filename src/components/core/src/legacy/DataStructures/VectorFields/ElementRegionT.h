/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTOBJECTT_H_
#define ELEMENTOBJECTT_H_

#include "managers/ObjectManagerBase.hpp"
#include "Common/Common.h"
//#include "ObjectDataStructureBaseT.h"
#include "StableTimeStep.h"
//#include "IO/ticpp/HierarchicalDataNode.h"
#include "ObjectManagers/EnergyT.h"

class IntegrationRuleT;
class MaterialBaseParameterDataT;
class MaterialBaseStateDataT;
class MaterialBaseT;
class NodeManager;
class FaceManagerT;
class EdgeManagerT;
class ElementRegionT;
class PhysicalDomainT;

class FiniteElementBase;
class QuadratureBase;
class BasisBase;

class MaterialBase;

namespace geosx
{
/**
 * Class to manage the data stored at the element level.
 */
class ElementRegionT : public ObjectDataStructureBaseT
{
public:



  static const char ElementObjectToElementManager[];
  static const char ElementToNode[];
  static const char ElementToFace[];
  static const char ElementToEdge[];

//  static const char ElementToCrackSurface[];
//  static const char ElementToCrackSurfaceVertex[];
//  static const char ElementToLocalVolume[];
//  static const char ElementToPhysicalNodes[];
//  static const char ElementToCrackToVertexNodes[];
//  static const char ElementToCracks[];


  ElementRegionT( ObjectManangerBase * const parent );


  ElementRegionT(const ElementRegionT& init);


  virtual ~ElementRegionT();

  void Initialize( const PhysicalDomainT& domain );

  virtual void DeserializeObjectField(const std::string& ifield, const array<real64>& field);
  virtual void DeserializeObjectFields(const array<string>& names, const array<array<real64> >& fields);

  void ReadXML( TICPP::HierarchicalDataNode* const hdn,
                const bool isRestart );

  globalIndex resize( const localIndex size, const bool assignGlobals=false );

  void SetGeometryBasedVariables();

  void AllocateElementLibrary( const int basis,
                               const int quadrature );

  void Initialize();

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL);
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL) { (void)referenceObject; }
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<gArray1d>& objectToCompositionObject  )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("ElementRegionT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }


  void ModifyToElementMapsFromSplit( const lSet& modifiedElements,
                                     NodeManager& nodeManager,
                                     FaceManagerT& faceManager );

  void UpdateExternalityFromSplit( const lSet& modifiedElements,
                                   NodeManager& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager );



  int CalculateShapeFunctionDerivatives( const NodeManager& nodeManager );

//  int CalculateShapeFunctionDerivativesCutElements(const NodeManagerT&
// nodeManager);


  int CalculateVelocityGradients( const NodeManager& nodeManager, const int calcGroup = 0 );

  int MaterialUpdate(const realT dt);

  int CalculateSmallDeformationNodalForces( NodeManager& nodeManager,
                                            StableTimeStep& timeStep,
                                            const realT dt );
  int CalculateNodalForces( NodeManager& nodeManager,
                            StableTimeStep& timeStep,
                            const realT dt );

  int CalculateNodalForcesFromOneElement(const localIndex nodeID,
                                         const localIndex elemID,
                                         NodeManager& nodeManager,
                                         R1Tensor& fNode );

  realT ElementGDivBeta(const localIndex elemID);

  void CalculateNodalForceFromStress(const localIndex elemID,
                                     const NodeManager& nodeManager,
                                     R2SymTensor& stress,
                                     array<R1Tensor>& fNode);

  int ProcessElements(  NodeManager& nodeManager,
                        StableTimeStep& timeStep,
                        const realT dt  );


  int CalculateNodalMasses( NodeManager& nodeManager );


  void SetIsAttachedToSendingGhostNode( const NodeManager& nodeManager );

  template< typename T_indices >
  unsigned int PackElements( bufvector& buffer,
                             lSet& sendnodes,
                             lSet& sendfaces,
                             const T_indices& elementList,
                             const NodeManager& nodeManager,
                             const FaceManagerT& faceManager,
                             const bool packConnectivityToGlobal,
                             const bool packFields,
                             const bool packMaps,
                             const bool packSets  ) const;

  unsigned int UnpackElements( const char*& buffer,
                               const NodeManager& nodeManager,
                               const FaceManagerT& faceManager,
                               lArray1d& elementRegionReceiveLocalIndices,
                               const bool unpackConnectivityToLocal,
                               const bool unpackFields,
                               const bool unpackMaps,
                               const bool unpackSets  );

  void ConnectivityFromGlobalToLocal( const lSet& list,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );


  template< typename T_indices >
  unsigned int PackFieldsIntoBuffer( bufvector& buffer,
                                     const array<string>& fieldNames,
                                     const T_indices& localIndices,
                                     const bool doBufferPacking=true ) const;

  template< typename T_indices >
  unsigned int PackFieldsIntoBuffer( char*& buffer,
                                     const array<string>& fieldNames,
                                     const T_indices& localIndices,
                                     const bool doBufferPacking=true ) const;

  unsigned int UnpackFieldsFromBuffer( const char*& buffer,
                                       const array<string>& fieldNames,
                                       const lArray1d& localIndices );

  template< typename T_indices >
  unsigned int PackAllFieldsIntoBuffer( bufvector& buffer,
                                        const T_indices& localIndices ) const;

  unsigned int UnpackAllFieldsFromBuffer( const char*& buffer,
                                          const lArray1d& localIndices );

  void UpdateElementFieldsWithGaussPointData();

  lArray2d& ElementToNodeMap();
  localIndex* ElementToNodeMap( const int elemNum );
  const localIndex* ElementToNodeMap( const int elemNum ) const;
  localIndex ElementToNodeMap( const int elemNum, const int nodeNum ) const;

  void GetFaceNodes( const localIndex elementIndex,
                     const localIndex localFaceIndex,
                     lArray1d& nodeIndicies ) const;

  void GetElementNeighbors(localIndex el, 
                           const FaceManagerT& faceManager, 
                           set<localIndex>& neighbors) const;
                     
  R1Tensor GetElementCenter(localIndex k, const NodeManager& nodeManager, const bool useReferencePos = false) const;


  void WriteSiloRegionMesh( SiloFile& siloFile,
                            const std::string& meshname,
                            const int cycleNum,
                            const realT problemTime,
                            const bool isRestart,
                            const std::string& regionName );

  void ReadSiloRegionMesh( const SiloFile& siloFile,
                           const std::string& meshname,
                           const int cycleNum,
                           const realT problemTime,
                           const bool isRestart,
                           const std::string& regionName );

  /// returns true if element id pair is inside region
  bool ContainsElement(const ElementIdPair& ep);

  array<integer> SiloNodeOrdering();

  bool SplitObject( const localIndex indexToSplit,
                    const int rank,
                    localIndex newIndices[2],
                    const bool forceSplit );


  std::string m_regionName;
  int m_regionNumber;
  const localIndex& m_numElems;
  unsigned int m_numNodesPerElem;
  unsigned int m_numIntegrationPointsPerElem;

  std::string m_elementType;
  std::string m_elementGeometryID;
  std::vector<std::string> m_parentFaceSetNames;

  unsigned int m_ElementDimension;


//  lArray1d& m_ElementObjectToElementManagerMap;
  FixedOneToManyRelation& m_toNodesRelation;
  FixedOneToManyRelation& m_toFacesRelation;

  //FixedOneToManyRelation& m_toEdgesRelation;  //This is neither filled nor
  // used.
//  OrderedVariableOneToManyRelation& m_toPhysicalNodes;
//  OrderedVariableOneToManyToManyRelation& m_toCrackToVertexNodes;
//  OrderedVariableOneToManyRelation& m_toCracks;

//  //////////////////////////////////////////////
//  // The following need to be deleted
//  UnorderedVariableOneToManyRelation& m_toCrackSurfacesRelation;
//  UnorderedVariableOneToManyRelation& m_toCrackSurfaceVerticesRelation;
//  UnorderedVariableOneToManyRelation& m_toLocalVolumeRelation;
//  /////////////////////////////////////////////

  array< Array2dT<R1Tensor> > m_dNdX;
  rArray2d m_detJ;
  rArray2d m_detJ_n;
  rArray2d m_detJ_np1;
  Array2dT< R2Tensor >  m_dUdX;
  Array2dT< R2Tensor > m_Finv;
  array<real64> m_Kregion;

  array< array<R2SymTensor> > m_Dadt;
  array< array<R2Tensor> >  m_Rot;


  array< rArray2d > m_Ke;
  array< rArray2d > m_matrixB;
  array< rArray2d > m_matrixE;



  unsigned int m_basis;
  unsigned int m_quadrature;

  FiniteElementBase* m_finiteElement;
  QuadratureBase* m_elementQuadrature;
  BasisBase* m_elementBasis;

//  MaterialManagerT m_material;

//  MaterialBaseT* m_materialComputations;


  //const int m_numFacesPerElement;
  //const int m_numNodesPerFace;
  int m_numFacesPerElement;
  int m_numNodesPerFace;

  EnergyT m_energy;

  realT m_hgDamp;
  realT m_hgStiff;
  realT m_failStress;

  bool m_plotMat;

#if USECPP11==1
  std::unique_ptr<MaterialBase> m_mat;
#else
  MaterialBase* m_mat;
#endif


private:
  ElementRegionT& operator=(const ElementRegionT& rhs);



};


inline lArray2d& ElementRegionT::ElementToNodeMap()
{
  return m_toNodesRelation;
}


inline localIndex* ElementRegionT::ElementToNodeMap( const int elemNum )
{
  return &(m_toNodesRelation(elemNum,0));
}

inline const localIndex* ElementRegionT::ElementToNodeMap( const int elemNum ) const
{
  return &(m_toNodesRelation(elemNum,0));
}

inline localIndex ElementRegionT::ElementToNodeMap( const int elemNum, const int nodeNum ) const
{
  return m_toNodesRelation(elemNum,nodeNum);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////



}



#endif /* ELEMENTOBJECTT_H_ */
