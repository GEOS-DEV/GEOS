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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SurfaceGenerator.hpp
 */
#ifndef SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_
#define SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{

struct ModifiedObjectLists
{
  set<localIndex> newNodes;
  set<localIndex> newEdges;
  set<localIndex> newFaces;
  set<localIndex> modifiedNodes;
  set<localIndex> modifiedEdges;
  set<localIndex> modifiedFaces;
  std::map< std::string, set<localIndex> > modifiedElements;
};


class SpatialPartition;

class NodeManager;
class EdgeManager;
class FaceManager;
class ExternalFaceManager;
class ElementRegionManager;
class ElementRegion;

/**
 * @class SurfaceGenerator
 *
 * This solver manages the mesh topology splitting methods.
 *
 */
class SurfaceGenerator : public SolverBase
{
public:
  SurfaceGenerator( const std::string& name,
                    ManagedGroup * const parent );
  ~SurfaceGenerator() override;


  static string CatalogName() { return "SurfaceGenerator"; }

  virtual void FillDocumentationNode() override;

  virtual void
  FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup ) override;

  virtual void FinalInitializationPreSubGroups( ManagedGroup * const problemManager ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual real64 SolverStep( real64 const& time_n,
                             real64 const& dt,
                             integer const cycleNumber,
                             DomainPartition * domain ) override;
//
//  virtual void ImplicitStepSetup( real64 const& time_n,
//                              real64 const& dt,
//                              DomainPartition * const domain,
//                              systemSolverInterface::EpetraBlockSystem * const blockSystem ) override;
//
//
//  virtual void AssembleSystem( DomainPartition * const domain,
//                               systemSolverInterface::EpetraBlockSystem * const blockSystem,
//                               real64 const time,
//                               real64 const dt ) override;
//
//  virtual void ApplyBoundaryConditions( DomainPartition * const domain,
//                                        systemSolverInterface::EpetraBlockSystem * const blockSystem,
//                                        real64 const time,
//                                        real64 const dt ) override;
//
//  virtual real64
//  CalculateResidualNorm( systemSolverInterface::EpetraBlockSystem const * const blockSystem ) override;

//  virtual void SolveSystem( systemSolverInterface::EpetraBlockSystem * const blockSystem,
//                            SystemSolverParameters const * const params ) override;

//  virtual void
//  ApplySystemSolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
//                       real64 const scalingFactor,
//                       DomainPartition * const domain ) override;

//  virtual void ResetStateToBeginningOfStep( DomainPartition * const domain ) override;
//
//  virtual  void ImplicitStepComplete( real64 const & time,
//                                      real64 const & dt,
//                                      DomainPartition * const domain ) override;
/**@}*/


  int SeparationDriver( NodeManager& nodeManager,
                        EdgeManager& edgeManager,
                        FaceManager& faceManager,
                        ElementRegionManager& elementManager,
                        SpatialPartition& partition,
                        const bool prefrac,
                        const realT time );
private:


private:


  /**
   * @brief Function to identify which faces are ready for rupture
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param partition
   * @param prefrac
   */
  void IdentifyRupturedFaces( NodeManager & nodeManager,
                              EdgeManager & edgeManager,
                              FaceManager & faceManager,
                              ElementRegionManager & elementManager,
                              SpatialPartition& partition,
                              const bool prefrac );

  /**
   * @brief
   * @param edgeID
   * @param trailFaceID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param vecTipNorm
   * @param vecTip
   * @return
   */
  realT CalculateEdgeSIF ( const localIndex edgeID,
                           localIndex& trailFaceID,
                           NodeManager & nodeManager,
                           EdgeManager & edgeManager,
                           FaceManager & faceManager,
                           ElementRegionManager & elementManager,
                           R1Tensor& vecTipNorm,
                           R1Tensor& vecTip );

  /**
   * @brief
   * @param edgeID
   * @param trailFaceID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param vecTipNorm
   * @param vecTip
   * @param modifiedObjects
   * @param edgeMode
   */
  void MarkRuptureFaceFromEdge ( const localIndex edgeID,
                                 localIndex& trailFaceID,
                                 NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager,
                                 R1Tensor& vecTipNorm,
                                 R1Tensor& vecTip,
                                 ModifiedObjectLists& modifiedObjects,
                                 const int edgeMode );

  /**
   *
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   */
  void PostUpdateRuptureStates( NodeManager & nodeManager,
                                EdgeManager & edgeManager,
                                FaceManager & faceManager,
                                ElementRegionManager & elementManager,
                                array1d<set<localIndex> >& nodesToRupturedFaces,
                                array1d<set<localIndex> >& edgesToRupturedFaces );

  /**
   *
   * @param edgeID
   * @param nodeManager
   * @param faceManager
   * @param edgeManager
   * @param prefrac
   * @return
   */
  int CheckEdgeSplitability( const localIndex edgeID,
                             NodeManager & nodeManager,
                             FaceManager & faceManager,
                             EdgeManager & edgeManager,
                             const bool prefrac );

  /**
   *
   * @param nodeID
   * @param nodeManager
   * @param faceManager
   * @param edgeManager
   * @param prefrac
   * @return
   */
  int CheckNodeSplitability( const localIndex nodeID,
                             NodeManager & nodeManager,
                             FaceManager & faceManager,
                             EdgeManager & edgeManager,
                             const bool prefrac );

//  void UpdatePathCheckingArrays();

  /**
   * @brief check and split node in mesh
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elemManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param elementManager
   * @param modifiedObjects
   * @param prefrac
   * @return
   */
  bool ProcessNode( const localIndex nodeID,
                    NodeManager & nodeManager,
                    EdgeManager & edgeManager,
                    FaceManager & faceManager,
                    ElementRegionManager & elemManager,
                    arrayView1d<set<localIndex> >& nodesToRupturedFaces,
                    arrayView1d<set<localIndex> >& edgesToRupturedFaces,
                    ElementRegionManager & elementManager,
                    ModifiedObjectLists& modifiedObjects,
                    const bool prefrac );

  /**
   * @brief Find a fracture path for surface generation
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elemManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param separationPathFaces
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool FindFracturePlanes( const localIndex nodeID,
                           const NodeManager & nodeManager,
                           const EdgeManager & edgeManager,
                           const FaceManager & faceManager,
                           ElementRegionManager & elemManager,
                           const arrayView1d<set<localIndex> >& nodesToRupturedFaces,
                           const arrayView1d<set<localIndex> >& edgesToRupturedFaces,
                           set<localIndex>& separationPathFaces,
                           map<localIndex, int>& edgeLocations,
                           map<localIndex, int>& faceLocations,
                           map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );


  /**
   * @brief given a fracture path, split the mesh according to the path.
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param modifiedObjects
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param separationPathFaces
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   */
  void PerformFracture( const localIndex nodeID,
                        NodeManager & nodeManager,
                        EdgeManager & edgeManager,
                        FaceManager & faceManager,
                        ElementRegionManager & elementManager,
                        ModifiedObjectLists& modifiedObjects,
                        arrayView1d<set<localIndex> >& nodesToRupturedFaces,
                        arrayView1d<set<localIndex> >& edgesToRupturedFaces,
                        const set<localIndex>& separationPathFaces,
                        const map<localIndex, int>& edgeLocations,
                        const map<localIndex, int>& faceLocations,
                        const map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );

  /**
   * @brief function to set which side of the fracture plane all objects are on
   * @param separationPathFaces
   * @param elemManager
   * @param faceManager
   * @param nodesToElements
   * @param localFacesToEdges
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool SetLocations( const set<localIndex>& separationPathFaces,
                     ElementRegionManager & elemManager,
                     const FaceManager & faceManager,
                     const set< std::pair<CellBlockSubRegion*, localIndex> >& nodesToElements,
                     const map< localIndex, std::pair<localIndex, localIndex> >& localFacesToEdges,
                     map<localIndex, int>& edgeLocations,
                     map<localIndex, int>& faceLocations,
                     map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );

  /**
   * @brief function to set which side of the fracture plane all objects are on
   * @param side
   * @param elem
   * @param separationPathFaces
   * @param elemManager
   * @param faceManager
   * @param nodesToElements
   * @param localFacesToEdges
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool SetElemLocations( const int side,
                         const std::pair<CellBlockSubRegion*, localIndex >& elem,
                         const set<localIndex>& separationPathFaces,
                         ElementRegionManager & elemManager,
                         const FaceManager & faceManager,
                         const set< std::pair<CellBlockSubRegion*, localIndex> >& nodesToElements,
                         const map< localIndex, std::pair<localIndex, localIndex> >& localFacesToEdges,
                         map<localIndex, int>& edgeLocations,
                         map<localIndex, int>& faceLocations,
                         map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );

  /**
   *
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  realT CalculateKinkAngle ( const localIndex nodeID,
                             const NodeManager & nodeManager,
                             EdgeManager & edgeManager,
                             FaceManager & faceManager );

  /**
   *
   * @param faceManager
   * @param edgeManager
   * @param nodeManager
   * @param modifiedObjects
   * @param prefrac
   */
  void CalculateKinkAngles ( FaceManager & faceManager,
                             EdgeManager & edgeManager,
                             NodeManager & nodeManager,
                             ModifiedObjectLists& modifiedObjects,
                             const bool prefrac );

  /**
   *
   * @param edgeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  realT MinimumToughnessOnEdge( const localIndex edgeID,
                                const NodeManager & nodeManager,
                                EdgeManager & edgeManager,
                                FaceManager & faceManager );

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto ruptureStateString = "ruptureState";
    constexpr static auto failCriterionString = "failCriterion";
    constexpr static auto degreeFromCrackString = "degreeFromCrack";
  }; //SurfaceGenViewKeys;

private:
  /// choice of failure criterion
  integer m_failCriterion=1;

  /// set of sepearable faces
  localIndex_set m_separableFaceSet;

  /// copy of the original node->face mapping prior to any separation
  array1d< set<localIndex> > m_originalNodetoFaces;

  /// copy of the original node->edge mapping prior to any separation
  array1d< set<localIndex> > m_originalNodetoEdges;

  /// copy of the original face->edge mapping prior to any separation
  array1d< array1d<localIndex> > m_originalFaceToEdges;

  /// collection of faces that have been used for separation of each node
  array1d< set<localIndex> > m_usedFacesForNode;

  /// copy of the original face->elemRegion mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemRegion;

  /// copy of the original face->elemSubRegion mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemSubRegion;

  /// copy of the original face->elemIndex mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemIndex;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_ */
