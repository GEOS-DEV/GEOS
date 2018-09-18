/*
 * SurfaceGenerator.hpp
 *
 *  Created on: Jul 3, 2018
 *      Author: settgast
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

  virtual void FinalInitialization( ManagedGroup * const problemManager ) override final;


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


  void IdentifyRupturedFaces( NodeManager & nodeManager,
                              EdgeManager & edgeManager,
                              FaceManager & faceManager,
                              ElementRegionManager & elementManager,
                              SpatialPartition& partition,
                              const bool prefrac  );

  realT CalculateEdgeSIF ( const localIndex edgeID,
                           localIndex& trailFaceID,
                         NodeManager & nodeManager,
                         EdgeManager & edgeManager,
                         FaceManager & faceManager,
                         ElementRegionManager & elementManager,
                         R1Tensor& vecTipNorm,
                         R1Tensor& vecTip);

  void MarkRuptureFaceFromEdge ( const localIndex edgeID,
                                 localIndex& trailFaceID,
                                 NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager,
                                 R1Tensor& vecTipNorm,
                                 R1Tensor& vecTip,
                                 ModifiedObjectLists& modifiedObjects,
                                 const int edgeMode);

//  void UpdateRuptureStates( NodeManager & nodeManager,
//                            EdgeManager & edgeManager,
//                            FaceManager & faceManager,
//                            ElementRegionManager & elementManager,
//                            array1d<set<localIndex>>& nodesToRupturedFaces,
//                            array1d<set<localIndex>>& edgesToRupturedFaces,
//                            const bool prefrac = false );
  void PostUpdateRuptureStates( NodeManager & nodeManager,
                                EdgeManager & edgeManager,
                                FaceManager & faceManager,
                                ElementRegionManager & elementManager,
                                array1d<set<localIndex>>& nodesToRupturedFaces,
                                array1d<set<localIndex>>& edgesToRupturedFaces );

  int CheckEdgeSplitability( const localIndex edgeID,
                                NodeManager & nodeManager,
                                FaceManager & faceManager,
                                EdgeManager & edgeManager,
                                const bool prefrac);

  int CheckNodeSplitability( const localIndex nodeID,
                                NodeManager & nodeManager,
                                FaceManager & faceManager,
                                EdgeManager & edgeManager,
                                const bool prefrac);

  void UpdatePathCheckingArrays();

  bool ProcessNode( const localIndex nodeID,
                    NodeManager & nodeManager,
                    EdgeManager & edgeManager,
                    FaceManager & faceManager,
                    ElementRegionManager & elemManager,
                    array1d<set<localIndex>>& nodesToRupturedFaces,
                    array1d<set<localIndex>>& edgesToRupturedFaces,
                    ElementRegionManager & elementManager,
                    ModifiedObjectLists& modifiedObjects,
                    const bool prefrac );

  bool FindFracturePlanes( const localIndex nodeID,
                           const NodeManager & nodeManager,
                           const EdgeManager & edgeManager,
                           const FaceManager & faceManager,
                           ElementRegionManager & elemManager,
                           const array1d<set<localIndex>>& nodesToRupturedFaces,
                           const array1d<set<localIndex>>& edgesToRupturedFaces,
                           set<localIndex>& separationPathFaces,
                           map<localIndex,int>& edgeLocations,
                           map<localIndex,int>& faceLocations,
                           map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations  );



  void PerformFracture( const localIndex nodeID,
                        NodeManager & nodeManager,
                        EdgeManager & edgeManager,
                        FaceManager & faceManager,
                        ElementRegionManager & elementManager,
                        ModifiedObjectLists& modifiedObjects,
                        array1d<set<localIndex>>& nodesToRupturedFaces,
                        array1d<set<localIndex>>& edgesToRupturedFaces,
                        const set<localIndex>& separationPathFaces,
                        const map<localIndex,int>& edgeLocations,
                        const map<localIndex,int>& faceLocations,
                        const map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );


bool SetLocations( const set<localIndex>& separationPathFaces,
                   ElementRegionManager & elemManager,
                   const FaceManager & faceManager,
                   const set< std::pair<CellBlockSubRegion*,localIndex> >& nodesToElements,
                   const map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                   map<localIndex,int>& edgeLocations,
                   map<localIndex,int>& faceLocations,
                   map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );

bool SetElemLocations( const int side,
                       const std::pair<CellBlockSubRegion*, localIndex >& elem,
                       const set<localIndex>& separationPathFaces,
                       ElementRegionManager & elemManager,
                       const FaceManager & faceManager,
                       const set< std::pair<CellBlockSubRegion*,localIndex> >& nodesToElements,
                       const map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                       map<localIndex,int>& edgeLocations,
                       map<localIndex,int>& faceLocations,
                       map< std::pair<CellBlockSubRegion*, localIndex >, int>& elemLocations );

  realT CalculateKinkAngle (const localIndex nodeID,
                            const NodeManager & nodeManager,
                            EdgeManager & edgeManager,
                            FaceManager & faceManager);

  void CalculateKinkAngles (FaceManager & faceManager,
                            EdgeManager & edgeManager,
                            NodeManager & nodeManager,
                            ModifiedObjectLists& modifiedObjects,
                            const bool prefrac);

  realT MinimumToughnessOnEdge( const localIndex edgeID,
                                const NodeManager & nodeManager,
                                EdgeManager & edgeManager,
                                FaceManager & faceManager );


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto ruptureStateString = "ruptureState";
    constexpr static auto failCriterionString = "failCriterion";
    constexpr static auto degreeFromCrackString = "degreeFromCrack";
  } SurfaceGenViewKeys;


  integer m_failCriterion=1;
  localIndex_set m_separableFaceSet;

  array1d< set<localIndex> > m_originalNodetoFaces;
  array1d< set<localIndex> > m_originalNodetoEdges;

  array1d< array1d<localIndex> > m_originalFaceToEdges;
  array1d< set<localIndex> > m_usedFacesForNode;

  array2d< localIndex > m_originalFacesToElemRegion;
  array2d< localIndex > m_originalFacesToElemSubRegion;
  array2d< localIndex > m_originalFacesToElemIndex;


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_SURFACEGENERATION_SURFACEGENERATOR_HPP_ */
