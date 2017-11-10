/**
 * @file Fractunator.h
 * @author settgast1
 * @date Jul 14, 2011
 */

#ifndef FRACTUNATORBASE_H_
#define FRACTUNATORBASE_H_

#include "Common/Common.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"


class SpatialPartition;
class ModifiedObjectLists;

class FractunatorBase
{
public:
  FractunatorBase();
  virtual ~FractunatorBase();

  virtual void Initialize( NodeManager& nodeManager,
                           EdgeManagerT& edgeManager,
                           FaceManagerT& faceManager,
                           ElementManagerT& elementManager);

  virtual void RegisterFieldsAndMaps( NodeManager& nodeManager,
                              EdgeManagerT& edgeManager,
                              FaceManagerT& faceManager );

  virtual int SeparationDriver( NodeManager& nodeManager,
                                EdgeManagerT& edgeManager,
                                FaceManagerT& faceManager,
                                ExternalFaceManagerT& externalFaceManager,
                                ElementManagerT& elementManager,
                                SpatialPartition& partition,
                                const bool prefrac,
                                const realT time) = 0;

  virtual void ReadXML( TICPP::HierarchicalDataNode& hdn );

  virtual void WriteSilo( SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart );

  virtual void ReadSilo( const SiloFile& siloFile,
                         const int cycleNum,
                         const realT problemTime,
                         const bool isRestart );

  virtual void PreexistingFracture2D( NodeManager& nodeManager,
                                      EdgeManagerT& edgeManager,
                                      FaceManagerT& faceManager,
                                      ExternalFaceManagerT& externalFaceManager,
                                      ElementManagerT& elementManager,
                                      SpatialPartition& partition,
                                      const bool prefrac)=0;





  int m_verbose;
  realT m_failstress;
  std::string m_separableNodeSet;
  std::string m_separableFaceSet;
  realT m_failgap;
  FaceManagerT m_virtualFaces;
  int m_checkInterval;
  realT m_rockToughness;
  realT m_maxTurnAngle;
  int m_failCriterion;  // =0 stress based; = 1 SIF based; =2 mixed.
  int m_allowVacuumFrac;  // If = 0, only allow fracture to propagate when the tip cell is saturated
  bool m_displacementBasedSIF;
  realT m_kJn;


  virtual void UpdateCohesiveZones( NodeManager& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager ){}


protected:
  virtual void UpdateRuptureStates( NodeManager& nodeManager,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager,
                            ElementManagerT& elementManager,
                            array<lSet>& nodesToRupturedFaces,
                            array<lSet>& edgesToRupturedFaces,
                            const bool prefrac = false );

  virtual bool ProcessNode( const localIndex nodeID,
                            NodeManager& nodeManager,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager,
                            ExternalFaceManagerT& externalFaceManager,
                            array<lSet>& nodesToRupturedFaces,
                            array<lSet>& edgesToRupturedFaces,
                            ElementManagerT& elementManager,
                            ModifiedObjectLists& modifiedObjects,
                            const bool prefrac);

  int CheckOrphanElement (FaceManagerT& faceManager,
                        localIndex iFace);

  void MarkBirthTime( FaceManagerT& faceManager,
                      ModifiedObjectLists& modifiedObjects,
                      const realT time);

  void CorrectSplitNodalMass (NodeManager& nodeManager,
                              localIndex node0,
                              localIndex node1);

  //The two functions are for 3D and 2D, respectively.
  realT MinimumToughnessOnEdge( const localIndex edgeID,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager);
  realT MinimumToughnessOnNode( const localIndex nodeID,
                                NodeManager& nodeManager,
                                FaceManagerT& faceManager);
  void MarkDiscreteFractureNetworkFaces(FaceManagerT& faceManager,
                                        SpatialPartition& partition);



  bool m_tounessSetByInitialCondition;
  bool m_failStressSetByInitialCondition;
  realT m_thresholdForEvaluateFace;

  std::string m_dfnPrefix;

private:


  virtual int CheckNodeSplitability( const localIndex nodeID,
                              NodeManager& nodeManager,
                              FaceManagerT& faceManager,
                              EdgeManagerT& edgeManager,
                              const bool prefrac) = 0;

  virtual int CheckEdgeSplitability( const localIndex edgeID,
                                 NodeManager& nodeManager,
                                 FaceManagerT& faceManager,
                                 EdgeManagerT& edgeManager,
                                 const bool prefrac) = 0;


  virtual bool FindFracturePlanes( const localIndex nodeID,
                                   const NodeManager& nodeManager,
                                   const EdgeManagerT& edgeManager,
                                   const FaceManagerT& faceManager,
                                   const array<lSet>& nodesToRupturedFaces,
                                   const array<lSet>& edgesToRupturedFaces,
                                   lSet& separationPathFaces,
                                   std::map<localIndex,int>& edgeLocations,
                                   std::map<localIndex,int>& faceLocations,
                                   std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations  ) = 0;



  virtual void PerformFracture( const localIndex nodeID,
                                NodeManager& nodeManager,
                                EdgeManagerT& edgeManager,
                                FaceManagerT& faceManager,
                                ExternalFaceManagerT& externalFaceManager,
                                ElementManagerT& elementManager,
                                ModifiedObjectLists& modifiedObjects,
                                array<lSet>& nodesToRupturedFaces,
                                array<lSet>& edgesToRupturedFaces,
                                const lSet& separationPathFaces,
                                const std::map<localIndex,int>& edgeLocations,
                                const std::map<localIndex,int>& faceLocations,
                                const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations ) = 0;



  virtual void WriteSiloDerived( SiloFile& siloFile,
                                 const int cycleNum,
                                 const realT problemTime,
                                 const bool isRestart ){}

  virtual void ReadSiloDerived( const SiloFile& siloFile,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart ){}


};





#endif /* FRACTUNATORBASE_H_ */
