/**
 * @file Fractunator.h
 * @author settgast1
 * @date Jul 14, 2011
 */

#ifndef FRACTUNATOR_H_
#define FRACTUNATOR_H_

#include "Common/Common.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"


class Fractunator
{
public:
  Fractunator();
  virtual ~Fractunator();

  void Initialize( const NodeManagerT& nodeManager ,
                   const EdgeManagerT& edgeManager ,
                   const FaceManagerT& faceManager ,
                   const ElementManagerT& elementManager );


  void RegisterFieldsAndMaps( NodeManagerT& nodeManager,
                              EdgeManagerT& edgeManager,
                              FaceManagerT& faceManager );

  void SeparationDriver( NodeManagerT& nodeManager,
                         EdgeManagerT& edgeManager,
                         FaceManagerT& faceManager,
                         ExternalFaceManagerT& externalFaceManager,
                         ElementManagerT& elementManager);

  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  int m_verbose;
  realT m_failstress;

private:
  void UpdateRuptureStates( NodeManagerT& nodeManager,
                            EdgeManagerT& edgeManager,
                            FaceManagerT& faceManager,
                            ElementManagerT& elementManager );

  void UpdatePathCheckingArrays();

  bool FindFracturePlanes( const localIndex nodeID,
                           const NodeManagerT& nodeManager,
                           const EdgeManagerT& edgeManager,
                           const FaceManagerT& faceManager,
                           lSet& separationPathFaces,
                           std::map<localIndex,int>& edgeLocations,
                           std::map<localIndex,int>& faceLocations,
                           std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations  );



  void PerformFracture( const localIndex nodeID,
                        NodeManagerT& nodeManager,
                        EdgeManagerT& edgeManager,
                        FaceManagerT& faceManager,
                        ElementManagerT& elementManager,
                        const lSet& separationPathFaces,
                        const std::map<localIndex,int>& edgeLocations,
                        const std::map<localIndex,int>& faceLocations,
                        const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations );

bool SetLocations( const int side,
                   const lSet& separationPathFaces,
                   const FaceManagerT& faceManager,
                   const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                   std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                   std::map<localIndex,int>& edgeLocations,
                   std::map<localIndex,int>& faceLocations,
                   std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations );

  void ApplyGapDamping( NodeManagerT& nodeManager,
                        const FaceManagerT& faceManager,
                        const realT dt  );
};



#endif /* FRACTUNATOR_H_ */
