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
 * @file ExternalFaceManagerT.h
 * @author Scott Johnson
 * @date Jun 15, 2011
 */

#ifndef EXTERNALFACEMANAGERT_H_
#define EXTERNALFACEMANAGERT_H_

#include "Common/Common.h"
#include "EnergyT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "FaceManagerT.h"
#include "DataStructures/VectorFields/StableTimeStep.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "DiscreteElementManagerT.h"

#include "ExternalFaceStructs.h"
#include "Utilities/GeometryUtilities.h"

#include "Contact/SpatialSorterFactory.h"
#include "Constitutive/Interface/InterfaceFactory.h"

#include "IO/ticpp/HierarchicalDataNode.h"

/**
 * @author Scott Johnson
 * @brief Class to manager the collection of external faces as references in a face manager collection
 */
class ExternalFaceManagerT : public ObjectDataStructureBaseT
{
public:
  ExternalFaceManagerT();
  ExternalFaceManagerT(FaceManagerT* fm, FaceManagerT* fmDiscreteElement);
  virtual ~ExternalFaceManagerT();

  void ReadXML(TICPP::HierarchicalDataNode* hdn);

  void Initialize(  ) {}

  virtual void DeserializeObjectField(const std::string& name, const rArray1d& field);
  virtual void DeserializeObjectFields(const sArray1d& names, const Array1dT<rArray1d>& fields);

  void erase( const localIndex i );
  globalIndex resize( const localIndex size, const bool assignGlobals = false );

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject  = NULL) {}
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL) {}
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager ,
                                                         Array1dT<gArray1d>& objectToCompositionObject  )
  {
    (void)compositionObjectManager;
    (void)objectToCompositionObject;
    throw GPException("ExternalFaceManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n");
  }

  void BuildExternalFaces();
  void BuildExternalFacesDiscreteElements();

  void SplitFace(const localIndex parentIndex,
                 const localIndex newFaceIndex,
                 const NodeManagerT& nodeManager);

  void SplitFace( const localIndex parentIndex,
                  const localIndex newFaceIndex1,
                  const localIndex newFaceIndex2,
                  const NodeManagerT& nodeManager);

  void WriteSiloExternalFaces( SiloFile& siloFile,
                               const std::string& siloDirName,
                               const std::string& meshname,
                               const int centering,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart,
                               const std::string& regionName = "none",
                               const lArray1d& mask = lArray1d() );

  void ReadSiloExternalFaces( const SiloFile& siloFile,
                              const std::string& siloDirName,
                              const std::string& meshname,
                              const int centering,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const std::string& regionName = "none",
                              const lArray1d& mask = lArray1d() );
protected:
  virtual globalIndex insert( const localIndex i, const bool assignGlobals = false );

private:
  inline size_t DiscreteElementFaceIndexOffset() const { return this->m_faceManager->DataLengths(); }

  inline void IncrementSearchRadius ( const realT rvel, realT& radius) const
  {
    radius *= (1.0 + this->m_tol.searchRadiusFactor);
    radius += (rvel + this->m_tol.maximumSeparation);
  }

  //FIXME: replace with Geodyn function!!!!
  inline int GeodynPartition(const R1Tensor& ) { return -1; }

public:
  bool RecalculateNeighborList(    const NodeManagerT& m_nodeManager,
                                   const NodeManagerT& m_discreteElementNodeManager,
                                   const DiscreteElementManagerT& m_discreteElementManager,
                                   const realT dt = 0., const bool forceRecalculate = false,
                                   const bool useReferencePosition = false);

#ifdef SRC_INTERNAL
  void GeodynCoupling( NodeManagerT& nodeManager );
  void GeodynCouplingParallel( NodeManagerT& nodeManager );
#endif

  inline bool IsFiniteElement(const localIndex kf) const
  {
    return kf < this->nfe;
  }

  inline localIndex IndexFEFirst() const { return 0;}
  inline localIndex IndexFEAfterLast() const { return this->nfe;}
  inline localIndex IndexDEFirst() const { return this->nfe;}
  inline localIndex IndexDEAfterLast() const { return this->nfe + this->nde;}

  inline localIndex FaceIndex(const localIndex kf, bool& fe) const
  {
    if(!m_faceManager)
      throw GPException("m_faceManager undefined in call to ExternalFaceManager::FaceIndex!");
    fe = IsFiniteElement(kf);
    return this->m_externalFaceToFaceMap[kf] - (fe ? 0 : this->DiscreteElementFaceIndexOffset());
  }
  
  void UpdateGeometricContactProperties(const realT dt,
                                          PhysicalDomainT& domain)
  {
    Array1dT<Array1dT<R1Tensor> > xs;
    xs.resize(this->DataLengths());
    UpdateGeometricContactProperties(dt, domain, xs);
  }

  void UpdateGeometricContactProperties(const realT dt,
                                        PhysicalDomainT& domain,
                                        Array1dT<Array1dT<R1Tensor> >& xs,
                                        const bool updateFESoln = true);

  void SetExcludeFromContact( const NodeManagerT& nodeManager, bool reset );

  void UpdateAndApplyContactStresses( StableTimeStep& maxdt,
                                      const realT dt,
                                      PhysicalDomainT& domain,
                                      const Array1dT<Array1dT<R1Tensor> >& xs);

  inline bool AutoContact() const { return m_autoContact; }

  //@annavarapusr1: Functions to evaluate weighting and stabilization parameter for Nitsche's method
  //                Also evaluate projection tensor to transform between xyz and n-tau1-tau2 planes
  void GetProjectionTensorAndWeightingAndStabilizationParameters( const int dim,
                                                                  const bool planeStress,
                                                                  PhysicalDomainT& domain);
  void GetMeasOmgAndNormD( const int dim,
                           const localIndex FaceID,
                           const PhysicalDomainT& domain,
                           const bool planeStress,
                           realT& measOmg,
                           realT& NormD);

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
private:

  void UpdateGeometricContactPropertiesSub(const realT dt,
                                        PhysicalDomainT& domain,
                                        const Array1dT<Array1dT<R1Tensor> >& xs);

  void SetCommonPlaneGeometryAsOverlap(const localIndex index,
                                       const localIndex kf1,
                                       const localIndex kf2,
                                       PhysicalDomainT& domain,
                                       const Array1dT<Array1dT<R1Tensor> >& xs);

  void SetGhostRank(    const NodeManagerT& m_nodeManager,
                        const NodeManagerT& m_discreteElementNodeManager,
                        const iArray1d& excludeFromContact,
                        iArray1d& faceAttachedToALocalNode);

  void SetGhostRank(    const NodeManagerT& m_nodeManager,
                        const iArray1d& excludeFromContact,
                        iArray1d& faceAttachedToALocalNode);

  void PostSortUpdate( const NodeManagerT& nodeManager,
                       const NodeManagerT& discreteElementNodeManager,
                       const lSet& toResort);

  void PostSortUpdate( const NodeManagerT& nodeManager);

  inline int  MechanicalToOpen() { return opencontact;}

  inline bool IsAnyContact(const int activeC){return activeC != nocontact;}
  inline bool IsDEMContact(const int activeC){return activeC == demcontact;}

  inline int  UpdateContactFlag(const int previousActiveC, const int newActiveC) { return (IsMechanicalContact(previousActiveC) && IsMechanicalContact(newActiveC)) ? persistentmechanicalcontact : newActiveC;}

  inline bool IsMechanicalContact(const int activeC){return activeC == mechanicalcontact || activeC == persistentmechanicalcontact;}
  inline bool IsNewMechanicalContact(const int activeC){return activeC == mechanicalcontact;}

  inline bool IsOpenContact(const int activeC){return activeC == opencontact;}

  /**
   * @brief In parallel, it is advantageous to prevent the ghost faces from contacting each other and to avoid self-contact
   * @author Scott Johnson
   */
  void RemoveInvalidPairsFromNeighborList(const iArray1d& faceAttachedToALocalNode,
                                          const lSet& toResort);

  void RemoveInvalidPairsFromNeighborListSub(const lArray1d& parentElement,
                                             const iArray1d& parentElementRegion,
                                             const iArray1d& faceAttachedToALocalNode,
                                             const localIndex a,
                                             lArray1d& current);

  static int CommonPlaneInterferenceGeometry(const R1Tensor& xfc1,
                                             const R1Tensor& nx1,
                                             const ExternalFaceStruct& efs1,

                                             const R1Tensor& xfc2,
                                             const R1Tensor& nx2,
                                             const ExternalFaceStruct& efs2,

                                             const ToleranceStruct& tol,

                                             R1Tensor& centerCommonPlane,
                                             R1Tensor& applicationPoint,
                                             R1Tensor& normalCommonPlane,
                                             realT& areaCommonPlane,
                                             Array1dT<R1Tensor>& pointsCommonPlane);

  static int CommonEdgeInterferenceGeometry(  const R1Tensor& xfc1,
                                              const R1Tensor& nx1,
                                              const ExternalFaceStruct& efs1,

                                              const R1Tensor& xfc2,
                                              const R1Tensor& nx2,
                                              const ExternalFaceStruct& efs2,

                                              const ToleranceStruct& tol,

                                              R1Tensor& centerCommonEdge,
                                              R1Tensor& normalCommonEdge,
                                              realT& lengthCommonEdge,
                                              Array1dT<R1Tensor>& pointCommonEdge);

  static realT ApplyStress(const R1Tensor& normal,
                           const R1Tensor& applicationPoint,
                           const Array1dT<R1Tensor>& xs,
                           const lArray1d& faceToNodes,
                           const Array1dT<lSet>&  nodeToFaces,
                           const rArray1d& mass,
                           const bool is_fe,
                           const R1Tensor& stress,
                           const realT area,
                           const realT tolParentSolution,
                           R1Tensor& faceParentSolution,
                           Array1dT<R1Tensor>& forces,
                           Array1dT<R1Tensor>& contactForces);

  static void FacePositionAndVelocityFE(const R1Tensor& normal,
                                        const R1Tensor& applicationPoint,
                                        const Array1dT<R1Tensor>& xs,
                                        const Array1dT<R1Tensor>& vs,
                                        const realT tolParentSolution,
                                        R1Tensor& faceParentSolution,
                                        R1Tensor& point,
                                        R1Tensor& velocity);

  static void FacePositionAndVelocityFE(const R1Tensor& normal,
                                        const R1Tensor& applicationPoint,
                                        const Array1dT<R1Tensor>& xs,
                                        const Array1dT<R1Tensor>& vs,
                                        R1Tensor& point,
                                        R1Tensor& velocity);

  static void FacePositionAndVelocityFE(const Array1dT<R1Tensor>& xs,
                                        const Array1dT<R1Tensor>& vs,
                                        R1Tensor& point,
                                        R1Tensor& velocity);

//  static realT CalculateStress(const R1Tensor& normal, const R1Tensor& dx,
//                               const realT normalApproach, realT& stressNormal,
//                               realT& stressDilation, R1Tensor& stressShear, R1Tensor& stress,
//                               R1Tensor& relativeShear,
//                               const realT normalApproachAtDilationInitiation,
//                               const realT normalApproachNormalYield, realT& normalGap,
//                               const realT porePressure, int& ifail,
//                               const realT relativeShearAtDilationInitiation,
//                               const realT relativeShearAtResidualStrength,
//                               const realT relativeShearAtDilationLimit,
//                               const realT normalStressAtDilationLimit,
//                               const realT tangentOfDilationAngleAtZeroStress,
//                               realT& tanFrictionCoefficientInitial,
//                               const realT& tanFrictionCoefficientResidual,
//                               const realT knCoefficientElastic, const realT knPlasticLoading,
//                               const realT knPlasticUnloading, const realT knDilation,
//                               const realT kShear, const realT cohesion, const realT dt);

  static bool Criterion1(const R1Tensor& nx1, const R1Tensor& nx2, const realT tolCosMin);
  static bool Criterion2(const R1Tensor& e1,
                                const R1Tensor& e2,
                                const Array1dT<R1Tensor>& xs1,
                                const Array1dT<R1Tensor>& xs2,
                                Array1dT<R1TensorT<2> >& xsl1,
                                Array1dT<R1TensorT<2> >& xsl2,
                                realT& minDim);
  static bool Criterion3(const Array1dT<R1TensorT<2> >& xsl1,
                                const Array1dT<R1TensorT<2> >& xsl2,
                                Array1dT<R1TensorT<2> >& xsl,
                                realT& area,
                                const realT positionTolerance,
                                const realT areaTolerance);

  void RemoveFromNeighborLists(const localIndex ilocal);
  void RemoveFromNeighborListInverse(const localIndex i0,
                                     const localIndex i1);

  static realT StableTimestep(const realT massFace1, const realT massFace2,
                              const realT area, const realT arealStiffness)
  {
    const realT mass = massFace2 > massFace1 ? massFace1 : massFace2;
    const realT stiffness = area * arealStiffness;
    const realT dttmp = isZero(stiffness, 10*std::numeric_limits<realT>::min()) ?
        std::numeric_limits<realT>::max() :
        0.3 * sqrt(mass/stiffness);
    return dttmp;
  }

public:
  FaceManagerT* m_faceManager;
  FaceManagerT* m_discreteElementFaceManager;

  Array1dT< lArray1d >& m_neighborList;
  Array1dT< lSet >& m_neighborListInverse;

  ///Flag for whether contact is turned on
  bool m_contactActive;

  ///Flag for whether to allow self-contact
  bool m_contactSelf;

  ///Flag for whether contact is turned on
  bool m_autoContact;

  //Flag for whether to use smoothed contact
  bool m_smoothedContact;

  //search and geometry attributes
  ToleranceStruct m_tol;



  //model attributes and model
  InterfaceBase* m_contact;

private:
  enum ContactType
  {
    nocontact            = 0,
    opencontact          = 1,
    mechanicalcontact    = 2,
    persistentmechanicalcontact = 3,
    demcontact = 4
  };

  lArray1d& m_externalFaceToFaceMap;
  localIndex nfe, nde;
  bool m_sorted;
  SpatialSorting::SpatialSorterBase* m_sorter;

};
#endif /* EXTERNALFACEMANAGERT_H_ */
