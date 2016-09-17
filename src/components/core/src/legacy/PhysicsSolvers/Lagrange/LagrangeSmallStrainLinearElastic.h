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
 * @file ImplicitMechanicsSolver.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef LAGRANGESMALLSTRAINLINEARELASTIC_H_
#define LAGRANGESMALLSTRAINLINEARELASTIC_H_

#include "LagrangeSmallStrain.h"





class LagrangeSmallStrainLinearElastic: public LagrangeSmallStrain
{
public:
  LagrangeSmallStrainLinearElastic(  const std::string& name,
                                     ProblemManagerT* const pm );
  ~LagrangeSmallStrainLinearElastic();


  static const char* SolverName()
  {
    return "LagrangeSmallStrainLinearElastic";
  }


  virtual void PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition );

//  virtual void PostProcess (PhysicalDomainT& domain,
//                    SpatialPartition& partition,
//                    const sArray1d& namesOfSolverRegions);

protected:
  //   @annavarapusr1: Functions used in Penalty and Nitsche stiffness calculations

  struct TrialTractions { Array1dT<rArray1d> tracPar, tracSib;};
  struct UpdatedTractions { Array1dT<rArray1d> tracPar, tracSib;};
  struct InitialModulii { Array1dT<rArray2d> nParDotDParBPar, nSibDotDSibBSib, nParDotDParBParPermute, nSibDotDSibBSibPermute;
                          rArray2d alphaPar, alphaSib;};
  struct UpdatedModulii { Array1dT<rArray2d> nParDotDParBPar, nSibDotDSibBSib, alphaParTimesNSub, alphaSibTimesNSub;};

  virtual void GetContactFacePairIndices ( const PhysicalDomainT& domain,
                                           const localIndex i,
                                           localIndex& kf1,
                                           localIndex& kf2);

  virtual void InsertGlobalIndices( PhysicalDomainT& domain);

  virtual void UpdateContactDataStructures( PhysicalDomainT& domain,
                                            const bool setActiveInit );

  virtual void GetContactStiffnessContribution( PhysicalDomainT& domain);

  bool IsContactActive( const PhysicalDomainT& domain,
                        const localIndex& iContFace);

  virtual void GetRowDofIndexColDofIndex( const localIndex i,
                                          const PhysicalDomainT& domain,
                                          const bool nitsche_flag,
                                          iArray1d& rowDofIndex,
                                          iArray1d& colDofIndex);

  virtual void GetRowDofIndexColDofIndex( const localIndex i,
                                          const PhysicalDomainT& domain,
                                          const bool nitsche_flag,
                                          Epetra_IntSerialDenseVector& rowDofIndex,
                                          Epetra_IntSerialDenseVector& colDofIndex);

  virtual void GetNodIndicesNotOnFace( const PhysicalDomainT& domain,
                                       const localIndex kf,
                                       const bool globFlag,
                                       iArray1d& NodIDNotOnFace);

  virtual void GetParentFaceNodesAndCoordsAndInterfaceGaussPoints( const localIndex i,
                                                                   const PhysicalDomainT& domain,
                                                                   const rArray1d& gauss,
                                                                   lArray1d& localParentFaceNodes,
                                                                   rArray1d& xe,
                                                                   rArray1d& psi,
                                                                   rArray1d& eta);

  virtual void GetLocalIndexOnInterface( const localIndex iC,
                                         const localIndex a,
                                         const PhysicalDomainT& domain,
                                         const lArray1d& localParentFaceNodes,
                                         localIndex& edgeLocalIndex,
                                         localIndex& sibLocalIndex);

  virtual void GetJacobianAndShapeFunctionsOnInterface( const localIndex numNodes,
                                                        const rArray1d& xe,
                                                        const realT& psi,
                                                        const realT& eta,
                                                        realT& jcob_sub,
                                                        rArray2d& N_sub);

  virtual void MultiplyArray( const rArray2d& A,
                              const rArray2d& B,
                              rArray2d& C);

  virtual void MultiplyArray( const rArray2d& A,
                              const rArray1d& B,
                              rArray1d& C);

  virtual void GetFaceConnToElemConnMap( const PhysicalDomainT& domain,
                                         const localIndex kf,
                                         const bool sibFlag,
                                         iArray1d& FaceConnToElemConn);

  virtual void GetNormalVoigt( const PhysicalDomainT& domain,
                               const localIndex kf,
                               rArray2d& normalVoigt);

  virtual void GetElasticityTensorVoigt( const PhysicalDomainT& domain,
                                         const localIndex kf,
                                         rArray2d& D);

  virtual void GetNitscheSpecificMatrices ( const PhysicalDomainT& domain,
                                            const ElementRegionT* elemRegion,
                                            const rArray2d& D,
                                            rArray2d& B,
                                            const rArray2d& normalVoigt,
                                            const rArray2d& P,
                                            const iArray1d& FaceConnToElemConn,
                                            const realT psi,
                                            const realT eta,
                                            const localIndex kf,
                                            const localIndex numNodesEle,
                                            const realT gamma,
                                            rArray2d& nDotDB,
                                            rArray2d& nDotDBPermute);

  virtual void GetShapeFunctionDerivativeMatrixConstStr( const PhysicalDomainT& domain,
                                                         const iArray1d& FaceConnToElemConn,
                                                         const localIndex kf,
                                                         rArray2d& B);

  virtual void GetShapeFunctionDerivativeMatrixQuadHex( const PhysicalDomainT& domain,
                                                        const iArray1d& FaceConnToElemConn,
                                                        const localIndex kf,
                                                        realT psi,
                                                        realT eta,
                                                        rArray2d& B);

  virtual void GetDBDotN( const rArray2d& D,
                          const rArray2d& B,
                          const rArray2d& normalVoigt,
                          const localIndex numNodesParentEle,
                          const realT gamma,
                          const PhysicalDomainT& domain,
                          const rArray2d& P,
                          rArray2d& nDotDB);

  virtual void GetPermutedDBDotN( const PhysicalDomainT& domain,
                                  const rArray2d& nDotDB,
                                  const iArray1d& FaceConnToElemConn,
                                  const localIndex kf,
                                  rArray2d& nDotDBPermute);

  virtual void GetPermutedNodalVector( const PhysicalDomainT& domain,
                                       const rArray1d& Vector,
                                       const iArray1d& FaceConnToElemConn,
                                       const localIndex kf,
                                       rArray1d& PermutedVector);

  virtual void GetInitialAlphaTensor( const PhysicalDomainT& domain,
                                      const localIndex kf,
                                      const rArray2d& P,
                                      rArray2d& alpha_XYZ);

  virtual void RotateMatrix( const rArray2d& Q,
                             const rArray2d& A,
                             rArray2d& APrime);

  virtual void RotateTensor( const rArray2d& Q,
                             const rArray2d& A,
                             rArray2d& APrime,
                             const bool transposeFlag);

  virtual void GetTransformationTensor( const PhysicalDomainT& domain,
                                        const localIndex kf,
                                        rArray2d& P);

  virtual void GetDisplacementEleNodes( const PhysicalDomainT& domain,
                                        const localIndex kf,
                                        rArray1d& uEle);

  virtual void GetDisplacementJumpInterface( const PhysicalDomainT& domain,
                                             const ElementRegionT& elemRegion,
                                             const localIndex kf1,
                                             const localIndex kf2,
                                             const rArray2d& N_sub,
                                             const localIndex iContFace,
                                             const lArray1d& localParentFaceNodes,
                                             rArray1d& uJump);

  virtual void GetTrialAndUpdatedTractions ( PhysicalDomainT& domain,
                                             const rArray2d& N_sub,
                                             const InitialModulii& initialModulus,
                                             const rArray2d& P_par,
                                             const rArray2d& P_sib,
                                             const lArray1d& localParentFaceNodes,
                                             const localIndex iGp,
                                             const localIndex kf1,
                                             const localIndex kf2,
                                             const localIndex iContFace,
                                             const localIndex numNodesParentEle,
                                             const localIndex numNodesSiblingEle,
                                             const std::string s_previous_base,
                                             const std::string s_current_base,
                                             iArray1d& stickGP_par,
                                             iArray1d& stickGP_sib,
                                             iArray1d& openingGP_par,
                                             iArray1d& openingGP_sib,
                                             TrialTractions& trialTractions,
                                             UpdatedTractions& updatedTractions);

  virtual void GetTractionInterface( const PhysicalDomainT& domain,
                                     const rArray2d& N_sub,
                                     const rArray2d& alphaPar,
                                     const rArray2d& alphaSib,
                                     const rArray2d& nParDotDParBParPermute,
                                     const rArray2d& nSibDotDSibBSibPermute,
                                     const lArray1d& localParentFaceNodes,
                                     const localIndex kf1,
                                     const localIndex kf2,
                                     const localIndex iContFace,
                                     rArray1d& tracPar,
                                     rArray1d& tracSib);

  virtual void GetTractionInterfaceSlidingTrial( const PhysicalDomainT& domain,
                                                 const rArray2d& N_sub,
                                                 const rArray2d& alphaPar,
                                                 const rArray2d& alphaSib,
                                                 const rArray2d& P_par,
                                                 const rArray2d& P_sib,
                                                 const rArray2d& nParDotDParBParPermute,
                                                 const rArray2d& nSibDotDSibBSibPermute,
                                                 const lArray1d& localParentFaceNodes,
                                                 const localIndex kf1,
                                                 const localIndex kf2,
                                                 const localIndex iContFace,
                                                 const std::string FieldName,
                                                 rArray1d& tracPar,
                                                 rArray1d& tracSib);

  virtual void GetTractionInterfaceSliding( const rArray2d& N_sub,
                                            const rArray2d& alphaPar,
                                            const rArray2d& alphaSib,
                                            const rArray2d& P_par,
                                            const rArray2d& P_sib,
                                            const rArray2d& nParDotDParBParPermute,
                                            const rArray2d& nSibDotDSibBSibPermute,
                                            const lArray1d& localParentFaceNodes,
                                            const localIndex kf1,
                                            const localIndex kf2,
                                            const localIndex iContFace,
                                            const int iGp,
                                            const std::string s_previous_base,
                                            const std::string s_current_base,
                                            PhysicalDomainT& domain,
                                            iArray1d& stickGP_par,
                                            iArray1d& stickGP_sib,
                                            iArray1d& openingGP_par,
                                            iArray1d& openingGP_sib,
                                            rArray1d& tracPar,
                                            rArray1d& tracSib,
                                            Array1dT<rArray1d>& tracParTrial,
                                            Array1dT<rArray1d>& tracSibTrial);

  virtual void GetTrialTangentModulus( const PhysicalDomainT& domain,
                                         const localIndex kf,
                                         const rArray2d& P,
                                         const iArray1d& FaceConnToElemConn,
                                         const realT psi,
                                         const realT eta,
                                         const rArray2d& D,
                                         const rArray2d& normalVoigt,
                                         const realT gamma,
                                         rArray2d& alpha,
                                         rArray2d& nDotDB);

  virtual void GetUpdatedModulusAtGaussPoints (const PhysicalDomainT& domain,
                                               const rArray2d& N_sub,
                                               const rArray2d& P,
                                               const iArray1d& FaceConnToElemConn,
                                               const localIndex kf,
                                               const localIndex iGp,
                                               const localIndex numNodesEle,
                                               const int openingGP,
                                               const int stickGP,
                                               const bool parFlag,
                                               const InitialModulii& initialModulus,
                                               const TrialTractions& trialTractions,
                                               UpdatedModulii& updatedModulus);

  virtual void UpdateTangentModulusForOpeningAndSliding(const PhysicalDomainT& domain,
                                                        const rArray2d& P,
                                                        const iArray1d& FaceConnToElemConn,
                                                        const localIndex kf,
                                                        rArray1d& trac,
                                                        const int openingGP,
                                                        const int stickGP,
                                                        const bool parFlag,
                                                        rArray2d& nDotDB,
                                                        rArray2d& alpha,
                                                        rArray2d& nDotDBPermute);

  virtual void GetUpdatedTangentModulus( const PhysicalDomainT& domain,
                                         const rArray2d& P,
                                         const iArray1d& FaceConnToElemConn,
                                         const localIndex kf,
                                         rArray1d& trac,
                                         rArray2d& alpha,
                                         rArray2d& nDotDB,
                                         const bool parFlag);

  virtual realT GetSlipStickState( const PhysicalDomainT& domain,
                                   const rArray1d& trac,
                                   const rArray2d& P,
                                   int& stick);

  virtual void GetNegativeArray( rArray1d& myVec);

  virtual void GetMatrixTranspose (const rArray2d& myMat,
                                   rArray2d& myMatTransposed);

  virtual void GetUpdatedStiffnessPenalty ( const PhysicalDomainT& domain,
                                            const rArray2d& P,
                                            rArray1d& trac,
                                            rArray2d& alphaXYZ,
                                            const bool parFlag);

  virtual void GetUpdatedStiffnessNitsche ( const PhysicalDomainT& domain,
                                            const rArray2d& P,
                                            rArray1d& trac,
                                            rArray2d& nDotDB);

  virtual void GetUpdatedTractionsAndPlasticSlip( const localIndex kf,
                                                  const rArray2d& P,
                                                  const realT& phiTrial,
                                                  const int stick,
                                                  const std::string FieldNamePrevious,
                                                  const std::string FieldNameCurrent,
                                                  const rArray2d& alphaXYZ,
                                                  PhysicalDomainT& domain,
                                                  rArray1d& tracXYZ);

  virtual void GetUpdatedTractionsOpeningMode( const rArray2d& P,
                                               const localIndex iGp,
                                               const realT tol,
                                               iArray1d& openingGP,
                                               rArray1d& trac);

  virtual void PostProcessFieldsForVisualizationAndConsistency( PhysicalDomainT& domain);

  virtual void StoreHistoryVariablesForCurrentLoadStepAndResetTheField( PhysicalDomainT& domain);

  void ApplyThermalStress( ElementRegionT& elemRegion,
                           NodeManager& nodeManager,
                           const localIndex& elementID,
                           Epetra_SerialDenseVector * rhs);

private:

  virtual void ProcessElementRegion( NodeManager& nodeManager,
                                     ElementRegionT& elemRegion,
                                     const realT dt );


  virtual realT CalculateElementResidualAndDerivative( const MaterialBaseParameterData& matParams,
                                                       const FiniteElementBase& fe,
                                                       const Array2dT<R1Tensor>& dNdX,
                                                       const realT* const detJ,
                                                       R2SymTensor const * const refStress,
                                                       Array1dT<R1Tensor> const & u,
                                                       Array1dT<R1Tensor> const & uhat,
                                                       Array1dT<R1Tensor> const & uhattilde,
                                                       Array1dT<R1Tensor> const & vtilde,
                                                       realT const dt,
                                                       Epetra_SerialDenseMatrix& dRdU,
                                                       Epetra_SerialDenseVector& R );

  void CalculateElementStresses( const NodeManager& nodeManager,
                                 ElementManagerT& elementManager );

};


#endif /* LAGRANGESMALLSTRAINLINEARELASTIC_H_ */
