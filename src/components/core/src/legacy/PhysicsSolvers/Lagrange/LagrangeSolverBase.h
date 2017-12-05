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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file ImplicitMechanicsSolver.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef LAGRANGESOLVERBASE_H_
#define LAGRANGESOLVERBASE_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#include "SurfaceGeneration/FractunatorBase.h"

#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif
/*
 #if GPAC_MPI
 #include "Epetra_MpiComm.h"
 #else
 #include "Epetra_SerialComm.h"
 #endif

 #include "Epetra_Map.h"
 #include "Epetra_FECrsGraph.h"
 #include "Epetra_FECrsMatrix.h"
 #include "Epetra_FEVector.h"
 #include "Epetra_SerialDenseMatrix.h"
 #include "Epetra_SerialDenseVector.h"
 #include "Epetra_IntSerialDenseVector.h"
 #include "Epetra_LinearProblem.h"

 #include "EpetraExt_RowMatrixOut.h"
 #include "Teuchos_RCP.hpp"
 #include "AztecOO.h"
 #include "ml_MultiLevelPreconditioner.h"
 */

class MaterialBaseParameterData;
//class XfemManager;


class LagrangeSolverBase : public SolverBase
{
public:
  LagrangeSolverBase(  const std::string& name,
                       ProblemManagerT* const pm );
  virtual ~LagrangeSolverBase();

  double TimeStep(const realT& time,
                  const realT& dt,
                  const int cycleNumber,
                  PhysicalDomainT& domain,
                  const array<string>& namesOfSolverRegions,
                  SpatialPartition& partition,
                  FractunatorBase* const fractunator);

  void TimeStepXFEM (const realT& time,
                     const realT& dt,
                     const int cycleNumber,
                     PhysicalDomainT& domain,
                     const array<string>& namesOfSolverRegions,
                     SpatialPartition& partition);

//  double TimeStep(const realT& time,
//                const realT& dt,
//                const int cycleNumber,
//                PhysicalDomainT& domain,
//                const array<string>& namesOfSolverRegions,
//                SpatialPartition& partition,
//                FractunatorBase* const fractunator,
//                XfemManager* elementSplitter);


  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeCommunications( PartitionBase& partition );

  virtual void RegisterFields(PhysicalDomainT& domain);

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn );

  void RegisterTemporaryFields( PhysicalDomainT& domain );
  void DeregisterTemporaryFields( PhysicalDomainT& domain );
  void FillTemporaryFields( PhysicalDomainT& domain );
  void OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain );

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const array<string>& namesOfSolverRegions);
  void ApplyTemperatureBC(PhysicalDomainT& domain, realT time);

  enum TimeIntegrationEnum
  {
    FirstTimeIntegrationEnum = 0,
    QuasiStatic = 0,
    ImplicitDynamic = 1,
    ExplicitDynamic = 2,
    numTimeIntegrationEnums = 3
  };

  enum TwoDOptions
  {
    PlaneStrain = 0,
    PlaneStress = 1,
    Axisym = 2,
    num2dOptions = 3
  };

  TimeIntegrationEnum IntToTimeIntegration( const int ival )
  {
    switch( ival )
    {
    case 0:
      return QuasiStatic;
    case 1:
      return ImplicitDynamic;
    case 2:
      return ExplicitDynamic;
    }
    throw GPException("Invalid value input into LagrangeSolverBase::IntToTimeIntegration()");
    return numTimeIntegrationEnums;
  }

  TwoDOptions IntToTwoDOptions( const int ival )
  {
    switch( ival )
    {
    case 0:
      return PlaneStrain;
    case 1:
      return PlaneStress;
    case 2:
      return Axisym;
    }
    throw GPException("Invalid value input into LagrangeSolverBase::IntToTwoDOptions()");
    return num2dOptions;
  }

  void SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                     SpatialPartition& partition,
                                     int& numLocalRows,
                                     int& numGlobalRows,
                                     array<integer>& localIndices,
                                     int offset );

  virtual void SetSparsityPattern( PhysicalDomainT& domain );

  void SetSparsityPatternXFEM ( PhysicalDomainT& domain);

  typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;



  virtual realT Assemble    ( PhysicalDomainT& domain,
                              const SpatialPartition& partition,
                              const realT time,
                              realT const dt );

  void EvaluateDispAtVertex(PhysicalDomainT&  domain,
                            const SpatialPartition& partition,
                            const realT time );

  void GetNodalForces( PhysicalDomainT&  domain,
                       const SpatialPartition& partition,
                       const realT dt );

  void GetInverseOfA3By3Mat(const rArray2d& mat,
                            rArray2d& invMat);

  void MultiplyArray(const rArray2d& A,
                     const array<real64>& B,
                     array<real64>& C);

  void CalculateNodalForces ( NodeManager& nodeManager,
                              ElementRegionT& elemRegion,
                              const localIndex elementID,
                              R2SymTensor& stress);

  virtual void PostProcessFieldsForVisualizationAndConsistency( PhysicalDomainT& domain)
  {
    (void)domain;
  }

  virtual void StoreHistoryVariablesForCurrentLoadStepAndResetTheField( PhysicalDomainT& domain)
  {
    (void)domain;
  }


  /*
     realT AssembleFluidPressureContributions( PhysicalDomainT& domain,
                                           const array<integer>&
                                              deformationTrilinosIndex,
                                           const array<integer>&
                                              flowTrilinosIndex,
                                           const array< array<real64> >& dwdu,
                                           const array<real64>& dwdw,
                                           const int flowDofOffset );
   */


//  std::string m_trilinosIndexStr;
  static int m_instances;

  realT DampingM() const { return m_dampingM; }
  realT DampingK() const { return m_dampingK; }



  int m_enableTimeIntegrationOption[numTimeIntegrationEnums];
  TimeIntegrationEnum m_timeIntegrationOption;

protected:

  virtual void InsertGlobalIndices( PhysicalDomainT& domain)
  {
    (void)domain;
  }

  virtual void UpdateContactDataStructures( PhysicalDomainT& domain,
                                            const bool setActiveInit )
  {
    (void)domain;
  }

  virtual void GetContactStiffnessContribution( PhysicalDomainT& domain)
  {
    (void)domain;
  }

  int dim;
  const unsigned this_mpi_process;
  const unsigned n_mpi_processes;
  const bool     verbose;


public:
  realT m_dampingM;
  realT m_dampingK;
  realT m_newmarkBeta;
  realT m_newmarkGamma;


  realT m_timeToSnapshotDisp;
  // Fu: If we specify an initial stress field, the stress field itself might
  // not be balanced.
  // We need to store the displacement after the model has settled.  We use the
  // difference between the total disp and this ref disp for the Displace
  // transformation in VisIt.
  bool m_writeNodalStress;

  realT m_refTemperature;
  realT m_defaultCTE;
  int m_useNodalTemperature;
  array<string> m_thermalRegionNames;
  int m_relax = 0;
//  int m_reassembleMatrix = 1;
//  array< std::tuple< localIndex,localIndex,realT > > m_originalMatrixValues;

  virtual void TimeStepImplicitSetup( const realT& time,
                                      const realT& dt,
                                      PhysicalDomainT& domain );

  virtual void TimeStepImplicitComplete( const realT& time,
                                         const realT& dt,
                                         PhysicalDomainT& domain );

  void CalculateAllNodalMassAndVolume( PhysicalDomainT& domain, SpatialPartition& partition );

protected:
  int m_tiedNodesFlag;
  array<lSet> m_KinematicConstraintNodes;
  realT m_tiedNodeNormalRuptureStress;
  realT m_tiedNodeShearRuptureStress;
  realT m_tiedNodeTolerance;

  array<integer> dummyDof;


  realT m_cfl;


  TwoDOptions m_2dOption;

  realT m_bulkQLinear;
  realT m_bulkQQuadratic;

public:
  R1Tensor m_gravityVector;

  int m_useVelocityQS = 0;
  realT m_displacementPenalty = 0.0;

private:
  virtual void SetupSystem ( PhysicalDomainT& domain, SpatialPartition& partition );

  virtual void TimeStepExplicitDynamic( const realT& time,
                                        const realT& dt,
                                        PhysicalDomainT& domain,
                                        const array<string>& namesOfSolverRegions,
                                        SpatialPartition& partition );



  virtual realT TimeStepImplicit( const realT& time,
                                  const realT& dt,
                                  PhysicalDomainT& domain,
                                  const array<string>& namesOfSolverRegions,
                                  SpatialPartition& partition,
                                  FractunatorBase* const fractunator );

public:  // Fu note: We need to access this from the new hydrofrac solver, which
         // is not a derived class.
  void ApplyForcesFromContact(PhysicalDomainT& domain,
                              StableTimeStep& timeStep,
                              const realT dt );


  virtual void ApplyThermalStress( ElementRegionT& elemRegion,
                                   NodeManager& nodeManager,
                                   const localIndex& elementID,
                                   Epetra_SerialDenseVector * rhs)
  {
    (void)elemRegion;
    (void)nodeManager;
    (void)elementID;
    (void)rhs;
  };


  virtual void ProcessElementRegions( NodeManager& nodeManager,
                                      ElementManagerT& elemManager,
                                      const array<string>& namesOfSolverRegions,
                                      const realT dt  );

  void SnapshotNodalDisplacement( NodeManager& nodeManager);

  void ApplyDisplacementPenalty( NodeManager& nodeManager );

  virtual void FixNodesBC( NodeManager const& nodeManager, const lSet& set )
  {
    (void)(nodeManager);
    (void)(set);
  }

private:
  virtual void ProcessElementRegion( NodeManager& nodeManager,
                                     ElementRegionT& elemRegion,
                                     const realT dt ) = 0;
  virtual void ApplyGravity( NodeManager& nodeManager,
                             ElementRegionT& elemRegion,
                             const realT dt );

  virtual void ProcessCohesiveZones( NodeManager& nodeManager,
                                     FaceManagerT& faceManager,
                                     const realT dt );

  //virtual void RegisterFields_Derived(PhysicalDomainT& domain);
  //virtual void ReadXML_Derived(TICPP::HierarchicalDataNode* hdn);


  virtual realT CalculateElementResidualAndDerivative( const MaterialBaseParameterData& matParams,
                                                       const FiniteElementBase& fe,
                                                       const Array2dT<R1Tensor>& dNdX,
                                                       const realT* const detJ,
                                                       R2SymTensor const * const refStress,
                                                       array<R1Tensor> const & u,
                                                       array<R1Tensor> const & uhat,
                                                       array<R1Tensor> const & uhattilde,
                                                       array<R1Tensor> const & vtilde,
                                                       realT const dt,
                                                       Epetra_SerialDenseMatrix& dRdU,
                                                       Epetra_SerialDenseVector& R) = 0;

  virtual void IntegrateCohesiveZoneContributions( const int numNodesInFace,
                                                   const realT area,
                                                   const R1Tensor& traction,
                                                   const R2Tensor& stiffness,
                                                   Epetra_SerialDenseMatrix& face_matrix,
                                                   Epetra_SerialDenseVector& face_rhs );



private:
  virtual void SetupMLPreconditioner( const PhysicalDomainT& domain, ML_Epetra::MultiLevelPreconditioner* MLPrec );


  // boundary conditions
  virtual void TractionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time)
  {
    (void)domain;
    (void)object;
    (void)bc;
    (void)set;
    (void)time;
  }

  virtual void PressureBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time)
  {
    (void)domain;
    (void)object;
    (void)bc;
    (void)set;
    (void)time;
  }

  virtual void DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                              BoundaryConditionBase* bc, const lSet& set, realT time)
  {
    (void)domain;
    (void)object;
    (void)bc;
    (void)set;
    (void)time;
  }

};



#endif /* ImplicitMechanicsSolver_H_ */
