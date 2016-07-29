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
 * @file ParallelPlateFlowSolverFV.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVER_IM_H_
#define PARALLELPLATEFLOWSOLVER_IM_H_

#include "../ParallelPlateFlowSolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Common/typedefs.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "../PhysicsSolverStrings.h"

#include <set>


#ifdef SRC_EXTERNAL
#include "BoundaryConditions/PerforatedCasedWellboreBoundaryCondition.h"
#include "BoundaryConditions/CavityPressureBoundaryCondition.h"
#endif

class ParallelPlateFlowSolverFV : public ParallelPlateFlowSolverBase
{
public:
  ParallelPlateFlowSolverFV(  const std::string& name,
                              ProblemManagerT* const pm );
  virtual ~ParallelPlateFlowSolverFV();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  void RegisterTemporaryFields( PhysicalDomainT& domain );
  void DeregisterTemporaryFields( PhysicalDomainT& domain );
  void FillTemporaryFields( PhysicalDomainT& domain );
  void OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain );

  virtual void TimeStepSetup( const realT& time,
                              const realT& dt,
                              PhysicalDomainT& domain,
                              SpatialPartition& partition,
                              const bool setupSystem  );

  virtual double TimeStepExecute( const realT& time,
                          const realT& dt,
                          PhysicalDomainT& domain,
                          SpatialPartition& partition );

  virtual void TimeStepCleanup( PhysicalDomainT& domain, const realT& dt );

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  void DefineFlowSets( PhysicalDomainT& domain );

  void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt );
  void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain );
  void UpdateFlux( const realT time,
                   const realT dt,
                   PhysicalDomainT& domain, SpatialPartition& partition);
  
  // Pressure controlled boundary condition
  void PressureBoundaryCondition(PhysicalDomainT& domain,
                                 ObjectDataStructureBaseT& object,
                                 BoundaryConditionBase* const bc,
                                 const lSet& set,
                                 const realT time,
                                 const realT dt,
                                 const int dofOffset );
                                  
  void PressureBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain,
                                                ObjectDataStructureBaseT& object ,
                                                BoundaryConditionBase* const bc,
                                                const lSet& set,
                                                const realT time,
                                                const realT dt);

  void FaceFluxBoundaryCondition( PhysicalDomainT& domain,ObjectDataStructureBaseT& object,
      BoundaryConditionBase* const bc, const lSet& aset, const realT time, const realT dt );

  virtual void MassBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                      BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void MassRateBC( PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                           BoundaryConditionBase* bc, const lSet& set, realT time, realT dt);

  virtual void VolumeRateBC( PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                             BoundaryConditionBase* bc, const lSet& set, realT time, realT dt);

  void SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName);

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ParallelPlateFlowSolverFV";};

  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,realT time,realT dt );

  // Boundary conditions
//  realT m_dt; // current dt

  // Flags
  bool m_doApertureUpdate;

  void SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                     SpatialPartition& partition,
                                     int& numLocalRows,
                                     int& numGlobalRows,
                                     iArray1d& localIndices,
                                     int offset );


  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt)
  {
    throw GPException("ParallelPlateFlowSolverBase::Assemble() not overridden");
  }


  realT Assemble    ( PhysicalDomainT& domain,
                     Epetra_System& epetraSystem,
                     const realT& time,
                     const realT& dt );

  realT AssembleCD( PhysicalDomainT& domain,
                    Epetra_System& epetraSystem,
                    const realT& time,
                    const realT& dt );

  using SolverBase::SetInitialGuess;
  virtual void SetInitialGuess( const PhysicalDomainT& domain,
                                realT* const local_solution );

  virtual realT CheckSolution( const realT* const local_solution,
                               const PhysicalDomainT& domain,
                               const localIndex dofOffset );

  virtual void PropagateSolution( const realT* const local_solution,
                                  const realT scalingFactor,
                                  PhysicalDomainT& domain,
                                  const localIndex dofOffset  );

  virtual void PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition );


/*
  void Solve ( PhysicalDomainT&  domain,
               SpatialPartition& partition,
               const realT time,
               const realT dt );
*/


  virtual void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain );


  virtual void CalculateCarterLeakOff( const realT time, const realT dt, PhysicalDomainT& domain );

  virtual void ApplyFluxBoundaryCondition( const realT time, const realT dt, const int cycleNumber, const int rank, PhysicalDomainT& domain );

  virtual void FlowControlledBoundaryCondition( PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
          BoundaryConditionBase* bc ,
          const lSet& set,
          realT time );


  virtual void CalculateApertureDerivatives( const FaceManagerT& faceManager,
                                             const NodeManagerT& nodeManager );
private:


  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeDensity( PhysicalDomainT& domain);
  void UpdateAperture(PhysicalDomainT&  domain);

  realT TwoFacePermeability_FV(const Array1dT<R1Tensor>& edgeCenters,
                            const rArray1d& edgeLengths,
                            const Array1dT<R1Tensor>& faceCenters,
                            const rArray1d& apertures,
                            const localIndex eg,
                            const localIndex kf,
                            const localIndex kfb,
                            const Array1dT<rArray1d>* const dwdu,
                            rArray1d* const dkdu_r,
                            rArray1d* const dkdu_s);
  
  lSet m_faceSet;
  localIndex m_numFaces;
  std::map<localIndex,localIndex> m_faceDofMap;
  
  realT m_negativePressureSlope;

  
  realT m_phi; // Mixed Euler parameter (0 = forward difference, 0.5 = central difference, 1.0 = backward difference)
  

  
  // MPI
  const int this_mpi_process;
  const int n_mpi_processes;
  
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
  
  static int m_instances;
  realT m_gamma;
  
  Array1dT<locallyIndexedReal> indexedFluxes;

  
};


#endif /* PARALLELPLATEFLOWSOLVER_IM_H_ */
