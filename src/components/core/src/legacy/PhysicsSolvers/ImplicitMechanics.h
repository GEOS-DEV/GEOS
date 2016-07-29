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

#ifndef ImplicitMechanicsSolver_H_
#define ImplicitMechanicsSolver_H_

#include "SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif

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




template <int dim>
class ImplicitMechanicsSolver: public SolverBase
{
public:
  ImplicitMechanicsSolver( const std::string& name,
                           ProblemManagerT* const pm );
  ~ImplicitMechanicsSolver();

  double
  TimeStep(const realT& time, const realT& dt,
           const int cycleNumber,
           PhysicalDomainT& domain,
           const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
           FractunatorBase* const fractunator);


  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition  );

  void InitializeCommunications( PartitionBase& partition );

  virtual void
  RegisterFields(PhysicalDomainT& domain);

  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn );

  static const char*
  SolverName()
  {
    return "ImplicitMechanicsSolver3D";
  };
  
  // boundary conditions
  virtual void TractionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                          BoundaryConditionBase* bc, const lSet& set, realT time);

  virtual void PressureBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time);
  
  virtual void DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                              BoundaryConditionBase* bc, const lSet& set, realT time);
                              
  
  
  virtual void NonpenetratingBC_NeighborUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                                               BoundaryConditionBase* bc, realT time);
  virtual void NonpenetratingBC_DetectContact(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                                              BoundaryConditionBase* bc, realT time);
  
  virtual void NonpenetratingBC_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                                         BoundaryConditionBase* bc, realT time);
                                         
  virtual void NonpenetratingBC_Apply(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                                      BoundaryConditionBase* bc, realT time);
                                      
  virtual void NonpenetratingBC_UpdateAperture(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
                                               BoundaryConditionBase* bc, realT time);
                      

private:

  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition);
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, realT time);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition, realT time);

  #if GPAC_MPI
    const Epetra_MpiComm & epetra_comm;
  #else
    const Epetra_SerialComm & epetra_comm;
  #endif

  const unsigned this_mpi_process;
  const unsigned n_mpi_processes;
  const bool     verbose;

  bool m_useMLPrecond;

  struct
  {
    double lambda; // Lame's constant
    double G;      // Shear modulus
  }
  equation_data;

  struct
  {
    double krylov_tol;
  }
  numerics;

  typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;

  Teuchos::RCP<Epetra_Map>         row_map;
  Teuchos::RCP<Epetra_FECrsGraph>  sparsity;
  Teuchos::RCP<Epetra_FECrsMatrix> matrix;
  Teuchos::RCP<Epetra_FEVector>    solution;
  Teuchos::RCP<Epetra_FEVector>    rhs;
  iArray1d dummyDof;
  
  std::string m_trilinosIndexStr;
  static int m_instances;


  realT m_cfl;
  realT m_nonContactModulus;
  bool m_doApertureUpdate;
  
  bool m_recordIncrementalDisplacement;
  
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
};

template <int dim>
int ImplicitMechanicsSolver<dim>::m_instances = 0;


#endif /* ImplicitMechanicsSolver_H_ */
