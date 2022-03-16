/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiResolutionHFSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIRESOLUTIONHFSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIRESOLUTIONHFSOLVER_HPP_

#include "physicsSolvers/multiphysics/PhaseFieldFractureSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsEmbeddedFractures.hpp"

namespace geosx
{

class SurfaceGenerator;


class MultiResolutionHFSolver : public SolverBase
{
public:
  MultiResolutionHFSolver( const string & name,
                       Group * const parent );

  ~MultiResolutionHFSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "MultiResolutionHF";
  }

  virtual void RegisterDataOnMesh( Group & MeshBodies );

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void applyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void updateNodeMap()

  virtual void setInitialCrackDamageBCs()

  virtual void prepareSubProblemBCs()
    
  virtual void setNextDt( real64 const & currentDt,
                          real64 & nextDt ) override;


  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override final;

  void updateDeformationForCoupling( DomainPartition & domain );

  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  void initializeNewFaceElements( DomainPartition const & domain );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * globalSolverNameString() { return "globalSolverName"; }

    constexpr static char const * localSolverNameString() { return "localSolverName"; }
        
    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * surfaceGeneratorNameString() { return "surfaceGeneratorName"; }

    constexpr static char const * maxNumResolvesString() { return "maxNumResolves"; }

  };

protected:
  virtual void postProcessInput() override final;

  virtual void
  initializePostInitialConditionsPreSubGroups() override final;

private:

  // name of the global efem solver
  string m_globalSolverName;

  //name of local phase-field fracture solver
  string m_localSolverName;

  // name of the contact relation
  string m_contactRelationName;

  /// name of the surface generator
  string m_surfaceGeneratorName;

  // list of damage dofs to be fixed in the subdomain boundary
  Array1d<Real64> m_dofListDamage;

  // list of disp dofs to be fixed in the subdomain boundary
  Array1d<Real64> m_dofListDisp;
  
  // list of displacement values to be prescribed in the subdomain boundary
  Array1d<Real64> m_fixedDispList;
  
  // map that translate subproblem's node numbers to outer problem's node numbers
  Array1d<int> m_nodeMap;

  // pointer to global efem solver
  SolidMechanicsEmbeddedFractures * m_globalSolver;
    
  // pointer to local phase-field fracture solver
  PhaseFieldFractureSolver * m_localSolver;

  /// pointer to the surface generator
  SurfaceGenerator * m_surfaceGenerator;

  std::unique_ptr< ParallelMatrix > m_blockDiagUU;

  integer m_maxNumResolves;
  integer m_numResolves[2];

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_MULTIRESOLUTIONHFSOLVER_HPP_ */
