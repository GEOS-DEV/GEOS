/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FracturedPoroelasticSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_POROELASTICSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_POROELASTICSOLVER_HPP_

#include "physicsSolvers/multiphysics/PoroelasticSolver.hpp"

namespace geosx
{

class SolidMechanicsEmbeddedFractures;

class FracturedPoroelasticSolver : public PoroelasticSolver
{
public:
  FracturedPoroelasticSolver( const std::string & name,
                     Group * const parent );
  ~FracturedPoroelasticSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "FracturedPoroelastic"; }

  virtual void SetupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            array1d< real64 > & localRhs,
                            array1d< real64 > & localSolution,
                            bool const setSparsity = true ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  AssembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  void
  AssembleCouplingTerms( DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs );

  virtual void
  ApplyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  void UpdateDeformationForCoupling( DomainPartition & domain );

  real64 SplitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto fracturesSolverNameString = "fracturesSolverName";
  } poroElasticSolverViewKeys;


  SolidMechanicsEmbeddedFractures * getFracturesSolver()
  {
    return this->getParent()->GetGroup( m_fracturesSolverName )->group_cast< SolidMechanicsEmbeddedFractures * >();
  }
  SolidMechanicsEmbeddedFractures const * getFracturesSolver() const
  {
    return this->getParent()->GetGroup( m_fracturesSolverName )->group_cast< SolidMechanicsEmbeddedFractures const * >();
  }

 protected:

  virtual void PostProcessInput() override final;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const problemManager ) override final;

private:

  void CreatePreconditioner();

  string m_fracturesSolverName;

  SolidMechanicsEmbeddedFractures * m_fracturesSolver;

};

ENUM_STRINGS( PoroelasticSolver::CouplingTypeOption, "FIM", "SIM_FixedStress" )

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_POROELASTICSOLVER_HPP_ */
