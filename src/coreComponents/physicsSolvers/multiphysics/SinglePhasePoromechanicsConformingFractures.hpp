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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsSolver.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/contact/LagrangianContactSolver.hpp"

namespace geosx
{

class SinglePhasePoromechanicsConformingFractures : public CoupledSolver< SinglePhasePoromechanicsSolver, LagrangianContactSolver >
{
public:

  using Base = CoupledSolver< SinglePhasePoromechanicsSolver, LagrangianContactSolver  >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Poromechanics = 0,
    Contact = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFractures"; }

  /**
   * @brief main constructor for SinglePhasePoromechanicsConformingFractures objects
   * @param name the name of this instantiation of SinglePhasePoromechanicsConformingFractures in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanicsConformingFractures
   */
  SinglePhasePoromechanicsConformingFractures( const string & name,
                                               Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanicsConformingFractures() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFractures object through the object
   * catalog.
   */
  static string catalogName() { return "SinglePhasePoromechanicsConformingFractures"; }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  LagrangianContactSolver * contactSolver() const
  {
    return std::get< toUnderlying( SolverType::Contact ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the poromechanics solver
   * @return a pointer to the flow solver
   */
  SinglePhasePoromechanicsSolver * poromechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::Poromechanics ) >( m_solvers );
  }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            ParallelVector & rhs,
                            ParallelVector & solution,
                            bool const setSparsity = true ) override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override final;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override final;

  virtual void updateState( DomainPartition & domain ) override final;


  // virtual void implicitStepComplete( real64 const & time_n,
  //                                    real64 const & dt,
  //                                    DomainPartition & domain ) override;

  bool resetConfigurationToDefault( DomainPartition & domain ) const override final;

  bool updateConfiguration( DomainPartition & domain ) override final;

  void initializePostInitialConditionsPostSubGroups() override final;

  /**@}*/

private:

  void assembleForceResidualDerivativeWrtPressure( MeshLevel const & mesh,
                                                   arrayView1d< string const > const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleFluidMassResidualDerivativeWrtDisplacement( MeshLevel const & mesh,
                                                           arrayView1d< string const > const & regionNames,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs );

  void assembleForceResidualDerivativeWrtPressure( MeshLevel & mesh,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void addTransmissibilityCouplingNNZ( DomainPartition const & domain,
                                       DofManager const & dofManager,
                                       arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addTransmissibilityCouplingPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const;

  /**
   * @brief Set up the Dflux_dApertureMatrix object
   *
   * @param domain
   * @param dofManager
   * @param localMatrix
   */
  void setUpDflux_dApertureMatrix( DomainPartition & domain,
                                   DofManager const & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix );

  /**
   * @brief
   *
   * @param domain
   */
  void updateHydraulicAperture( DomainPartition & domain );


  std::unique_ptr< CRSMatrix< real64, localIndex > > & getRefDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture;
  }

  CRSMatrixView< real64, localIndex const > getDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture->toViewConstSizes();
  }

  CRSMatrixView< real64 const, localIndex const > getDerivativeFluxResidual_dAperture() const
  {
    return m_derivativeFluxResidual_dAperture->toViewConst();
  }

  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;

  string const m_pressureKey = SinglePhaseBase::viewKeyStruct::elemDofFieldString();
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_ */
