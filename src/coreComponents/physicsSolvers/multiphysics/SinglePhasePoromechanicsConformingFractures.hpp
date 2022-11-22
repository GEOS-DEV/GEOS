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

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/LagrangianContactSolver.hpp"

namespace geosx
{

class SinglePhasePoromechanicsConformingFractures : public CoupledSolver< LagrangianContactSolver,
                                                                          SinglePhaseBase >
{
public:

  using Base = CoupledSolver< LagrangianContactSolver, SinglePhaseBase >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Contact = 0,
    Flow = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanics"; }

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
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  SinglePhaseBase * flowSolver() const
  {
    return std::get< toUnderlying( SolverType::Flow ) >( m_solvers );
  }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

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
                               arrayView1d< real64 > const & localRhs ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void updateState( DomainPartition & domain ) override;


  virtual void implicitStepComplete( real64 const & time_n,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  /**@}*/

protected:

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * porousMaterialNamesString() { return "porousMaterialNames"; }
  };

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;

private:

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
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_ */
