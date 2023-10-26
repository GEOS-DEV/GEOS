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
 * @file SinglePhasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICS_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geos
{

class SinglePhasePoromechanics : public CoupledSolver< SinglePhaseBase,
                                                       SolidMechanicsLagrangianFEM >
{
public:

  using Base = CoupledSolver< SinglePhaseBase, SolidMechanicsLagrangianFEM >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  enum class SolverType : integer
  {
    Flow = 0,
    SolidMechanics = 1
  };

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanics"; }

  /**
   * @brief main constructor for SinglePhasePoromechanics objects
   * @param name the name of this instantiation of SinglePhasePoromechanics in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanics
   */
  SinglePhasePoromechanics( const string & name,
                            Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanics() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanics object through the object catalog.
   */
  static string catalogName() { return "SinglePhasePoromechanics"; }
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override  { return catalogName(); }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
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

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override;

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

  virtual void updateState( DomainPartition & domain ) override;

  /*
   * @brief Utility function to set the stress initialization flag
   * @param[in] performStressInitialization true if the solver has to initialize stress, false otherwise
   */
  void setStressInitialization( integer const performStressInitialization )
  { m_performStressInitialization = performStressInitialization; }

  /**@}*/

  virtual void mapSolutionBetweenSolvers( DomainPartition & Domain, integer const idx ) override final;

  struct viewKeyStruct : Base::viewKeyStruct
  {
    /// Names of the porous materials
    constexpr static char const * porousMaterialNamesString() { return "porousMaterialNames"; }

    /// Flag to indicate that the simulation is thermal
    constexpr static char const * isThermalString() { return "isThermal"; }

    /// Flag to indicate that the solver is going to perform stress initialization
    constexpr static char const * performStressInitializationString() { return "performStressInitialization"; }
  };

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;

  void assembleElementBasedTerms( real64 const time_n,
                                  real64 const dt,
                                  DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );
  /// flag to determine whether or not this is a thermal simulation
  integer m_isThermal;

private:

  /**
   * @brief Helper function to recompute the bulk density
   * @param[in] subRegion the element subRegion
   */
  void updateBulkDensity( ElementSubRegionBase & subRegion );

  /**
   * @brief Helper function to average the mean stress increment
   * @param[in] domain the domain partition
   */
  void averageMeanStressIncrement( DomainPartition & domain );

  void createPreconditioner();

  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  real64 assemblyLaunch( MeshLevel & mesh,
                         DofManager const & dofManager,
                         arrayView1d< string const > const & regionNames,
                         string const & materialNamesString,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         real64 const dt,
                         PARAMS && ... params );

  /// Flag to indicate that the solver is going to perform stress initialization
  integer m_performStressInitialization;
};

template< typename CONSTITUTIVE_BASE,
          typename KERNEL_WRAPPER,
          typename ... PARAMS >
real64 SinglePhasePoromechanics::assemblyLaunch( MeshLevel & mesh,
                                                 DofManager const & dofManager,
                                                 arrayView1d< string const > const & regionNames,
                                                 string const & materialNamesString,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs,
                                                 real64 const dt,
                                                 PARAMS && ... params )
{
  GEOS_MARK_FUNCTION;

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const dofKey = dofManager.getKey( fields::solidMechanics::totalDisplacement::key() );
  arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  KERNEL_WRAPPER kernelWrapper( dofNumber,
                                dofManager.rankOffset(),
                                localMatrix,
                                localRhs,
                                dt,
                                gravityVectorData,
                                std::forward< PARAMS >( params )... );

  return finiteElement::
           regionBasedKernelApplication< parallelDevicePolicy< >,
                                         CONSTITUTIVE_BASE,
                                         CellElementSubRegion >( mesh,
                                                                 regionNames,
                                                                 solidMechanicsSolver()->getDiscretizationName(),
                                                                 materialNamesString,
                                                                 kernelWrapper );
}


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICS_HPP_ */
