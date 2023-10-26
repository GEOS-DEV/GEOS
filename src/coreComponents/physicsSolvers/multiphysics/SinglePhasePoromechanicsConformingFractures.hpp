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

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_

#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/contact/LagrangianContactSolver.hpp"

namespace geos
{

class SinglePhasePoromechanicsConformingFractures : public CoupledSolver< SinglePhasePoromechanics, LagrangianContactSolver >
{
public:

  using Base = CoupledSolver< SinglePhasePoromechanics, LagrangianContactSolver >;
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
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override  { return catalogName(); }

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
  SinglePhasePoromechanics * poromechanicsSolver() const
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

  virtual void updateState( DomainPartition & domain ) override final;


  // virtual void implicitStepComplete( real64 const & time_n,
  //                                    real64 const & dt,
  //                                    DomainPartition & domain ) override;

  bool resetConfigurationToDefault( DomainPartition & domain ) const override final;

  bool updateConfiguration( DomainPartition & domain ) override final;

  void initializePostInitialConditionsPostSubGroups() override final;

  void outputConfigurationStatistics( DomainPartition const & domain ) const override final;

  /**@}*/

  struct viewKeyStruct : Base::viewKeyStruct
  {
    /// Flag to indicate that the simulation is thermal
    constexpr static char const * isThermalString() { return "isThermal"; }
  };

private:

  void assembleCellBasedContributions( real64 const time_n,
                                       real64 const dt,
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs );

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override final;

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

  /**
   * @brief
   *
   * @param domain
   */
  void updateHydraulicApertureAndFracturePermeability( DomainPartition & domain );


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

  /// flag to determine whether or not this is a thermal simulation
  integer m_isThermal;
};

template< typename CONSTITUTIVE_BASE,
          typename KERNEL_WRAPPER,
          typename ... PARAMS >
real64 SinglePhasePoromechanicsConformingFractures::assemblyLaunch( MeshLevel & mesh,
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
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  KERNEL_WRAPPER kernelWrapper( dispDofNumber,
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
                                                                 contactSolver()->getSolidSolver()->getDiscretizationName(),
                                                                 materialNamesString,
                                                                 kernelWrapper );
}


} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP_ */
