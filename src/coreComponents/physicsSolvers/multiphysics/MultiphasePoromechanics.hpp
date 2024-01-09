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
 * @file MultiphasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsSolver.hpp"


namespace geos
{

namespace stabilization
{
enum class StabilizationType : integer
{
  None,
  Global,
  Local,
};

ENUM_STRINGS( StabilizationType,
              "None",
              "Global",
              "Local" );
}

template< typename FLOW_SOLVER >
class MultiphasePoromechanics : public PoromechanicsSolver< FLOW_SOLVER >
{
public:

  using Base = PoromechanicsSolver< FLOW_SOLVER >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;

  /**
   * @brief main constructor for MultiphasePoromechanics Objects
   * @param name the name of this instantiation of MultiphasePoromechanics in the repository
   * @param parent the parent group of this instantiation of MultiphasePoromechanics
   */
  MultiphasePoromechanics( const string & name,
                           dataRepository::Group * const parent );

  /// Destructor for the class
  ~MultiphasePoromechanics() override {};

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new MultiphasePoromechanics object through the object catalog.
   */
  static string catalogName();
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  SolidMechanicsLagrangianFEM * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( Base::SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  FLOW_SOLVER * flowSolver() const
  {
    return std::get< toUnderlying( Base::SolverType::Flow ) >( m_solvers );
  }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void postProcessInput() override;

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;

  virtual void setupCoupling( DomainPartition const & domain,
                              DofManager & dofManager ) const override;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void updateState( DomainPartition & domain ) override;

  /**@}*/

  /*
   * @brief Utility function to update the stabilization parameters at each time step
   * @param[in] domain the domain partition
   */
  void updateStabilizationParameters( DomainPartition & domain ) const;

  /*
   * @brief Utility function to set the stress initialization flag
   * @param[in] performStressInitialization true if the solver has to initialize stress, false otherwise
   */
  void setStressInitialization( integer const performStressInitialization )
  { m_performStressInitialization = performStressInitialization; }

  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType ) override final;

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;

  struct viewKeyStruct : Base::viewKeyStruct
  {
    /// Names of the porous materials
    constexpr static char const * porousMaterialNamesString() { return "porousMaterialNames"; }

    /// Type of stabilization used in the simulation
    constexpr static char const * stabilizationTypeString() { return "stabilizationType"; }

    /// Names of the regions where the stabilization is applied
    constexpr static char const * stabilizationRegionNamesString() { return "stabilizationRegionNames"; }

    /// Multiplier on stabilization
    constexpr static char const * stabilizationMultiplierString() { return "stabilizationMultiplier"; }

    /// Flag to determine whether or not this is aa thermal simulation
    constexpr static char const * isThermalString() { return "isThermal"; }

    /// Flag to indicate that the solver is going to perform stress initialization
    constexpr static char const * performStressInitializationString() { return "performStressInitialization"; }

  };

private:

  /**
   * @brief Helper function to recompute the bulk density
   * @param[in] subRegion the element subRegion
   */
  void updateBulkDensity( ElementSubRegionBase & subRegion );

  /**
   * @brief Helper function to average the mean total stress increment over quadrature points
   * @param[in] domain the domain partition
   */
  void averageMeanTotalStressIncrement( DomainPartition & domain );

  template< typename CONSTITUTIVE_BASE,
            typename KERNEL_WRAPPER,
            typename ... PARAMS >
  real64 assemblyLaunch( MeshLevel & mesh,
                         DofManager const & dofManager,
                         string_array const & regionNames,
                         string const & materialNamesString,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         real64 const dt,
                         PARAMS && ... params );

  /// Type of stabilization used in the simulation
  stabilization::StabilizationType m_stabilizationType;

  /// Names of the regions where the stabilization is applied
  string_array m_stabilizationRegionNames;

  /// Multiplier on stabilization constant
  real64 m_stabilizationMultiplier;

  /// flag to determine whether or not this is a thermal simulation
  integer m_isThermal;

  /// Flag to indicate that the solver is going to perform stress initialization
  integer m_performStressInitialization;

};

template< typename FLOW_SOLVER >
template< typename CONSTITUTIVE_BASE,
          typename KERNEL_WRAPPER,
          typename ... PARAMS >
real64 MultiphasePoromechanics< FLOW_SOLVER >::assemblyLaunch( MeshLevel & mesh,
                                                               DofManager const & dofManager,
                                                               string_array const & regionNames,
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

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( SolverBase::gravityVector() );

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

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_HPP_ */
