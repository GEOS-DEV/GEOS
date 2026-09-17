/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsConformingFracturesALM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_

#include "common/logger/Logger.hpp"
// #include "physicsSolvers/multiphysics/MultiphasePoromechanicsConformingFracturesALM.hpp"
#include "physicsSolvers/solidMechanics/contact/SolidMechanicsAugmentedLagrangianContact.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"

namespace geos
{

template< typename FLOW_SOLVER = SinglePhaseBase >
class SinglePhasePoromechanicsConformingFracturesALM : public PoromechanicsConformingFractures<  SinglePhasePoromechanics, FLOW_SOLVER, SolidMechanicsAugmentedLagrangianContact >
{
public:

  using Base = PoromechanicsConformingFractures< SinglePhasePoromechanics, FLOW_SOLVER , SolidMechanicsAugmentedLagrangianContact >;
  using Base::m_solvers;
  using Base::m_dofManager;
  using Base::m_localMatrix;
  using Base::m_rhs;
  using Base::m_solution;
  using Base::m_maxFaceNodes;

  using Base::m_derivativeFluxResidual_dAperture;
  using Base::m_derivativeFluxResidual_dApertureOffsets;

  /// True when the flow solver carries well degrees of freedom.
  static constexpr bool hasWells = std::is_same_v< FLOW_SOLVER, SinglePhaseReservoirAndWells<> >;

  static_assert( hasWells || std::is_same_v< FLOW_SOLVER, SinglePhaseBase >,
                 "SinglePhasePoromechanicsConformingFracturesALM supports only the SinglePhaseBase and "
                 "SinglePhaseReservoirAndWells<> flow solvers. Both setMGRStrategy and assembleSystem branch "
                 "on hasWells, so a new instantiation must be handled in both places." );

  /// String used to form the solverName used to register solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "poromechanicsConformingFracturesALM"; }

  /**
   * @brief main constructor for SinglePhasePoromechanicsConformingFracturesALM objects
   * @param name the name of this instantiation of SinglePhasePoromechanicsConformingFracturesALM in the repository
   * @param parent the parent group of this instantiation of SinglePhasePoromechanicsConformingFracturesALM
   */
  SinglePhasePoromechanicsConformingFracturesALM( const string & name,
                                                  dataRepository::Group * const parent );

  /// Destructor for the class
  ~SinglePhasePoromechanicsConformingFracturesALM() override {}

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new SinglePhasePoromechanicsConformingFracturesALM object through the
   * object
   * catalog.
   */
  static string catalogName()
  {
    if constexpr ( std::is_same_v< FLOW_SOLVER, SinglePhaseBase > )
    {
      return "SinglePhasePoromechanicsConformingFracturesALM";
    }
    else
    {
      return FLOW_SOLVER::catalogName() + "PoromechanicsConformingFracturesALM";
    }
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/


  virtual void setSparsityPattern( DomainPartition & domain,
                                   DofManager & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                   SparsityPattern< globalIndex > & pattern ) override final;


  virtual void updateState( DomainPartition & domain ) override final;

  virtual void setMGRStrategy() override final
  {
    LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();
    if( linearSolverParameters.preconditionerType != LinearSolverParameters::PreconditionerType::mgr )
    {
      return;
    }

    // Wells contribute their own dof labels and need an extra reduction level
    // to keep the well block out of the coarse grid, so they get a separate
    // strategy.
    if (this->m_isThermal)
    {
      if( this->m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::mgr )
      GEOS_ERROR( GEOS_FMT( "{}: MGR strategy is not implemented for {}", this->getName(), this->getCatalogName() ) );
    }

    if constexpr ( hasWells )
    {
      linearSolverParameters.mgr.strategy =
        LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFracturesALMReservoirFVM;
    }
    else
    {
      linearSolverParameters.mgr.strategy =
        LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanicsConformingFracturesALM;
    }
  
    linearSolverParameters.mgr.separateComponents = true;

    GEOS_LOG_LEVEL_RANK_0( logInfo::LinearSolver,
                           GEOS_FMT( "{}: MGR strategy set to {}", this->getName(),
                                     EnumStrings< LinearSolverParameters::MGR::StrategyType >::toString( linearSolverParameters.mgr.strategy ) ) );
  }

  /**@}*/

protected:

  virtual void initializePreSubGroups() override
  {
    Base::initializePreSubGroups();
  }

private:

  struct viewKeyStruct : public Base::viewKeyStruct
  {};


 

  void assembleForceResidualDerivativeWrtPressure( string const & meshName,
                                                   MeshLevel const & mesh,
                                                   string_array const & regionNames,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleFluidMassResidualDerivativeWrtDisplacement( string const & meshName,
                                                           MeshLevel const & mesh,
                                                           string_array const & regionNames,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs ) override final;




  void assembleMatrixPressureBubbleContribution( real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs ) override;


  void updateHydraulicApertureAndFracturePermeability( DomainPartition & domain );

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURESALM_HPP_ */
