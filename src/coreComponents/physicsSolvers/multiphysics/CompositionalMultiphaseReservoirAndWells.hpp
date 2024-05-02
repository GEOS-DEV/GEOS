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
 * @file CompositionalMultiphaseReservoirAndWells.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASERESERVOIRANDWELLS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASERESERVOIRANDWELLS_HPP_

#include "physicsSolvers/multiphysics/CoupledReservoirAndWellsBase.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"

namespace geos
{

template< typename COMPOSITIONAL_RESERVOIR_SOLVER >
class CompositionalMultiphaseReservoirAndWells : public CoupledReservoirAndWellsBase< COMPOSITIONAL_RESERVOIR_SOLVER,
                                                                                      CompositionalMultiphaseWell >
{
public:

  using Base = CoupledReservoirAndWellsBase< COMPOSITIONAL_RESERVOIR_SOLVER,
                                             CompositionalMultiphaseWell >;
  using Base::getLogLevel;
  using Base::m_solvers;
  using Base::m_linearSolverParameters;

  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  CompositionalMultiphaseReservoirAndWells( const string & name,
                                            dataRepository::Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseReservoirAndWells() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName();
  /**
   * @copydoc SolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void addCouplingSparsityPattern( DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           SparsityPatternView< globalIndex > const & pattern ) const override;

  virtual void assembleCouplingTerms( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition const & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs ) override;

  void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const
  { flowSolver()->assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );  }
  void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const
  { flowSolver()->assembleStabilizedFluxTerms( dt, domain, dofManager, localMatrix, localRhs );  }

  void keepFlowVariablesConstantDuringInitStep( bool const keepFlowVariablesConstantDuringInitStep )
  { flowSolver()->keepFlowVariablesConstantDuringInitStep( keepFlowVariablesConstantDuringInitStep ); }

  real64 updateFluidState( ElementSubRegionBase & subRegion ) const
  { return flowSolver()->updateFluidState( subRegion ); }
  void updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const
  { flowSolver()->updatePorosityAndPermeability( subRegion ); }
  void updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const
  { flowSolver()->updateSolidInternalEnergyModel( dataGroup ); }

  integer & isThermal() { return flowSolver()->isThermal(); }
  integer useSimpleAccumulation() const { return flowSolver()->useSimpleAccumulation(); }
  integer useTotalMassEquation() const { return flowSolver()->useTotalMassEquation(); }
  integer numFluidPhases() { return flowSolver()->numFluidPhases(); }
  integer numFluidComponents() { return flowSolver()->numFluidComponents(); }

  void enableFixedStressPoromechanicsUpdate()
  { flowSolver()->enableFixedStressPoromechanicsUpdate();  }

  void enableJumpStabilization()
  {flowSolver()->enableJumpStabilization(); }

  virtual void saveSequentialIterationState( DomainPartition & domain ) override final { flowSolver()->saveSequentialIterationState( domain ); }

protected:

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  CompositionalMultiphaseBase * flowSolver() const;

  void setMGRStrategy();

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COMPOSITIONALMULTIPHASERESERVOIRANDWELLS_HPP_ */
