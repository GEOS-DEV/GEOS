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
 * @file FlowSolverBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

/**
 * @class FlowSolverBase
 *
 * Base class for finite volume fluid flow solvers.
 * Provides some common features
 */
class FlowSolverBase : public SolverBase
{
public:

  /// String used to form the solverName used to register single-physics solvers in CoupledSolver
  static string coupledSolverAttributePrefix() { return "flow"; }

/**
 * @brief main constructor for Group Objects
 * @param name the name of this instantiation of Group in the repository
 * @param parent the parent group of this instantiation of Group
 */
  FlowSolverBase( const string & name,
                  Group * const parent );


  /// deleted default constructor
  FlowSolverBase() = delete;

  /// deleted copy constructor
  FlowSolverBase( FlowSolverBase const & ) = delete;

  /// default move constructor
  FlowSolverBase( FlowSolverBase && ) = default;

  /// deleted assignment operator
  FlowSolverBase & operator=( FlowSolverBase const & ) = delete;

  /// deleted move operator
  FlowSolverBase & operator=( FlowSolverBase && ) = delete;

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  localIndex numDofPerCell() const { return m_numDofPerCell; }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // input data
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * permeabilityString() { return "permeability"; }

    // misc inputs
    static constexpr char const * fluidNamesString() { return "fluidNames"; }
    static constexpr char const * solidNamesString() { return "solidNames"; }
    static constexpr char const * permeabilityNamesString() { return "permeabilityNames"; }
    static constexpr char const * isThermalString() { return "isThermal"; }
    static constexpr char const * solidInternalEnergyNamesString() { return "solidInternalEnergyNames"; }
    static constexpr char const * inputFluxEstimateString() { return "inputFluxEstimate"; }

  };

  void updatePorosityAndPermeability( CellElementSubRegion & subRegion ) const;

  virtual void updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion ) const;


  /**
   * @brief For each equilibrium initial condition, loop over all the target cells and compute the min/max elevation
   * @param[in] domain the domain partition
   * @param[in] equilNameToEquilId the map from the name of the initial condition to the initial condition index (used in min/maxElevation)
   * @param[out] maxElevation the max elevation for each initial condition
   * @param[out] minElevation the min elevation for each initial condition
   */
  void findMinMaxElevationInEquilibriumTarget( DomainPartition & domain, // cannot be const...
                                               std::map< string, localIndex > const & equilNameToEquilId,
                                               arrayView1d< real64 > const & maxElevation,
                                               arrayView1d< real64 > const & minElevation ) const;

  /**
   * @brief For each source flux boundary condition, loop over all the target cells and sum the owned cells
   * @param[in] time the time at the beginning of the time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   * @param[in] bcNameToBcId the map from the name of the boundary condition to the boundary condition index
   * @param[out] bcAllSetsSize the total number of owned cells for each source flux boundary condition
   */
  void computeSourceFluxSizeScalingFactor( real64 const & time,
                                           real64 const & dt,
                                           DomainPartition & domain, // cannot be const...
                                           std::map< string, localIndex > const & bcNameToBcId,
                                           arrayView1d< globalIndex > const & bcAllSetsSize ) const;


protected:

  /**
   * @brief Increment the cumulative flux from each aquifer
   * @param[in] time the time at the beginning of the time step
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   *
   * For now this function is here because it can be used for both single-phase flow and multiphase flow
   * This may have to be revisited when aquifer BC is implemented for hybrid FVM
   */
  virtual void saveAquiferConvergedState( real64 const & time,
                                          real64 const & dt,
                                          DomainPartition & domain );

  virtual void precomputeData( MeshLevel & mesh,
                               arrayView1d< string const > const & regionNames );

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const override;

  /// the number of Degrees of Freedom per cell
  integer m_numDofPerCell;

  /// flag to determine whether or not this is a thermal simulation
  integer m_isThermal;

  real64 m_fluxEstimate;


private:
  virtual void setConstitutiveNames( ElementSubRegionBase & subRegion ) const override;


};


}

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_
