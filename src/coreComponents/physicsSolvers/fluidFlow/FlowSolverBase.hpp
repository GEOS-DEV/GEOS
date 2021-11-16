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

  /**
   * @brief default destructor
   */
  virtual ~FlowSolverBase() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  void setPoroElasticCoupling() { m_poroElasticFlag = 1; }

  void setReservoirWellsCoupling() { m_coupledWellsFlag = 1; }

  arrayView1d< string const > fluidModelNames() const { return m_fluidModelNames; }

  arrayView1d< string const > permeabilityModelNames() const { return m_permeabilityModelNames; }

  virtual std::vector< string > getConstitutiveRelations( string const & regionName ) const override;


  localIndex numDofPerCell() const { return m_numDofPerCell; }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // input data
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * permeabilityString() { return "permeability"; }

    // gravity term precomputed values
    static constexpr char const * gravityCoefString() { return "gravityCoefficient"; }

    // misc inputs
    static constexpr char const * fluidNamesString() { return "fluidNames"; }
    static constexpr char const * solidNamesString() { return "solidNames"; }
    static constexpr char const * permeabilityNamesString() { return "permeabilityNames"; }
    static constexpr char const * pressureString() { return "pressure"; }
    static constexpr char const * deltaPressureString() { return "deltaPressure"; }
    static constexpr char const * initialPressureString() { return "initialPressure"; }
    static constexpr char const * deltaVolumeString() { return "deltaVolume"; }
    static constexpr char const * aperture0String() { return "aperture_n"; }
    static constexpr char const * hydraulicApertureString() { return "hydraulicAperture"; }
    static constexpr char const * conductivity0String() { return "conductivity0"; }
    static constexpr char const * inputFluxEstimateString() { return "inputFluxEstimate"; }
  };

  void updatePorosityAndPermeability( CellElementSubRegion & subRegion,
                                      localIndex const targetIndex ) const;

  virtual void updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion,
                                              localIndex const targetIndex ) const;

  /**
   * @brief Setup stored views into domain data for the current step
   * @param[in] mesh the mesh level object
   */
  virtual void resetViews( MeshLevel & mesh );

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

  virtual void precomputeData( MeshLevel & mesh );

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /// name of the fluid constitutive model
  array1d< string > m_fluidModelNames;

  /// name of the solid constitutive model
  array1d< string > m_solidModelNames;

  /// name of the permeability constituive model
  array1d< string > m_permeabilityModelNames;

  /// flag to determine whether or not coupled with solid solver
  integer m_poroElasticFlag;

  /// flag to determine whether or not coupled with wells
  integer m_coupledWellsFlag;

  /// the number of Degrees of Freedom per cell
  integer m_numDofPerCell;

  real64 m_fluxEstimate;

  /// views into pressure fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaPressure;

  /// views into constant data fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > m_elemGhostRank;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_volume;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_gravCoef;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >  m_permeability;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >  m_dPerm_dPressure;

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  m_elementSeparationCoefficient;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  m_element_dSeparationCoefficient_dAperture;
#endif

};


}

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_
