/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseProppantBase.hpp"

namespace geosx
{


/**
 * @class SinglePhaseFVM
 *
 * class to perform a single phase finite volume solve
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
template< typename BASE = SinglePhaseBase >
class SinglePhaseFVM : public BASE
{
public:


  // Aliasing public/protected members/methods of Group so we don't
  // have to use this->member etc.
  using BASE::getLogLevel;

  // Aliasing public/protected members/methods of SolverBase so we don't
  // have to use this->member etc.
  using BASE::m_systemSolverParameters;
  using BASE::m_cflFactor;
  using BASE::m_maxStableDt;
  using BASE::m_nextDt;
  using BASE::m_discretizationName;
  using BASE::targetRegionNames;
  using BASE::forTargetRegions;
  using BASE::forTargetRegionsComplete;
  using BASE::forTargetSubRegions;
  using BASE::forTargetSubRegionsComplete;
  using BASE::m_dofManager;
  using BASE::m_matrix;
  using BASE::m_rhs;
  using BASE::m_solution;
  using BASE::m_linearSolverParameters;
  using BASE::m_nonlinearSolverParameters;

  // Aliasing public/protected members/methods of FlowSolverBase so we don't
  // have to use this->member etc.
  using BASE::m_fluidModelNames;
  using BASE::m_solidModelNames;
  using BASE::m_poroElasticFlag;
  using BASE::m_coupledWellsFlag;
  using BASE::m_numDofPerCell;
  using BASE::m_derivativeFluxResidual_dAperture;
  using BASE::m_fluxEstimate;
  using BASE::m_elemGhostRank;
  using BASE::m_volume;
  using BASE::m_gravCoef;
  using BASE::m_porosityRef;
  using BASE::m_elementArea;
  using BASE::m_elementAperture0;
  using BASE::m_elementAperture;
  using BASE::m_effectiveAperture;


  // Aliasing public/protected members/methods of SinglePhaseBase so we don't
  // have to use this->member etc.
  using BASE::m_pressure;
  using BASE::m_deltaPressure;
  using BASE::m_deltaVolume;
  using BASE::m_porosity;
  using BASE::m_mobility;
  using BASE::m_dMobility_dPres;
  using BASE::m_porosityOld;
  using BASE::m_densityOld;
  using BASE::m_pvMult;
  using BASE::m_dPvMult_dPres;
  using BASE::m_density;
  using BASE::m_dDens_dPres;
  using BASE::m_viscosity;
  using BASE::m_dVisc_dPres;
  using BASE::m_totalMeanStressOld;
  using BASE::m_totalMeanStress;
  using BASE::m_bulkModulus;
  using BASE::m_biotCoefficient;
  using BASE::m_poroMultiplier;
  using BASE::m_transTMultiplier;

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseFVM( const std::string & name,
                  dataRepository::Group * const parent );


  /// deleted default constructor
  SinglePhaseFVM() = delete;

  /// deleted copy constructor
  SinglePhaseFVM( SinglePhaseFVM const & ) = delete;

  /// default move constructor
  SinglePhaseFVM( SinglePhaseFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseFVM & operator=( SinglePhaseFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseFVM & operator=( SinglePhaseFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  template< typename _BASE=BASE >
  static
  typename std::enable_if< std::is_same< _BASE, SinglePhaseBase >::value, string >::type
  CatalogName()
  {
    return "SinglePhaseFVM";
  }

  template< typename _BASE=BASE >
  static
  typename std::enable_if< std::is_same< _BASE, SinglePhaseProppantBase >::value, string >::type
  CatalogName()
  {
    return "SinglePhaseProppantFVM";
  }

  /**
   * @defgroup Solver Interface Function
   *
   * This function provides the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  SetupDofs( DomainPartition const * const domain,
             DofManager & dofManager ) const override;

  virtual void
  SetupSystem( DomainPartition * const domain,
               DofManager & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition * const domain,
                           DofManager const & dofManager,
                           ParallelMatrix & matrix,
                           ParallelVector & rhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const * const domain,
                         DofManager const & dofManager,
                         ParallelVector const & rhs ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       ParallelVector const & solution,
                       real64 const scalingFactor,
                       DomainPartition * const domain ) override;

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleFluxTerms( real64 const time_n,
                     real64 const dt,
                     DomainPartition const * const domain,
                     DofManager const * const dofManager,
                     ParallelMatrix * const matrix,
                     ParallelVector * const rhs ) override;

  /**@}*/

  struct viewKeyStruct : SinglePhaseBase::viewKeyStruct
  {} viewKeysSinglePhaseFVM;

  viewKeyStruct & viewKeys()
  { return viewKeysSinglePhaseFVM; }

  viewKeyStruct const & viewKeys() const
  { return viewKeysSinglePhaseFVM; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysSinglePhaseFVM;

  groupKeyStruct & groupKeys()
  { return groupKeysSinglePhaseFVM; }

  groupKeyStruct const & groupKeys() const
  { return groupKeysSinglePhaseFVM; }

private:

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param time_n current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void ApplyFaceDirichletBC_implicit( real64 const time_n,
                                      real64 const dt,
                                      DofManager const * const dofManager,
                                      DomainPartition * const domain,
                                      ParallelMatrix * const matrix,
                                      ParallelVector * const rhs );

  // no data needed here, see SinglePhaseBase

};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_
