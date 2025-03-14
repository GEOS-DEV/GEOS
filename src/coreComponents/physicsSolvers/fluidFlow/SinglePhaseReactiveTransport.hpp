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
 * @file SinglePhaseReactiveTransport.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVETRANSPORT_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVETRANSPORT_HPP_

#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseReactiveTransportFields.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/AccumulationKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/FluidUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/reactive/ResidualNormKernel.hpp"
#include "constitutive/fluid/singlefluid/reactive/ReactiveSingleFluid.hpp"


namespace geos
{

/**
 * @class SinglePhaseReactiveTransport
 *
 * A solver for single phase reactive transport
 */
class SinglePhaseReactiveTransport : public SinglePhaseBase
{

public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseReactiveTransport( const string & name,
                                dataRepository::Group * const parent );

  SinglePhaseReactiveTransport() = delete;

  /// deleted copy constructor
  SinglePhaseReactiveTransport( SinglePhaseReactiveTransport const & ) = delete;

  /// default move constructor
  SinglePhaseReactiveTransport( SinglePhaseReactiveTransport && ) = default;

  /// deleted assignment operator
  SinglePhaseReactiveTransport & operator=( SinglePhaseReactiveTransport const & ) = delete;

  /// deleted move operator
  SinglePhaseReactiveTransport & operator=( SinglePhaseReactiveTransport && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseReactiveTransport() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "SinglePhaseReactiveTransport";
  }

  /**
   * @copydoc PhysicsSolverBase::getCatalogName()
   */
  string getCatalogName() const override { return catalogName(); }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  assembleSystem( real64 const time_n,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( real64 const & time_n,
                         real64 const & dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       real64 const dt,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void saveConvergedState( ElementSubRegionBase & subRegion ) const override final;

  virtual void
  updateState ( DomainPartition & domain ) override final;
  
  void updateSpeciesAmount( ElementSubRegionBase & subRegion ) const;
  
  virtual void updateFluidModel( ObjectManagerBase & dataGroup ) const override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializeFluidState( MeshLevel & mesh, arrayView1d< string const > const & regionNames ) override;

  /**
   * @brief assembles the accumulation terms in total mass balance and primary species amount equation for all cells
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  void assembleAccumulationTermsInMassBalanceAndSpeciesAmountEqs( real64 const dt,
                                                                  DomainPartition & domain,
                                                                  DofManager const & dofManager,
                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                  arrayView1d< real64 > const & localRhs ) const;
  
  /**
   * @brief assembles the flux terms for all cells
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void
  applyPresSpeciesDirichletBC( real64 const time_n,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Utility function that encapsulates the call to FieldSpecificationBase::applyFieldValue in BC application
   * @param[in] time_n the time at the beginning of the step
   * @param[in] dt the time step
   * @param[in] mesh the mesh level object
   * @param[in] logMessage the log message issued by the solver if the bc is called
   * @param[in] fieldKey the key of the field specified in the xml file
   * @param[in] boundaryFieldKey the key of the boundary field
   */
  template< typename OBJECT_TYPE >
  void applyFieldValue( real64 const & time_n,
                        real64 const & dt,
                        MeshLevel & mesh,
                        char const logMessage[],
                        string const fieldKey,
                        string const boundaryFieldKey ) const;

  virtual void
  applyAquiferBC( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) const override;
  
  virtual void
  assembleEDFMFluxTerms( real64 const time_n,
                         real64 const dt,
                         DomainPartition const & domain,
                         DofManager const & dofManager,
                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         arrayView1d< real64 > const & localRhs,
                         string const & jumpDofKey ) override final;
  
  virtual void
  assembleStabilizedFluxTerms( real64 const dt,
                               DomainPartition const & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

protected:

  virtual FluidPropViews getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const override;

  /// the number of primary species in the fluid
  integer m_numPrimarySpecies;

  /// name of the reactive fluid constitutive model
  string m_reactiveFluidModelName;
};

template< typename OBJECT_TYPE >
void SinglePhaseReactiveTransport::applyFieldValue( real64 const & time_n,
                                                    real64 const & dt,
                                                    MeshLevel & mesh,
                                                    char const logMessage[],
                                                    string const fieldKey,
                                                    string const boundaryFieldKey ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< OBJECT_TYPE >( time_n + dt,
                                  mesh,
                                  fieldKey,
                                  [&]( FieldSpecificationBase const & fs,
                                       string const & setName,
                                       SortedArrayView< localIndex const > const & lset,
                                       OBJECT_TYPE & targetGroup,
                                       string const & )
  {
    if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( lset.size() );
      GEOS_LOG_RANK_0( GEOS_FMT( logMessage,
                                 getName(), time_n+dt, fs.getCatalogName(), fs.getName(),
                                 setName, targetGroup.getName(), fs.getScale(), numTargetElems ) );
    }

    // Specify the bc value of the field
    fs.applyFieldValue< FieldSpecificationEqual,
                        parallelDevicePolicy<> >( lset,
                                                  time_n + dt,
                                                  targetGroup,
                                                  boundaryFieldKey );
  } );
}

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEREACTIVETRANSPORT_HPP_
