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
 * @file PoromechanicsSolver.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSSOLVER_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSSOLVER_HPP_

#include "physicsSolvers/multiphysics/CoupledSolver.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

template< typename FLOW_SOLVER, typename MECHANICS_SOLVER = SolidMechanicsLagrangianFEM >
class PoromechanicsSolver : public CoupledSolver< FLOW_SOLVER, MECHANICS_SOLVER >
{
public:

  using Base = CoupledSolver< FLOW_SOLVER, MECHANICS_SOLVER >;
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
   * @brief main constructor for CoupledSolver Objects
   * @param name the name of this instantiation of CoupledSolver in the repository
   * @param parent the parent group of this instantiation of CoupledSolver
   */
  PoromechanicsSolver( const string & name,
                       dataRepository::Group * const parent )
    : Base( name, parent ),
    m_isThermal( 0 )
  {
    this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
      setApplyDefaultValue( 0 ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "Flag indicating whether the problem is thermal or not. Set isThermal=\"1\" to enable the thermal coupling" );

    this->registerWrapper( viewKeyStruct::performStressInitializationString(), &m_performStressInitialization ).
      setApplyDefaultValue( false ).
      setInputFlag( dataRepository::InputFlags::FALSE ).
      setDescription( "Flag to indicate that the solver is going to perform stress initialization" );
  }

  virtual void initializePostInitialConditionsPreSubGroups() override
  {
    Base::initializePostInitialConditionsPreSubGroups();

    GEOS_THROW_IF( this->m_isThermal && !this->flowSolver()->isThermal(),
                   GEOS_FMT( "{} {}: The attribute `{}` of the flow solver `{}` must be set to 1 since the poromechanics solver is thermal",
                             this->getCatalogName(), this->getName(), FlowSolverBase::viewKeyStruct::isThermalString(), this->flowSolver()->getName() ),
                   InputError );
  }

  virtual void initializePreSubGroups() override
  {
    Base::initializePreSubGroups();

    DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );

    this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                 MeshLevel & mesh,
                                                                                 arrayView1d< string const > const & regionNames )
    {
      ElementRegionManager & elementRegionManager = mesh.getElemManager();
      elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                         [&]( localIndex const,
                                                                              ElementSubRegionBase & subRegion )
      {
        string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
        porousName = this->template getConstitutiveName< constitutive::CoupledSolidBase >( subRegion );
        GEOS_THROW_IF( porousName.empty(),
                       GEOS_FMT( "{} {} : Solid model not found on subregion {}",
                                 this->getCatalogName(), this->getDataContext().toString(), subRegion.getName() ),
                       InputError );

        string & porosityModelName = subRegion.getReference< string >( constitutive::CoupledSolidBase::viewKeyStruct::porosityModelNameString() );
        porosityModelName = this->template getConstitutiveName< constitutive::PorosityBase >( subRegion );
        GEOS_THROW_IF( porosityModelName.empty(),
                       GEOS_FMT( "{} {} : Porosity model not found on subregion {}",
                                 this->catalogName(), this->getDataContext().toString(), subRegion.getName() ),
                       InputError );

        if( subRegion.hasField< fields::poromechanics::bulkDensity >() )
        {
          // get the solid model to know the number of quadrature points and resize the bulk density
          constitutive::CoupledSolidBase const & solid = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, porousName );
          subRegion.getField< fields::poromechanics::bulkDensity >().resizeDimension< 1 >( solid.getDensity().size( 1 ) );
        }
      } );
    } );
  }

  virtual void registerDataOnMesh( dataRepository::Group & meshBodies ) override
  {
    SolverBase::registerDataOnMesh( meshBodies );

    if( this->getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
    {
      // to let the solid mechanics solver that there is a pressure and temperature RHS in the mechanics solve
      solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate();
      // to let the flow solver that saving pressure_k and temperature_k is necessary (for the fixed-stress porosity terms)
      flowSolver()->enableFixedStressPoromechanicsUpdate();
    }

    SolverBase::forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                [&]( localIndex const,
                                                                     ElementSubRegionBase & subRegion )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::porousMaterialNamesString() ).
          setPlotLevel( dataRepository::PlotLevel::NOPLOT ).
          setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 );

        // This is needed by the way the surface generator currently does things.
        subRegion.registerWrapper< string >( constitutive::CoupledSolidBase::viewKeyStruct::porosityModelNameString() ).
          setPlotLevel( dataRepository::PlotLevel::NOPLOT ).
          setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).
          setSizedFromParent( 0 );

        if( this->getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
        {
          // register the bulk density for use in the solid mechanics solver
          // ideally we would resize it here as well, but the solid model name is not available yet (see below)
          subRegion.registerField< fields::poromechanics::bulkDensity >( this->getName() );
        }
      } );
    } );
  }

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override
  {
    flowSolver()->keepFlowVariablesConstantDuringInitStep( m_performStressInitialization );
    Base::implicitStepSetup( time_n, dt, domain );
  }

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override
  {
    // note that the order of operations matters a lot here (for instance for the MGR labels)
    // we must set up dofs for solid mechanics first, and then for flow
    // that's the reason why this function is here and not in CoupledSolvers.hpp
    solidMechanicsSolver()->setupDofs( domain, dofManager );
    flowSolver()->setupDofs( domain, dofManager );

    this->setupCoupling( domain, dofManager );
  }

  virtual bool checkSequentialConvergence( int const & iter,
                                           real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition & domain ) override
  {
    // always force outer loop for initialization
    auto & subcycling = this->getNonlinearSolverParameters().m_subcyclingOption;
    auto const subcycling_orig = subcycling;
    if( m_performStressInitialization )
      subcycling = 1;

    bool isConverged = Base::checkSequentialConvergence( iter, time_n, dt, domain );

    // restore original
    subcycling = subcycling_orig;

    return isConverged;
  }

  /**
   * @brief accessor for the pointer to the solid mechanics solver
   * @return a pointer to the solid mechanics solver
   */
  MECHANICS_SOLVER * solidMechanicsSolver() const
  {
    return std::get< toUnderlying( SolverType::SolidMechanics ) >( m_solvers );
  }

  /**
   * @brief accessor for the pointer to the flow solver
   * @return a pointer to the flow solver
   */
  FLOW_SOLVER * flowSolver() const
  {
    return std::get< toUnderlying( SolverType::Flow ) >( m_solvers );
  }

  /*
   * @brief Utility function to set the stress initialization flag
   * @param[in] performStressInitialization true if the solver has to initialize stress, false otherwise
   */
  void setStressInitialization( integer const performStressInitialization )
  {
    m_performStressInitialization = performStressInitialization;
  }

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

  /* Implementation of Nonlinear Acceleration (Aitken) of averageMeanTotalStressIncrement */

  void recordAverageMeanTotalStressIncrement( DomainPartition & domain,
                                              array1d< real64 > & averageMeanTotalStressIncrement )
  {
    averageMeanTotalStressIncrement.resize( 0 );
    SolverBase::forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                             MeshLevel & mesh,
                                                                             arrayView1d< string const > const & regionNames ) {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion ) {
        // get the solid model (to access stress increment)
        string const solidName = subRegion.template getReference< string >( "porousMaterialNames" );
        constitutive::CoupledSolidBase & solid = SolverBase::getConstitutiveModel< constitutive::CoupledSolidBase >(
          subRegion, solidName );

        arrayView1d< const real64 > const & averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();
        for( localIndex k = 0; k < localIndex( averageMeanTotalStressIncrement_k.size()); k++ )
        {
          averageMeanTotalStressIncrement.emplace_back( averageMeanTotalStressIncrement_k[k] );
        }
      } );
    } );
  }

  void applyAcceleratedAverageMeanTotalStressIncrement( DomainPartition & domain,
                                                        array1d< real64 > & averageMeanTotalStressIncrement )
  {
    integer i = 0;
    SolverBase::forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                             MeshLevel & mesh,
                                                                             arrayView1d< string const > const & regionNames ) {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion ) {
        // get the solid model (to access stress increment)
        string const solidName = subRegion.template getReference< string >( "porousMaterialNames" );
        constitutive::CoupledSolidBase & solid = SolverBase::getConstitutiveModel< constitutive::CoupledSolidBase >(
          subRegion, solidName );
        auto & porosityModel = dynamic_cast< constitutive::BiotPorosity const & >( solid.getBasePorosityModel());
        arrayView1d< real64 > const & averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();
        for( localIndex k = 0; k < localIndex( averageMeanTotalStressIncrement_k.size()); k++ )
        {
          porosityModel.updateAverageMeanTotalStressIncrement( k, averageMeanTotalStressIncrement[i] );
          i++;
        }
      } );
    } );
  }

  real64 computeAitkenRelaxationFactor( array1d< real64 > const & s0,
                                        array1d< real64 > const & s1,
                                        array1d< real64 > const & s1_tilde,
                                        array1d< real64 > const & s2_tilde,
                                        real64 const omega0 )
  {
    array1d< real64 > r1 = axpy( s1_tilde, s0, -1.0 );
    array1d< real64 > r2 = axpy( s2_tilde, s1, -1.0 );

    // diff = r2 - r1
    array1d< real64 > diff = axpy( r2, r1, -1.0 );

    real64 const denom = dot( diff, diff );
    real64 const numer = dot( r1, diff );

    real64 omega1 = 1.0;
    if( !isZero( denom ))
    {
      omega1 = -1.0 * omega0 * numer / denom;
    }
    return omega1;
  }

  array1d< real64 > computeUpdate( array1d< real64 > const & s1,
                                   array1d< real64 > const & s2_tilde,
                                   real64 const omega1 )
  {
    return axpy( scale( s1, 1.0 - omega1 ),
                 scale( s2_tilde, omega1 ),
                 1.0 );
  }

  void startSequentialIteration( integer const & iter,
                                 DomainPartition & domain ) override
  {
    if( this->getNonlinearSolverParameters().m_nonlinearAccelerationType == NonlinearSolverParameters::NonlinearAccelerationType::Aitken )
    {
      if( iter == 0 )
      {
        recordAverageMeanTotalStressIncrement( domain, m_s1 );
      }
      else
      {
        m_s0 = m_s1;
        m_s1 = m_s2;
        m_s1_tilde = m_s2_tilde;
        m_omega0 = m_omega1;
      }
    }
  }

  void finishSequentialIteration( integer const & iter,
                                  DomainPartition & domain ) override
  {
    if( this->getNonlinearSolverParameters().m_nonlinearAccelerationType == NonlinearSolverParameters::NonlinearAccelerationType::Aitken )
    {
      if( iter == 0 )
      {
        m_s2 = m_s2_tilde;
        m_omega1 = 1.0;
      }
      else
      {
        m_omega1 = computeAitkenRelaxationFactor( m_s0, m_s1, m_s1_tilde, m_s2_tilde, m_omega0 );
        m_s2 = computeUpdate( m_s1, m_s2_tilde, m_omega1 );
        applyAcceleratedAverageMeanTotalStressIncrement( domain, m_s2 );
      }
    }
  }

  virtual void mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType ) override
  {
    GEOS_MARK_FUNCTION;

    /// After the flow solver
    if( solverType == static_cast< integer >( SolverType::Flow ) )
    {
      this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                                  MeshLevel & mesh,
                                                                                  arrayView1d< string const > const & regionNames )
      {

        mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                              auto & subRegion )
        {
          // update the bulk density
          // TODO: ideally, we would not recompute the bulk density, but a more general "rhs" containing the body force and the
          // pressure/temperature terms
          updateBulkDensity( subRegion );
        } );
      } );
    }

    /// After the solid mechanics solver
    if( solverType == static_cast< integer >( SolverType::SolidMechanics ) )
    {
      // compute the average of the mean total stress increment over quadrature points
      averageMeanTotalStressIncrement( domain );

      this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                                  MeshLevel & mesh,
                                                                                  arrayView1d< string const > const & regionNames )
      {

        mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                              auto & subRegion )
        {
          // update the porosity after a change in displacement (after mechanics solve)
          // or a change in pressure/temperature (after a flow solve)
          flowSolver()->updatePorosityAndPermeability( subRegion );
          // update bulk density to reflect porosity change into mechanics
          updateBulkDensity( subRegion );
        } );
      } );
    }

    // needed to perform nonlinear acceleration
    if( solverType == static_cast< integer >( SolverType::SolidMechanics ) &&
        this->getNonlinearSolverParameters().m_nonlinearAccelerationType== NonlinearSolverParameters::NonlinearAccelerationType::Aitken )
    {
      recordAverageMeanTotalStressIncrement( domain, m_s2_tilde );
    }
  }

  /**
   * @brief Helper function to average the mean total stress increment over quadrature points
   * @param[in] domain the domain partition
   */
  void averageMeanTotalStressIncrement( DomainPartition & domain )
  {
    this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                                MeshLevel & mesh,
                                                                                arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        // get the solid model (to access stress increment)
        string const solidName = subRegion.template getReference< string >( viewKeyStruct::porousMaterialNamesString() );
        constitutive::CoupledSolidBase & solid = this->template getConstitutiveModel< constitutive::CoupledSolidBase >( subRegion, solidName );

        arrayView2d< real64 const > const meanTotalStressIncrement_k = solid.getMeanTotalStressIncrement_k();
        arrayView1d< real64 > const averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();

        finiteElement::FiniteElementBase & subRegionFE =
          subRegion.template getReference< finiteElement::FiniteElementBase >( solidMechanicsSolver()->getDiscretizationName() );

        // determine the finite element type
        finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
        dispatch3D( subRegionFE, [&] ( auto const finiteElement )
        {
          using FE_TYPE = decltype( finiteElement );

          // call the factory and launch the kernel
          AverageOverQuadraturePoints1DKernelFactory::
            createAndLaunch< CellElementSubRegion,
                             FE_TYPE,
                             parallelDevicePolicy<> >( mesh.getNodeManager(),
                                                       mesh.getEdgeManager(),
                                                       mesh.getFaceManager(),
                                                       subRegion,
                                                       finiteElement,
                                                       meanTotalStressIncrement_k,
                                                       averageMeanTotalStressIncrement_k );
        } );
      } );
    } );
  }

  virtual void updateBulkDensity( ElementSubRegionBase & subRegion ) = 0;

  virtual void validateNonlinearAcceleration() override
  {
    if( MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 )
    {
      GEOS_ERROR( "Nonlinear acceleration is not implemented for MPI runs" );
    }
  }

  /// Flag to determine whether or not this is a thermal simulation
  integer m_isThermal;

  /// Flag to indicate that the solver is going to perform stress initialization
  integer m_performStressInitialization;

  /// Member variables needed for Nonlinear Acceleration ( Aitken ). Naming convention follows ( Jiang & Tchelepi, 2019 )
  array1d< real64 > m_s0; // Accelerated averageMeanTotalStresIncrement @ outer iteration v ( two iterations ago )
  array1d< real64 > m_s1; // Accelerated averageMeanTotalStresIncrement @ outer iteration v + 1 ( previous iteration )
  array1d< real64 > m_s1_tilde; // Unaccelerated averageMeanTotalStresIncrement @ outer iteration v + 1 ( previous iteration )
  array1d< real64 > m_s2; // Accelerated averageMeanTotalStresIncrement @ outer iteration v + 2 ( current iteration )
  array1d< real64 > m_s2_tilde; // Unaccelerated averageMeanTotalStresIncrement @ outer iteration v + 1 ( current iteration )
  real64 m_omega0; // Old Aitken relaxation factor
  real64 m_omega1; // New Aitken relaxation factor

};

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSSOLVER_HPP_
