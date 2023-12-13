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
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

template< typename FLOW_SOLVER >
class PoromechanicsSolver : public CoupledSolver< FLOW_SOLVER,
                                                  SolidMechanicsLagrangianFEM >
{
public:

  using Base = CoupledSolver< FLOW_SOLVER,
                              SolidMechanicsLagrangianFEM >;
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
    : Base( name, parent )
  {}

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
        geos::constitutive::CoupledSolidBase & solid = SolverBase::getConstitutiveModel< geos::constitutive::CoupledSolidBase >(
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
        geos::constitutive::CoupledSolidBase & solid = SolverBase::getConstitutiveModel< geos::constitutive::CoupledSolidBase >(
          subRegion, solidName );
        auto & porosityModel = dynamic_cast< geos::constitutive::BiotPorosity const & >( solid.getBasePorosityModel());
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
    if( solverType == static_cast< integer >( SolverType::SolidMechanics ) &&
        this->getNonlinearSolverParameters().m_nonlinearAccelerationType== NonlinearSolverParameters::NonlinearAccelerationType::Aitken )
    {
      recordAverageMeanTotalStressIncrement( domain, m_s2_tilde );
    }
  }

  virtual void validateNonlinearAcceleration() override
  {
    if( MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 )
    {
      GEOS_ERROR( "Nonlinear acceleration is not implemented for MPI runs" );
    }
  }

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
