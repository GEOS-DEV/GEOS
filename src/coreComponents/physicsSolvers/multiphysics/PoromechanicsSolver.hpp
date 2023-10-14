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
  {
    this->registerWrapper( viewKeyStruct::useNAString(), &m_useNA ).
      setApplyDefaultValue( 0 ).
      setInputFlag( dataRepository::InputFlags::OPTIONAL ).
      setDescription( "Flag to indicate that the solver is going to use nonlinear acceleration" );
  }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    /// Flag to indicate that the solver is going to use nonlinear acceleration
    constexpr static char const * useNAString() { return "useNonlinearAcceleration"; }
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
                                                        array1d< real64 > & s )
  {
    // s denotes averageMeanTotalStressIncrement
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
          porosityModel.updateAverageMeanTotalStressIncrement( k, s[i] );
          i++;
        }
      } );
    } );
  }

  array1d< real64 > addTwoVecs( array1d< real64 > const & vec1,
                                array1d< real64 > const & vec2,
                                real64 const sign )
  {
    GEOS_ASSERT( vec1.size() == vec2.size());
    array1d< real64 > result;
    const localIndex N = vec1.size();
    for( localIndex i = 0; i < N; i++ )
    {
      result.emplace_back( vec1[i] + sign * vec2[i] );
    }
    return result;
  }

  array1d< real64 > scalarMultiplyAVec( array1d< real64 > const & vec,
                                        real64 const scalarMult )
  {
    array1d< real64 > result;
    const localIndex N = vec.size();
    for( localIndex i = 0; i < N; i++ )
    {
      result.emplace_back( scalarMult * vec[i] );
    }
    return result;
  }

  real64 dotTwoVecs( array1d< real64 > const & vec1,
                     array1d< real64 > const & vec2 )
  {
    GEOS_ASSERT( vec1.size() == vec2.size());
    real64 result = 0;
    const localIndex N = vec1.size();
    for( localIndex i = 0; i < N; i++ )
    {
      result += vec1[i] * vec2[i];
    }
    return result;
  }

  real64 computeAitkenRelaxationFactor( array1d< real64 > const & s0,
                                        array1d< real64 > const & s1,
                                        array1d< real64 > const & s1_tilde,
                                        array1d< real64 > const & s2_tilde,
                                        real64 const omega0 )
  {
    array1d< real64 > r1 = addTwoVecs( s1_tilde, s0, -1.0 );
    array1d< real64 > r2 = addTwoVecs( s2_tilde, s1, -1.0 );

    // diff = r2 - r1
    array1d< real64 > diff = addTwoVecs( r2, r1, -1.0 );

    const real64 denom = dotTwoVecs( diff, diff );
    const real64 numer = dotTwoVecs( r1, diff );

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
    return addTwoVecs( scalarMultiplyAVec( s1, 1.0 - omega1 ),
                       scalarMultiplyAVec( s2_tilde, omega1 ),
                       1.0 );
  }

  void startSequentialIteration( integer const & iter,
                                 DomainPartition & domain ) override
  {
    if( m_useNA )
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
    if( m_useNA )
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
    if( m_useNA && solverType == static_cast< integer >( SolverType::SolidMechanics ) )
    {
      recordAverageMeanTotalStressIncrement( domain, m_s2_tilde );
    }
  }

  virtual void
  postProcessInput() override
  {
    Base::postProcessInput();

    if( m_useNA && MpiWrapper::commSize( MPI_COMM_GEOSX ) > 1 )
    {
      GEOS_ERROR( "Nonlinear acceleration is not implemented for MPI runs" );
    }
  }

  /// member variables needed for Nonlinear Acceleration (Aitken)
  integer m_useNA;
  array1d< real64 > m_s0;
  array1d< real64 > m_s1;
  array1d< real64 > m_s1_tilde;
  array1d< real64 > m_s2;
  array1d< real64 > m_s2_tilde;
  real64 m_omega0;
  real64 m_omega1;

};

} /* namespace geos */

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSSOLVER_HPP_
