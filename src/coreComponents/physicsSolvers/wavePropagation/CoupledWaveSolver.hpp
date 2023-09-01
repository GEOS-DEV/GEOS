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
 * @file CoupledWaveSolver.hpp
 *
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_COUPLEDWAVESOLVER_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_COUPLEDWAVESOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

#include <tuple>

namespace geos
{

template< typename ... SOLVERS >
class CoupledWaveSolver : public SolverBase
{

public:

  /**
   * @brief main constructor for CoupledWaveSolver Objects
   * @param name the name of this instantiation of CoupledWaveSolver in the repository
   * @param parent the parent group of this instantiation of CoupledWaveSolver
   */
  CoupledWaveSolver( const string & name,
                     Group * const parent )
    : SolverBase( name, parent )
  {
    forEachArgInTuple( m_solvers, [&]( auto solver, auto idx )
    {
      using SolverType = TYPEOFPTR( solver );
      string const key = SolverType::coupledSolverAttributePrefix() + "SolverName";
      registerWrapper( key, &m_names[idx()] ).
        setInputFlag( dataRepository::InputFlags::REQUIRED ).
        setDescription( "Name of the " + SolverType::coupledSolverAttributePrefix() + " solver used by the coupled solver" );
    } );
  }

  /// deleted copy constructor
  CoupledWaveSolver( CoupledWaveSolver const & ) = delete;

  /// default move constructor
  CoupledWaveSolver( CoupledWaveSolver && ) = default;

  /// deleted assignment operator
  CoupledWaveSolver & operator=( CoupledWaveSolver const & ) = delete;

  /// deleted move operator
  CoupledWaveSolver & operator=( CoupledWaveSolver && ) = delete;

  virtual void
  postProcessInput() override
  {
    std::cout << "\t[CoupledWaveSolver::postProcessInput]" << std::endl;
    SolverBase::postProcessInput();

    forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
    {
      using SolverPtr = TYPEOFREF( solver );
      using SolverType = TYPEOFPTR( SolverPtr {} );
      solver = this->getParent().template getGroupPointer< SolverType >( m_names[idx()] );
      GEOS_THROW_IF( solver == nullptr,
                     GEOS_FMT( "Could not find solver '{}' of type {}",
                               m_names[idx()], LvArray::system::demangleType< SolverType >() ),
                     InputError );
      // solver->postProcessInput();  // error: function "geos::AcousticWaveEquationSEM::postProcessInput" (declared at line 160 of
      // AcousticWaveEquationSEM.hpp) is inaccessible
    } );
  }

  /*
     virtual void
     cleanup( real64 const time_n,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           DomainPartition & domain ) override
     {
     std::cout << "\t[CoupledWaveSolver::cleanup]" << std::endl;

     // call the base class cleanup (for reporting purposes)
     SolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

     forEachArgInTuple( m_solvers, [&]( auto & solver, auto idx )
     {
      using SolverPtr = TYPEOFREF( solver );
      using SolverType = TYPEOFPTR( SolverPtr {} );
      solver = this->getParent().template getGroupPointer< SolverType >( m_names[idx()] );
      solver->cleanup(time_n, cycleNumber, eventCounter, eventProgress, domain);
     } );
     }
   */

protected:

  /// Pointers of the single-physics solvers
  std::tuple< SOLVERS *... > m_solvers;

  /// Names of the single-physics solvers
  std::array< string, sizeof...( SOLVERS ) > m_names;
};

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_COUPLEDWAVESOLVER_HPP_ */
