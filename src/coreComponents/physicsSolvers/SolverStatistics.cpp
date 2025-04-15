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

/*
 * @file SolverStatistics.cpp
 */

#include "SolverStatistics.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "fileIO/Outputs/OutputBase.hpp"

namespace geos
{

using namespace dataRepository;

SolverStatistics::SolverStatistics( string const & name, Group * const parent )
  : Group( name, parent ),
  m_currentNumOuterLoopIterations( 0 ),
  m_currentNumNonlinearIterations( 0 ),
  m_currentNumLinearIterations( 0 ),
  m_currentNewtonIter( 0 ),
  m_outputDir( joinPath( OutputBase::getOutputDirectory(), "convergence" ))
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_numTimeSteps ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time steps" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_numTimeStepCuts ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time step cuts" );

  registerWrapper( viewKeyStruct::numSuccessfulOuterLoopIterationsString(), &m_numSuccessfulOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful outer loop iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(), &m_numSuccessfulNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful nonlinear iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(), &m_numSuccessfulLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful linear iterations" );


  registerWrapper( viewKeyStruct::numDiscardedOuterLoopIterationsString(), &m_numDiscardedOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded outer loop iterations" );

  registerWrapper( viewKeyStruct::numDiscardedNonlinearIterationsString(), &m_numDiscardedNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded nonlinear iterations" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_numDiscardedLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded linear iterations" );

  makeDirsForPath( m_outputDir );
}

void SolverStatistics::prepareResidualTableLayout( bool isThermal )
{
  using TableLayoutArgs = std::initializer_list< std::variant< string_view, TableLayout::Column > >;

  m_nonLinearNormsLayout = isThermal
    ? std::make_unique< TableLayout >( TableLayoutArgs{
      std::variant< string_view, TableLayout::Column >{"Time-steps"},
      std::variant< string_view, TableLayout::Column >{"Newton Iter"},
      std::variant< string_view, TableLayout::Column >{"Rmass"},
      std::variant< string_view, TableLayout::Column >{"RVol"},
      std::variant< string_view, TableLayout::Column >{"REnergy"},
      std::variant< string_view, TableLayout::Column >{"Residual norm"}
    } )
    : std::make_unique< TableLayout >( TableLayoutArgs{
      std::variant< string_view, TableLayout::Column >{"Time-steps"},
      std::variant< string_view, TableLayout::Column >{"Newton Iter"},
      std::variant< string_view, TableLayout::Column >{"Rmass"},
      std::variant< string_view, TableLayout::Column >{"RVol"},
      std::variant< string_view, TableLayout::Column >{"Residual norm"}
    } );
}

void SolverStatistics::initializeTimeStepStatistics()
{
  // the time step begins, we reset the individual-timestep counters
  m_currentNumOuterLoopIterations = 0;
  m_currentNumNonlinearIterations = 0;
  m_currentNumLinearIterations = 0;
}

void SolverStatistics::logNonlinearIteration( integer const numLinearIterations )
{
  // we have just performed a Newton iteration, so we increment the individual-timestep counters
  m_currentNumNonlinearIterations++;
  m_currentNumLinearIterations += numLinearIterations;
}

void SolverStatistics::logNonlinearIteration()
{
  // we have just performed an outer iteration, so we increment the individual-timestep counter (number of outer iteration)
  m_currentNumNonlinearIterations++;
}

void SolverStatistics::logOuterLoopIteration()
{
  // we have just performed an outer loop iteration, so we increment the individual-timestep counter for outer loop iterations
  m_currentNumOuterLoopIterations++;
}


void SolverStatistics::logTimeStepCut()
{
  // we have just cut the time step, so we increment the cumulative counters for discarded timesteps
  m_numDiscardedOuterLoopIterations += m_currentNumOuterLoopIterations;
  m_numDiscardedNonlinearIterations += m_currentNumNonlinearIterations;
  m_numDiscardedLinearIterations += m_currentNumLinearIterations;
  m_numTimeStepCuts++;

  // we are going to restart the timestep from the previous converged time step, so we have to re-initialize the statistics
  initializeTimeStepStatistics();
}

void SolverStatistics::logNewtonIter( integer currentNewtonIter )
{ m_currentNewtonIter = currentNewtonIter; }

void SolverStatistics::logResidualNorm( real64 residualNorm,
                                        real64 residualMass,
                                        real64 residualVol )
{
  m_residualNorm = residualNorm;
  m_residualMass = residualMass;
  m_residualVol = residualVol;
}

void SolverStatistics::logThermalResidualNorm( real64 residualEnergy )
{ m_residualEnergy = residualEnergy; }

void SolverStatistics::saveTimeStepStatistics()
{
  // the timestep has converged, so we increment the cumulative counters for successful timesteps
  m_numSuccessfulOuterLoopIterations += m_currentNumOuterLoopIterations;
  m_numSuccessfulNonlinearIterations += m_currentNumNonlinearIterations;
  m_numSuccessfulLinearIterations += m_currentNumLinearIterations;
  m_numTimeSteps++;
}



void SolverStatistics::registerResidualNormToTable()
{
  m_nonLinearNormsData.addRow( getNumTimeSteps(),
                               getNewtonIter(),
                               getResidualNorm(),
                               getResidualMass(),
                               getResidualVol());
}

void SolverStatistics::registerThermalResidualNormToTable()
{
  m_nonLinearNormsData.addRow( getNumTimeSteps(),
                               getNewtonIter(),
                               getResidualNorm(),
                               getResidualMass(),
                               getResidualVol(),
                               getResidualEnergy());
}

void SolverStatistics::registerStatsToTable()
{
  m_nonLinearData.addRow( getNumTimeSteps(),
                          getNumTimeStepCuts(),
                          getNumSuccessfulOuterLoopIterations(),
                          getNumSuccessfulNonlinearIterations(),
                          getNumSuccessfulLinearIterations(),
                          getNumDiscardedOuterLoopIterations(),
                          getNumDiscardedNonlinearIterations(),
                          getNumDiscardedLinearIterations());
}

void SolverStatistics::outputStatistics( bool writeCSV )
{
  { // output to log
    GEOS_LOG_RANK_0( GEOS_FMT( "{}, number of Time-steps: {}", getParent().getName(), m_numTimeSteps ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "{}, number of Time steps cut: {}", getParent().getName(), m_numTimeStepCuts ) );
    TableLayout statsLayout( GEOS_FMT( "{} statistics", getParent().getName()),
                             { TableLayout::Column()
                                 .setName( "Components" )
                                 .setValuesAlignment( TableLayout::Alignment::left ),
                               "Iter"} );
    TableTextFormatter const statsFormatter( statsLayout );

    TableData nonLinearDataLog;
    nonLinearDataLog.addRow( "Successful outer loop", m_numSuccessfulOuterLoopIterations );
    nonLinearDataLog.addRow( "Successful nonlinear", m_numSuccessfulNonlinearIterations );
    nonLinearDataLog.addRow( "Successful linear", m_numSuccessfulLinearIterations );
    nonLinearDataLog.addRow( "Discarded outer loop", m_numDiscardedOuterLoopIterations );
    nonLinearDataLog.addRow( "Discarded nonlinear", m_numDiscardedNonlinearIterations );
    nonLinearDataLog.addRow( "Discarded linear", m_numDiscardedLinearIterations );

    std::cout << statsFormatter.toString( nonLinearDataLog ) << std::endl;
  }

  if( writeCSV )
  {
    string const fileName = GEOS_FMT( "{}/statistics_cycle_{}.csv", m_outputDir, m_numTimeSteps );
    std::ofstream logStream( fileName );
    std::cout << "fileName : " << fileName <<  std::endl;

    TableLayout statsLayout( {  "Time-steps", "Time steps cut",
                                "Successful outer loop", "Successful nonlinear", "Successful linear",
                                "Discarded outer loop", "Discarded nonlinear", "Discarded linear"} );
    TableCSVFormatter const csvOutput( statsLayout );
    std::cout << csvOutput.toString( m_nonLinearData ) << std::endl;

    logStream << csvOutput.toString( m_nonLinearData );
    logStream.close();
    m_nonLinearNormsData.clear();
  }
}

void SolverStatistics::outputResidualNorm( bool writeCSV )
{
  std::ofstream logStream( m_residualNormsFileName );

  if( writeCSV )
  {
    TableCSVFormatter const csvOutput( *m_nonLinearNormsLayout );
    std::cout << csvOutput.toString( m_nonLinearNormsData ) << std::endl;
    logStream << csvOutput.toString( m_nonLinearNormsData );
  }
  logStream.close();
  m_nonLinearNormsData.clear();

}
} // namespace geos
