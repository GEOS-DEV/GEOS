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

#include "fileIO/Outputs/OutputBase.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

using namespace dataRepository;

SolverStatistics::SolverStatistics( string const & name, Group * const parent )
  : Group( name, parent ),
  m_iterationsStats( groupKeyStruct::IterationsStatisticsString(), this ),
  m_convergenceStats()
{
  m_outputDir =  joinPath( OutputBase::getOutputDirectory(), m_directoryName );
}

void SolverStatistics::setOutputFilesName( string_view solverName )
{
  m_iterationsStats.setFilename( GEOS_FMT( "{}/{}_iterations.csv", m_outputDir, solverName ));
  m_convergenceStats.setFilename( GEOS_FMT( "{}/{}_convergence.csv", m_outputDir, solverName ));
}

IterationsStatistics::IterationsStatistics( string const & name, Group * const parent )
  : Group( name, parent ),
  m_currentNumConfigIterations( 0 ),
  m_currentNumNonlinearIterations( 0 ),
  m_currentNumLinearIterations( 0 )
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_numTimeSteps ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time steps" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_numTimeStepCuts ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time step cuts" );

  registerWrapper( viewKeyStruct::numSuccessfulConfigIterationsString(),
                   &m_numSuccessfulConfigIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful configuration iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(),
                   &m_numSuccessfulNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful nonlinear iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(),
                   &m_numSuccessfulLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful linear iterations" );

  registerWrapper( viewKeyStruct::numDiscardedConfigIterationsString(),
                   &m_numDiscardedConfigIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded configuration iterations" );

  registerWrapper( viewKeyStruct::numDiscardedNonlinearIterationsString(),
                   &m_numDiscardedNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded nonlinear iterations" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(),
                   &m_numDiscardedLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded linear iterations" );

  m_iterationCSVLayout.setTitle( GEOS_FMT( "{} iterations", getParent().getName()));
  m_iterationCSVLayout.addColumns( { "Number of time steps", "Number of time step cuts",
                                     "Successful configuration", "Successful nonlinear", "Successful linear",
                                     "Discarded configuration", "Discarded nonlinear", "Discarded linear",
                                     "Setup time", "Solve time"} );
}

void IterationsStatistics::resetCurrentTimeStepStatistics()
{
  m_currentNumConfigIterations = 0;
  m_currentNumNonlinearIterations = 0;
  m_currentNumLinearIterations = 0;
}

void IterationsStatistics::incrementConfigIteration()
{
  m_currentNumConfigIterations++;
}

void IterationsStatistics::updateNonlinearIteration( integer const numLinearIterations )
{
  m_currentNumNonlinearIterations++;
  m_currentNumLinearIterations += numLinearIterations;
}

void IterationsStatistics::accumulateSolverLinearTime( real64 setupTime, real64 solveTime )
{
  m_setupTime += setupTime;
  m_solveTime += solveTime;
}

void IterationsStatistics::resetSolverLinearTime()
{
  m_setupTime = 0.0;
  m_solveTime = 0.0;
}

void IterationsStatistics::iterateTimeStepStatistics()
{
  m_numSuccessfulConfigIterations += m_currentNumConfigIterations;
  m_numSuccessfulNonlinearIterations += m_currentNumNonlinearIterations;
  m_numSuccessfulLinearIterations += m_currentNumLinearIterations;
  m_numTimeSteps++;
}

void IterationsStatistics::updateTimeStepCut()
{
  m_numDiscardedConfigIterations += m_currentNumConfigIterations;
  m_numDiscardedNonlinearIterations += m_currentNumNonlinearIterations;
  m_numDiscardedLinearIterations += m_currentNumLinearIterations;
  m_numTimeStepCuts++;

  // we are going to restart the timestep from the previous converged time step, so we have to re-initialize the statistics
  resetCurrentTimeStepStatistics();
}

void IterationsStatistics::writeIterationStatsToTable()
{
  if( m_numTimeSteps == 0 || !m_logOutputRequest || MpiWrapper::commRank() != 0 )
    return;

  m_iterationData.addRow( m_numTimeSteps,
                          m_numTimeStepCuts,
                          m_numSuccessfulConfigIterations,
                          m_numSuccessfulNonlinearIterations,
                          m_numSuccessfulLinearIterations,
                          m_numDiscardedConfigIterations,
                          m_numDiscardedNonlinearIterations,
                          m_numDiscardedLinearIterations,
                          m_setupTime,
                          m_solveTime );

  if( !m_CSVOutputOpened )
  {
    m_logStream.open( m_iterationsFilename );
    m_iterationCSVFormatter = std::make_unique< TableCSVFormatter >( m_iterationCSVLayout );
    m_logStream << m_iterationCSVFormatter->headerToString( );
    m_CSVOutputOpened  = true;
  }

  m_logStream << m_iterationCSVFormatter->dataToString( m_iterationData );
  m_logStream.flush();
  m_iterationData.clear();

  resetSolverLinearTime();
}

void IterationsStatistics::outputStatistics() const
{
  // no statistics to output when no time steps have been recorded
  if( m_numTimeSteps == 0 || !m_logOutputRequest )
    return;

  TableLayout iterationLogLayout ( GEOS_FMT( "{}", m_tableIterationName ), {"Statistics", "Value"} );

  TableTextFormatter const statsFormatter( iterationLogLayout );

  TableData iterationDataLog;
  iterationDataLog.addRow( "Time steps", m_numTimeSteps );
  iterationDataLog.addRow( "Time step cuts", m_numTimeStepCuts );
  if( m_numSuccessfulConfigIterations > 0 ) // only print when solver did any configuration iterations
    iterationDataLog.addRow( "Successful configuration iterations", m_numSuccessfulConfigIterations );
  iterationDataLog.addRow( "Successful nonlinear iterations", m_numSuccessfulNonlinearIterations );
  if( m_numSuccessfulLinearIterations > 0 ) // only print when solver did any linear iterations
    iterationDataLog.addRow( "Successful linear iterations", m_numSuccessfulLinearIterations );
  if( m_numSuccessfulConfigIterations > 0 ) // only print when solver did any configuration iterations
    iterationDataLog.addRow( "Discarded configuration iterations", m_numDiscardedConfigIterations );
  iterationDataLog.addRow( "Discarded nonlinear iterations", m_numDiscardedNonlinearIterations );
  if( m_numSuccessfulLinearIterations > 0 ) // only print when solver did any linear iterations
    iterationDataLog.addRow( "Discarded linear iterations", m_numDiscardedLinearIterations );

  GEOS_LOG_RANK_0( statsFormatter.toString( iterationDataLog ));
}

ConvergenceStatistics::ConvergenceStatistics():
  m_convergenceLayout()
{}

void ConvergenceStatistics::writeConvergenceStatsToTable()
{
  if( !m_CSVOutputRequest || MpiWrapper::commRank() != 0 )
    return;
  stdVector< TableData::CellData > residualsNormCells;

  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_cycleNumber )} ));
  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_time_n )} ));
  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_dt )} ));
  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_iteration )} ));

  for( auto const & residual : m_residuals )
  {
    residualsNormCells.emplace_back( TableData::CellData(
      {
        CellType::Value,
        !std::isnan( residual.second ) ? GEOS_FMT( "{}", residual.second ) : "0"
      } ));
  }

  m_convergenceData.addRow( residualsNormCells );

  if( !m_CSVOutputOpened )
  {
    string_array header = {"Cycle number", "time_n (s)", "dt (s)", "iteration"};
    for( auto const & residual : m_residuals )
    {
      header.emplace_back( residual.first );
    }
    m_convergenceLayout.addColumns( header );

    m_logStream.open( m_convergenceFilename );
    m_convergenceFormatter = std::make_unique< TableCSVFormatter >( m_convergenceLayout );
    m_logStream << m_convergenceFormatter->headerToString( );
    m_CSVOutputOpened  = true;
  }

  m_logStream << m_convergenceFormatter->dataToString( m_convergenceData );
  m_logStream.flush();
  m_convergenceData.clear();
}

void ConvergenceStatistics::updateSolverStep( real64 const & time_n, real64 const & dt,
                                              integer const cycleNumber, integer const newtonIter )
{
  m_time_n = time_n;
  m_dt = dt;
  m_cycleNumber = cycleNumber;
  m_iteration = newtonIter;
}

} // namespace geos
