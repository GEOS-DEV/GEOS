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
  m_currentNewtonIter( 0 ),
  m_currentNumOuterLoopIterations( 0 ),
  m_currentNumNonlinearIterations( 0 ),
  m_currentNumLinearIterations( 0 )
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_numTimeSteps ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time steps" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_numTimeStepCuts ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time step cuts" );

  registerWrapper( viewKeyStruct::numSuccessfulOuterLoopIterationsString(),
                   &m_numSuccessfulOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful outer loop iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(),
                   &m_numSuccessfulNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful nonlinear iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(),
                   &m_numSuccessfulLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful linear iterations" );


  registerWrapper( viewKeyStruct::numDiscardedOuterLoopIterationsString(),
                   &m_numDiscardedOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded outer loop iterations" );

  registerWrapper( viewKeyStruct::numDiscardedNonlinearIterationsString(),
                   &m_numDiscardedNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded nonlinear iterations" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(),
                   &m_numDiscardedLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded linear iterations" );

  m_iterationCSVLayout = std::make_unique< TableLayout >();
  m_iterationCSVLayout->setTitle( GEOS_FMT( "{} iterations", getParent().getName()));
  m_iterationCSVLayout->addColumns( { "numTimeSteps", "newton Iteration", "numTimeStepCuts",
                                      "setup time", "solve time",
                                      "Successful outer loop", "Successful nonlinear", "Successful linear",
                                      "Discarded outer loop", "Discarded nonlinear", "Discarded linear"} );
  //m_iterationCSVLayout->addColumn( "Iter" );

}

void IterationsStatistics::resetCurrentTimeStepStatistics()
{
  // the time step begins, we reset the individual-timestep counters
  if( m_isIterativeSolver )
  {
    m_currentNumOuterLoopIterations = 0;
    m_currentNumNonlinearIterations = 0;
    m_currentNumLinearIterations = 0;
  }
}

void IterationsStatistics::updateNonlinearIteration( integer const numLinearIterations )
{
  if( m_isIterativeSolver )
  {
    // we have just performed a Newton iteration, so we increment the individual-timestep counters
    m_currentNumNonlinearIterations++;
    m_currentNumLinearIterations += numLinearIterations;
  }
}

void IterationsStatistics::updateNonlinearIteration()
{
  if( m_isIterativeSolver )
  {
    // we have just performed an outer iteration, so we increment the individual-timestep counter (number of outer iteration)
    m_currentNumNonlinearIterations++;
  }
}

void IterationsStatistics::accumulateSolverLinearTime( real64 setupTime, real64 solveTime )
{
  if( m_isIterativeSolver )
  {
    m_setupTime += setupTime;
    m_solveTime += solveTime;
  }
}

void IterationsStatistics::resetSolverLinearTime()
{
  if( m_isIterativeSolver )
  {
    m_setupTime = 0.0;
    m_solveTime = 0.0;
  }
}

void IterationsStatistics::incrementNonlinearIteration()
{
  if( m_isIterativeSolver )
  {
    // we have just performed an outer loop iteration, so we increment the individual-timestep counter for outer loop iterations
    m_currentNumOuterLoopIterations++;
  }
}

void IterationsStatistics::iterateTimeStepStatistics()
{
  if( m_isIterativeSolver )
  {
    // the timestep has converged, so we increment the cumulative counters for successful timesteps
    m_numSuccessfulOuterLoopIterations += m_currentNumOuterLoopIterations;
    m_numSuccessfulNonlinearIterations += m_currentNumNonlinearIterations;
    m_numSuccessfulLinearIterations += m_currentNumLinearIterations;
    m_numTimeSteps++;
  }
}

void IterationsStatistics::updateTimeStepCut()
{
  if( m_isIterativeSolver )
  {
    // we have just cut the time step, so we increment the cumulative counters for discarded timesteps
    m_numDiscardedOuterLoopIterations += m_currentNumOuterLoopIterations;
    m_numDiscardedNonlinearIterations += m_currentNumNonlinearIterations;
    m_numDiscardedLinearIterations += m_currentNumLinearIterations;
    m_numTimeStepCuts++;

    // we are going to restart the timestep from the previous converged time step, so we have to re-initialize the statistics
    resetCurrentTimeStepStatistics();
  }
}

void IterationsStatistics::writeStatsToTable()
{
  if( !m_isIterativeSolver || !m_csvOutput )
    return;

  m_iterationData.addRow( m_numTimeSteps,
                          m_currentNewtonIter,
                          m_numTimeStepCuts,
                          m_setupTime,
                          m_solveTime,
                          m_numSuccessfulOuterLoopIterations,
                          m_numSuccessfulNonlinearIterations,
                          m_numSuccessfulLinearIterations,
                          m_numDiscardedOuterLoopIterations,
                          m_numDiscardedNonlinearIterations,
                          m_numDiscardedLinearIterations );

  if( !logStream.is_open() )
  {
    logStream.open( m_iterationsFilename );
    m_iterationCSVFormatter.reset( new TableCSVFormatter( *m_iterationCSVLayout ));
    logStream << m_iterationCSVFormatter->headerToString( );
  }

  logStream << m_iterationCSVFormatter->dataToString( m_iterationData );
  logStream.flush();
  m_iterationData.clear();
}

void IterationsStatistics::outputStatistics()
{
  if( !m_isIterativeSolver || !m_csvOutput )
    return;

  {
    GEOS_LOG_RANK_0( GEOS_FMT( "{}, number of Time-steps: {}", getParent().getName(), m_numTimeSteps ) );
    GEOS_LOG_RANK_0( GEOS_FMT( "{}, number of Time steps cut: {}", getParent().getName(), m_numTimeStepCuts ) );
    TableLayout iterationLogLayout ( GEOS_FMT( "{} iterations", getParent().getName() ),
                                     { "Components", "Iter  "} );

    TableTextFormatter const statsFormatter( iterationLogLayout );

    TableData iterationDataLog;
    iterationDataLog.addRow( "Successful outer loop", m_numSuccessfulOuterLoopIterations );
    iterationDataLog.addRow( "Successful nonlinear", m_numSuccessfulNonlinearIterations );
    iterationDataLog.addRow( "Successful linear", m_numSuccessfulLinearIterations );
    iterationDataLog.addRow( "Discarded outer loop", m_numDiscardedOuterLoopIterations );
    iterationDataLog.addRow( "Discarded nonlinear", m_numDiscardedNonlinearIterations );
    iterationDataLog.addRow( "Discarded linear", m_numDiscardedLinearIterations );

    GEOS_LOG_RANK_0( statsFormatter.toString( iterationDataLog ));
  }

  logStream.close();
}

ConvergenceStatistics::ConvergenceStatistics()
{
  m_convergenceLayout = std::make_unique< TableLayout >();

  m_convergenceLayout->addColumns( {"Cycle number", "time_n (s)", "dt (s)",
                                    "RMass", "RVol", "REnergy",
                                    "RFlow", "RBubbleDisp", "RFrac",
                                    "Rstick", "Rslip", "Ropen",
                                    "RSolid", "RContact", "RProppant",
                                    "RWell", "RDamage", "RTotal", "R"} );
}

void ConvergenceStatistics::writeResidualNormToTable()
{

  if( !m_csvOutput )
    return;

  stdVector< TableData::CellData > residualsNormCells;

  struct ResidualInfo
  {
    const real64 & value;
    string_view columnName;
  };

  ResidualInfo residuals[] = {
    { m_residualMass, "RMass" },
    { m_residualVol, "RVol" },
    { m_residualEnergy, "REnergy" },
    { m_residualFlow, "RFlow" },
    { m_residualBubbleDisp, "RBubbleDisp" },
    { m_residualFracture, "RFrac" },
    { m_residualStick, "Rstick" },
    { m_residualSlip, "Rslip" },
    { m_residualOpen, "Ropen" },
    { m_residualSolid, "RSolid" },
    { m_residualContact, "RContact" },
    { m_residualProppant, "RProppant" },
    { m_residualWell, "RWell" },
    { m_residualDamage, "RDamage" },
    { m_totalResidual, "RTotal" },
    { m_residualNormT, "R" }
  };

  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_cycleNumber )} ));
  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_time_n )} ));
  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_dt )} ));

  for( auto const & residual : residuals )
  {
    residualsNormCells.emplace_back( TableData::CellData(
      {
        CellType::Value,
        !std::isnan( residual.value ) ? GEOS_FMT( "{}", residual.value ) : "0"
      } ));
  }

  m_convergenceData.addRow( residualsNormCells );

  if( !logStream.is_open() )
  {
    logStream.open( m_convergenceFilename );
    m_convergenceFormatter.reset( new TableCSVFormatter( *m_convergenceLayout ));
    logStream << m_convergenceFormatter->headerToString( );
  }

  logStream << m_convergenceFormatter->dataToString( m_convergenceData );
  logStream.flush();
  m_convergenceData.clear();
}

void ConvergenceStatistics::updateSolverStep( real64 const & time_n, real64 const & dt, integer const cycleNumber )
{
  m_time_n = time_n;
  m_dt = dt;
  m_cycleNumber = cycleNumber;
}
void ConvergenceStatistics::resetResidualsValue()
{
  m_residualSolid = 0;
  m_totalResidual=0;
  m_residualNormT=0;
}

void ConvergenceStatistics::outputResidualNorm()
{
  if( m_csvOutput )
  {
    logStream << m_convergenceFormatter->dataToString( m_convergenceData );
  }
  logStream.close();

}

} // namespace geos
