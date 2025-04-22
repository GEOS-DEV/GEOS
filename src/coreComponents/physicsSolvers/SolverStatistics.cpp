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
#include "common/format/table/TableTypes.hpp"
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
  m_outputDir( joinPath( OutputBase::getOutputDirectory(), "convergence" ))// TODO DANS LE HPP
{
  makeDirsForPath( m_outputDir );

  using TableLayoutArgs = std::initializer_list< std::variant< string_view, TableLayout::Column > >;

  m_convergenceStats.m_nonLinearNormsLayout =  std::make_unique< TableLayout >(
    TableLayoutArgs{
      std::variant< string_view, TableLayout::Column >{"Time-steps"},
      std::variant< string_view, TableLayout::Column >{"Newton Iter"}
    } );
}

SolverStatistics::IterationsStatistics::IterationsStatistics()
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_iterationsStats.m_numTimeSteps ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time steps" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_iterationsStats.m_numTimeStepCuts ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of time step cuts" );

  registerWrapper( viewKeyStruct::numSuccessfulOuterLoopIterationsString(),
                   &m_iterationsStats.m_numSuccessfulOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful outer loop iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(),
                   &m_iterationsStats.m_numSuccessfulNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful nonlinear iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(),
                   &m_iterationsStats.m_numSuccessfulLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of successful linear iterations" );


  registerWrapper( viewKeyStruct::numDiscardedOuterLoopIterationsString(),
                   &m_iterationsStats.m_numDiscardedOuterLoopIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded outer loop iterations" );

  registerWrapper( viewKeyStruct::numDiscardedNonlinearIterationsString(),
                   &m_iterationsStats.m_numDiscardedNonlinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded nonlinear iterations" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(),
                   &m_iterationsStats.m_numDiscardedLinearIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Cumulative number of discarded linear iterations" );
}

void SolverStatistics::IterationsStatistics::initializeTimeStepStatistics()
{
  // the time step begins, we reset the individual-timestep counters
  m_currentNumOuterLoopIterations = 0;
  m_currentNumNonlinearIterations = 0;
  m_currentNumLinearIterations = 0;
}

void SolverStatistics::IterationsStatistics::logNonlinearIteration( integer const numLinearIterations )
{
  // we have just performed a Newton iteration, so we increment the individual-timestep counters
  m_currentNumNonlinearIterations++;
  m_currentNumLinearIterations += numLinearIterations;
}

void SolverStatistics::IterationsStatistics::logNonlinearIteration()
{
  // we have just performed an outer iteration, so we increment the individual-timestep counter (number of outer iteration)
  m_currentNumNonlinearIterations++;
}

void SolverStatistics::IterationsStatistics::logOuterLoopIteration()
{
  // we have just performed an outer loop iteration, so we increment the individual-timestep counter for outer loop iterations
  m_currentNumOuterLoopIterations++;
}

void SolverStatistics::IterationsStatistics::saveTimeStepStatistics()
{
  // the timestep has converged, so we increment the cumulative counters for successful timesteps
  m_numSuccessfulOuterLoopIterations += m_currentNumOuterLoopIterations;
  m_numSuccessfulNonlinearIterations += m_currentNumNonlinearIterations;
  m_numSuccessfulLinearIterations += m_currentNumLinearIterations;
  m_numTimeSteps++;
}

void SolverStatistics::IterationsStatistics::logTimeStepCut()
{
  // we have just cut the time step, so we increment the cumulative counters for discarded timesteps
  m_numDiscardedOuterLoopIterations += m_currentNumOuterLoopIterations;
  m_numDiscardedNonlinearIterations += m_currentNumNonlinearIterations;
  m_numDiscardedLinearIterations += m_currentNumLinearIterations;
  m_numTimeStepCuts++;

  // we are going to restart the timestep from the previous converged time step, so we have to re-initialize the statistics
  initializeTimeStepStatistics();
}

void SolverStatistics::IterationsStatistics::registerStatsToTable()
{
  m_nonLinearData.addRow( m_numTimeSteps,
                          m_numTimeStepCuts,
                          m_numSuccessfulOuterLoopIterations,
                          m_numSuccessfulNonlinearIterations,
                          m_numSuccessfulLinearIterations,
                          m_numDiscardedOuterLoopIterations,
                          m_numDiscardedNonlinearIterations,
                          m_numDiscardedLinearIterations );
}

void SolverStatistics::IterationsStatistics::outputStatistics( bool writeCSV )
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
    m_nonLinearData.clear();
  }
}

SolverStatistics::ConvergenceStatistics::ConvergenceStatistics()
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_currentNewtonIter )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum number of current Newton iterations" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_residualMass )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual mass" );

  registerWrapper( viewKeyStruct::numSuccessfulOuterLoopIterationsString(), &m_residualVol )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual volume" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(), &m_residualEnergy )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual energy" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(), &m_residualFlow )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual flow" );

  registerWrapper( viewKeyStruct::numDiscardedOuterLoopIterationsString(), &m_residualBubbleDisp )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual bubble displacement" );

  registerWrapper( viewKeyStruct::numDiscardedNonlinearIterationsString(), &m_residualFracture )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual fracture" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualStick )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual stick" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualSlip )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual slip" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualOpen )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual open" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualSolid )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual solid" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualContact )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual contact" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualProppant )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual proppant" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualWell )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual well" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_residualDamage )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum value for residual damage" );

  registerWrapper( viewKeyStruct::numDiscardedLinearIterationsString(), &m_totalResidual )
    .setApplyDefaultValue( std::numeric_limits< real64 >::max())
    .setDescription( "Maximum total residual value" );
}


void SolverStatistics::ConvergenceStatistics::removeInvalidResidualNorms()
{
  for( int i = 0; i <= m_currentNewtonIter; i++ )
    m_nonLinearNormsData.getTableDataRows().pop_back();
}

void SolverStatistics::ConvergenceStatistics::logNewtonIter( integer currentNewtonIter )
{ m_currentNewtonIter = currentNewtonIter; }



void SolverStatistics::ConvergenceStatistics::registerResidualNormToTable()
{
  std::vector< TableData::CellData > residualsNormCells;

  auto hasValue = []( real64 value )
  {
    return std::fabs( value - std::numeric_limits< real64 >::max()) > std::numeric_limits< real64 >::epsilon();
  };

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
    { m_totalResidual, "RTotal" }
  };


  if( m_numTimeSteps == 0 && m_currentNewtonIter == 0 )
  {
    for( auto const & residual : residuals )
    {
      if( hasValue( residual.value ) )
      {
        m_nonLinearNormsLayout->addColumn( residual.columnName );
      }
    }
  }

  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_numTimeSteps )} ));
  residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                         GEOS_FMT( "{}", m_currentNewtonIter )} ));

  for( auto const & residual : residuals )
  {
    if( hasValue( residual.value ) )
    {
      residualsNormCells.emplace_back( TableData::CellData( {CellType::Value,
                                                             GEOS_FMT( "{}", residual.value )} ));
    }
  }

  m_nonLinearNormsData.addRow( residualsNormCells );
}

void SolverStatistics::ConvergenceStatistics::outputResidualNorm( bool writeCSV )
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
