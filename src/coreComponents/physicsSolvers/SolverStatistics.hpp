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
 * @file SolverStatistics.hpp
 */


#ifndef GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
#define GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP

#include "dataRepository/Group.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/format/table/TableTypes.hpp"

namespace geos
{

/**
 * @brief Class containing solver iterations data for a time-step
 */
class IterationsStatistics : public dataRepository::Group
{
public:

/**
 * @brief Constructor for SolverStatistics Objects.
 * @param[in] name the name of this instantiation of SolverStatistics in the repository
 * @param[in] parent the parent group of this instantiation of SolverStatistics
 * @note Currently, we register only the iteration statistics as the convergence value will be lost during the solving
 */
  IterationsStatistics( string const & name,
                        dataRepository::Group * const parent );

  /// State of log output. True when writeStatistics is set to 1
  bool m_logOutput = false;

  /// State of csv output. True when writeStatistics is set to 2
  bool m_csvOutput = false;

  /// Number of time steps
  integer m_numTimeSteps = 0;

  /// Number of time step cuts
  integer m_numTimeStepCuts = 0;

  /// Number of configuration iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumConfigIterations = 0;

  /// Number of nonlinear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumNonlinearIterations = 0;

  /// Number of linear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumLinearIterations = 0;

  /// Cumulative number of successful configuration iterations
  integer m_numSuccessfulConfigIterations = 0;

  /// Cumulative number of successful nonlinear iterations
  integer m_numSuccessfulNonlinearIterations = 0;

  /// Cumulative number of successful linear iterations
  integer m_numSuccessfulLinearIterations = 0;

  /// Cumulative number of discarded configuration iterations
  integer m_numDiscardedConfigIterations  = 0;

  /// Cumulative number of discarded nonlinear iterations
  integer m_numDiscardedNonlinearIterations = 0;

  /// Cumulative number of discarded linear iterations
  integer m_numDiscardedLinearIterations = 0;

  /// Linear solver setup
  real64 m_setupTime = 0.0;

  /// Linear solver solve
  real64 m_solveTime = 0.0;

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the number of time steps
    static constexpr char const * numTimeStepsString() { return "numTimeSteps"; }
    /// String key for the number of time step cuts
    static constexpr char const * numTimeStepCutsString() { return "numTimeStepCuts"; }

    /// String key for the successful number of configuration iterations
    static constexpr char const * numSuccessfulConfigIterationsString() { return "numSuccessfulConfigIterations"; }
    /// String key for the successful number of nonlinear iterations
    static constexpr char const * numSuccessfulNonlinearIterationsString() { return "numSuccessfulNonlinearIterations"; }
    /// String key for the successful number of linear iterations
    static constexpr char const * numSuccessfulLinearIterationsString() { return "numSuccessfulLinearIterations"; }

    /// String key for the discarded number of configuration iterations
    static constexpr char const * numDiscardedConfigIterationsString() { return "numDiscardedConfigIterations"; }
    /// String key for the discarded number of nonlinear iterations
    static constexpr char const * numDiscardedNonlinearIterationsString() { return "numDiscardedNonlinearIterations"; }
    /// String key for the discarded number of linear iterations
    static constexpr char const * numDiscardedLinearIterationsString() { return "numDiscardedLinearIterations"; }
  };

  /**
   * @brief Set the csv state output
   * @param state The csv state
   */
  void setCSVOutput( bool state )
  { m_csvOutput = state; }

  /**
   * @brief Set the log state output
   * @param state The log state
   */
  void setLogOutput( bool state )
  { m_logOutput = state; }


  /**
   * @brief Initialize the counters used for an individual time step
   */
  void resetCurrentTimeStepStatistics();

  /**
   * @brief Tell the solverStatistics that we are doing a configuration iteration
   */
  void incrementConfigIteration();

  /**
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @param[in] numLinearIterations the number of linear iterations done by the linear solver
   * @detail This function is well suited for Newton's method, or for single-physics solvers in sequential schemes
   */
  void updateNonlinearIteration( integer const numLinearIterations );

  /**
   * @brief Accumulate the setupTime & solveTime result over each newton iteration
   * @param setupTime The linear setup time
   * @param solveTime The linear solve time
   */
  void accumulateSolverLinearTime( real64 setupTime, real64 solveTime );

  /**
   * @brief Reset the setupTime & solveTime to 0 at the end of each cycle
   */
  void resetSolverLinearTime();

  /**
   * @brief Tell the solverStatistics that there is a time step cut
   */
  void updateTimeStepCut();

  /**
   * @brief Save the statistics for the individual time step and increment the cumulative stats
   */
  void iterateTimeStepStatistics();

  /**
   * @brief Write all the iteration statistics into the ouput stream
   */
  void writeIterationStatsToTable();

  /**
   * @brief  Set the filename output file.
   * @param filename The filename as a string_view.
   */
  void setFilename( string_view filename )
  { m_iterationsFilename = filename; }

  /**
   * @brief Set the filename output file.
   * @param name The filename as a string_view.
   */
  void setTableName( string_view name )
  { m_tableIterationName = name; }

  /**
   * @brief Output the statistics to the console in table format
   */
  void outputStatistics();

private:
  /// Stream output for the iteration statistics
  std::ofstream m_logStream;
  /// Table Layout contenaning header for both CSV and log
  std::unique_ptr< TableLayout > m_iterationCSVLayout;
  /// Contain the iteration data for both CSV and log output
  TableData m_iterationData;
  /// Format the iteration statistics for the CSV file
  std::unique_ptr< TableCSVFormatter > m_iterationCSVFormatter;
  /// String name for the log table interation
  string m_tableIterationName;
  /// Filename for the iteration CSV file.
  string m_iterationsFilename;
};

/**
 * @brief Class containing convergence information given a time-step
 */
class ConvergenceStatistics
{
public:

  /**
   * @brief Construct a new Convergence Statistics object
   */
  ConvergenceStatistics();

  /// State of csv output. True when writeSolver is set to 2
  bool m_csvOutput = false;

  /// The time at the beginning of the step
  real64 m_time_n = 0.0;

  /// The desired timestep
  real64 m_dt = 0.0;

  /// Current cycle number
  integer m_cycleNumber = 0;

  /// Current newton iteration
  integer m_iteration = 0;

  /// Residuals with their names
  std::map< string, real64 > m_residuals;

  /**
   * @brief Set the csv state output
   * @param state csv state
   */
  void setCSVOutput( bool state )
  { m_csvOutput = state; }

  /**
   * @brief Write all the convergence statistics into the ouput stream
   */
  void writeConvergenceStatsToTable();

  /**
   * @brief Update the solver step with the time informations
   * @param time_n The time at the beginning of the step
   * @param dt The desired timestep
   * @param cycleNumber The current cycle number
   * @param iterNumber The current iteration number
   */
  void updateSolverStep( real64 const & time_n,
                         real64 const & dt,
                         integer const cycleNumber,
                         integer const newtonIter );

  /**
   * @brief Reset the solid residuals value.
   * Call by SolidMechanicsStateReset.
   */
  void resetResidualsValue();

  /**
   * @brief Set the filename output file.
   * @param filename The filename as a string_view.
   */
  void setFilename( string_view filename )
  { m_convergenceFilename = filename; }

  /**
   * @brief Set the filename output file.
   * @param filename The filename as a string_view.
   */
  void setTableName( string_view name )
  { m_tableConvergenceName = name; }

private:
  /// Stream output for the convergence statistics
  std::ofstream m_logStream;
  /// Contain the layout for both the CSV and log output.
  /// For a solver, output all residuals residuals name available
  std::unique_ptr< TableLayout > m_convergenceLayout;
  /// Contain the convergence data for both CSV and log output
  TableData m_convergenceData;
  /// Format the convergence statistics for the CSV file
  std::unique_ptr< TableCSVFormatter > m_convergenceFormatter;
  /// String name for the log table convergence
  string m_tableConvergenceName;
  /// Filename for the solver CSV convergence file.
  string m_convergenceFilename;
};

/**
 * @class SolverStatistics
 * @brief This class records solver statistics for each time step.
 */
class SolverStatistics : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for SolverStatistics Objects.
   * @param[in] name the name of this instantiation of SolverStatistics in the repository
   * @param[in] parent the parent group of this instantiation of SolverStatistics
   * @note Currently we register only the iteration statistics as the convergence value will be lost during the solving
   */
  SolverStatistics( string const & name,
                    dataRepository::Group * const parent );

  /**
   * @brief Group key associated with IterationsStatistics.
   */
  struct groupKeyStruct
  {
    /// @return string for the IterationsStatistics wrapper
    static constexpr char const * IterationsStatisticsString() { return "IterationsStatistics"; }
  };

  /**
   * @brief Create a convergence directory if we enable the csv for iteration or convergence statistics
   * @param writeConvergence Boolean for convergencce CSV output
   * @param writeIteration Boolean for iteration CSV output
   */
  void makeDir( bool writeSolverIteration )
  { if( writeSolverIteration ) makeDirsForPath( m_outputDir ); }

  /**
   * @brief Set the Residual Norms filename
   * @param solverName The solverName as a string_view.
   */
  void setOutputFilesName( string_view solverName );

  /// Contain iteration data given a time step
  IterationsStatistics m_iterationsStats;
  /// Contain convergence data given a time step
  ConvergenceStatistics m_convergenceStats;

protected:
  /// Output directory for solver statistics (CSV), passed in the constructor.
  string m_outputDir;

private:
  /// Name of the directory containing solvers statistics csv
  string m_directoryName = "convergence";
};

} //namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLVERSTATISTICS_HPP
