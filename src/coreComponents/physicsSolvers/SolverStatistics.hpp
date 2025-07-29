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

  /// indicate if the containing solver does non-linear iterations (and so, produces iterations statistics)
  bool m_isIterativeSolver = true;

  /// State of log output. True when writeSolverStatistics is set to 1
  bool m_logOutput = false;

  /// State of csv output. True when writeSolverStatistics is set to 2
  bool m_csvOutput = false;

  /// Number of time steps
  integer m_numTimeSteps = 0;

  /// Number of time step cuts
  integer m_numTimeStepCuts = 0;

  /// Linear solver setup
  real64 m_setupTime = 0.0;

  /// Linear solver solve
  real64 m_solveTime = 0.0;

  /// Maximum number of current Newton iterations.
  integer m_currentNewtonIter = 0;


  /// Number of outer loop iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumOuterLoopIterations = 0;

  /// Number of nonlinear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumNonlinearIterations = 0;

  /// Number of linear iterations in the current time step (utility variable constantly overwritten)
  integer m_currentNumLinearIterations = 0;


  /// Cumulative number of successful outer loop iterations
  integer m_numSuccessfulOuterLoopIterations = 0;

  /// Cumulative number of successful nonlinear iterations
  integer m_numSuccessfulNonlinearIterations = 0;

  /// Cumulative number of successful linear iterations
  integer m_numSuccessfulLinearIterations = 0;


  /// Cumulative number of discarded outer loop iterations
  integer m_numDiscardedOuterLoopIterations = 0;

  /// Cumulative number of discarded nonlinear iterations
  integer m_numDiscardedNonlinearIterations = 0;

  /// Cumulative number of discarded linear iterations
  integer m_numDiscardedLinearIterations = 0;

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

    /// String key for the successful number of outer loop iterations
    static constexpr char const * numSuccessfulOuterLoopIterationsString() { return "numSuccessfulOuterLoopIterations"; }
    /// String key for the successful number of nonlinear iterations
    static constexpr char const * numSuccessfulNonlinearIterationsString() { return "numSuccessfulNonlinearIterations"; }
    /// String key for the successful number of linear iterations
    static constexpr char const * numSuccessfulLinearIterationsString() { return "numSuccessfulLinearIterations"; }

    /// String key for the discarded number of outer loop iterations
    static constexpr char const * numDiscardedOuterLoopIterationsString() { return "numDiscardedOuterLoopIterations"; }
    /// String key for the discarded number of nonlinear iterations
    static constexpr char const * numDiscardedNonlinearIterationsString() { return "numDiscardedNonlinearIterations"; }
    /// String key for the discarded number of linear iterations
    static constexpr char const * numDiscardedLinearIterationsString() { return "numDiscardedLinearIterations"; }
  };

  /**
   * @brief Set whether the solver is iterative.
   * @param isIterative The iterative state.
   */
  void setIterativeSolver( bool isIterative )
  {
    m_isIterativeSolver = isIterative;
  }

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
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @param[in] numLinearIterations the number of linear iterations done by the linear solver
   * @detail This function is well suited for Newton's method, or for single-physics solvers in sequential schemes
   */
  void updateNonlinearIteration( integer const numLinearIterations );

  /**
   * @brief Save the current newton iteration
   * @param currentNewtonIter The current newton iteration performed by the the linear solver
   */
  void updateNewtonIter( integer currentNewtonIter )
  { if( m_isIterativeSolver ) m_currentNewtonIter = currentNewtonIter; }

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
   * @brief Tell the solverStatistics that we are doing a nonlinear iteration
   * @detail This function is well suited for the outer loop in sequential schemes
   */
  void updateNonlinearIteration();

  /**
   * @brief Tell the solverStatistics that we are doing an outer loop iteration
   */
  void incrementNonlinearIteration();

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
  { if( m_isIterativeSolver ) m_iterationsFilename = filename; }

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

  /// State of csv output. True when writeSolverStatistics is set to 2
  bool m_csvOutput = false;

  /// The time at the beginning of the step
  real64 m_time_n = 0.0;

  /// The desired timestepe
  real64 m_dt = 0.0;

  /// Current cycle number
  integer m_cycleNumber = 0;

  /// Maximum value for residual mass.
  real64 m_residualMass = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual volume.
  real64 m_residualVol = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual energy.
  real64 m_residualEnergy = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual flow.
  real64 m_residualFlow = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual bubble displacement.
  real64 m_residualBubbleDisp = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual fracture.
  real64 m_residualFracture = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual stick.
  real64 m_residualStick = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual slip.
  real64 m_residualSlip = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual open.
  real64 m_residualOpen = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual solid.
  real64 m_residualSolid = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual contact.
  real64 m_residualContact = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual proppant.
  real64 m_residualProppant = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual damage.
  real64 m_residualDamage = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum total residual value.
  real64 m_totalResidual = std::numeric_limits< real64 >::quiet_NaN();

  /// Maximum value for residual norm at the end of line search
  real64 m_residualNormT = std::numeric_limits< real64 >::quiet_NaN();

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
   */
  void updateSolverStep( real64 const & time_n, real64 const & dt, integer const cycleNumber );

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
   * @noteCurrently we register only the iteration statistics as the convergence value will be lost during the solving
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

  /**
   * @return Return the string directory where all CSV are generated
   */
  string_view getDirectory() const
  { return m_directoryName;}

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
